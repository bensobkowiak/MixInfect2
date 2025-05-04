##########################################################################################################################
###################### MixInfect2 - Identify mixed samples from VCF file #################################################
##########################################################################################################################

if (!require("mclust")){install.packages("mclust")}
if (!require("stringr")){install.packages("stringr")}
if (!require("optparse")){install.packages("optparse")}
if (!require("foreach")){install.packages("foreach")}
if (!require("doMC")){install.packages("doMC")}

library(mclust)
library(stringr)
library(optparse)
library(foreach)
library(doMC)

MixInfect2 <- function(VCFfile, prefix = "output", maskFile = NULL, useFilter = TRUE, 
                       minQual = 20, LowCov = 10, minDepth = 5, 
                       popFreq_threshold = 0.1, SNPwindow = 100, n_threads = 4) {
  
  # Register the number of threads for parallel processing
  registerDoMC(cores = n_threads)
  
  # Read VCF file
  vcf <- read.table(VCFfile, stringsAsFactors = FALSE)
  header_input <- as.matrix(read.table(VCFfile, comment.char=" ", sep="\n", stringsAsFactors = FALSE))
  end_head <- which(grepl("#CHROM", header_input) == TRUE)
  header <- as.data.frame(header_input[1:end_head-1])
  names <- unlist(strsplit(header_input[end_head], "\t"))
  colnames(vcf) <- names
  format <- which(names == "FORMAT")
  head_start <- names[1:format]
  names <- names[(format + 1):length(names)]
  rm(header_input)
  
  ## Remove indels
  ind <- lapply(1:nrow(vcf), function(i) {
    length(unlist(strsplit(vcf[i, 4], ""))) > 1 || length(unlist(strsplit(unlist(strsplit(vcf[i, 5], ","))[1], ""))) > 1
  })
  vcf <- vcf[which(ind==FALSE), ]
  
  # Remove filtered variants
  if (useFilter) {
    vcf <- vcf[vcf[, which(head_start == "FILTER")] == "PASS", ]
  }
  
  # Remove low quality variants
  vcf <- vcf[vcf[, which(head_start == "QUAL")] > minQual, ]
  
  # Remove spanning variants
  vcf <- vcf[grep("*", vcf[, which(head_start == "ALT")], fixed = TRUE, invert = TRUE), ]
  
  # Remove SNPs in repeat regions
  if (!is.null(maskFile)) {
    maskFile <- read.csv(maskFile)
    res1 <- vector()
    for (i in 1:nrow(maskFile)) {
      res1 <- c(res1, c(maskFile[i, 1]:maskFile[i, 2]))
    }
    vcf <- vcf[which(!vcf[, 2] %in% res1), ]
  }
  
  # Determine AD field and create matrix of AD and GT
  AD <- which(unlist(str_split(vcf[1, format], ":")) == 'AD')
  GT <- which(unlist(str_split(vcf[1, format], ":")) == 'GT')
  DP <- which(unlist(str_split(vcf[1, format], ":")) == 'DP')
  
  # Make new matrices of separated GT, DP, and AD fields
  GT_mat <- matrix(ncol = length((format + 1):ncol(vcf)), nrow = nrow(vcf))
  AD_mat <- matrix(ncol = length((format + 1):ncol(vcf)), nrow = nrow(vcf))
  DP_mat <- matrix(ncol = length((format + 1):ncol(vcf)), nrow = nrow(vcf))
  for (i in 1:ncol(GT_mat)) {
    GT_mat[, i] <- sapply(str_split(vcf[, i + format], ":"), "[[", GT)
    newDP <- str_split(vcf[, i + format], ":")
    newDP <- lapply(newDP, function(x) if (length(x) < 3) c(x, 0,0,0) else x)
    DP_mat[, i] <- sapply(newDP, "[[", DP)
    AD_mat[, i] <- sapply(newDP, "[[", AD)
  }

  DP_mat[which(DP_mat == ".")] <- "0"
  GT_mat[which(as.numeric(DP_mat) < LowCov)] <- "?"
  
  # Create output file
  outfile <- as.data.frame(matrix(NA, nrow = length(names), ncol = 7))
  colnames(outfile) <- c("SampleName", "Mix.Non-mix", "hSNPs", "Total.SNPs", 
                         "Proportion.hSNPs_totalSNPs", "No.strains", "Major.strain.proportion")
  outfile[, 1] <- names
  
  # Identify mixed calls
  mixed_calls <- c("0/1", "0/2", "0/3", "1/2", "1/3", "2/3",
                   "0|1", "0|2", "0|3", "1|2", "1|3", "2|3")
  alt_calls <- c("1/1", "2/2", "3/3", "1|1", "2|2", "3|3")
  
  # Process AD fields to identify mixed calls
  for (col in 1:ncol(AD_mat)) {
    ADmix <- str_split(AD_mat[, col], ",")
    for (m in 1:length(ADmix)) {
      AD_site <- as.numeric(unlist(ADmix[m]))
      AD_site <- AD_site[AD_site != 0]
      if (length(AD_site) > 1) {
        AD_site <- AD_site[order(AD_site, decreasing = TRUE)]
        if (AD_site[2] >= minDepth) {
          GT_mat[m, col] <- "0/1"
        }
      }
    }
  }
  
  # Keep loci with an alternative or mixed call
  keep <- apply(as.data.frame(GT_mat), 1, function(row) {
    any(row %in% c(mixed_calls, alt_calls))
  })
  GT_mat <- GT_mat[keep, , drop=F]
  vcf <- vcf[keep, ,drop = F]
  
  # Mask sites with mixed frequency over popFreq_threshold
  propMix <- numeric()
  for (i in 1:nrow(GT_mat)) {
    propMix[i] <- length(which(GT_mat[i, ] %in% mixed_calls)) / ncol(GT_mat)
  }
  propMix <- which(propMix > popFreq_threshold)
  if (length(propMix) > 0) {
    GT_mat <- GT_mat[-propMix, ]
    vcf <- vcf[-propMix, ]
  }
  
  # Calculate hSNPs, total SNPs, and proportions
  mixes <- matrix(0, ncol = ncol(GT_mat), nrow = 4)
  for (i in 1:ncol(GT_mat)) {
    mixes[1, i] <- length(which(GT_mat[, i] %in% mixed_calls))
  }
  for (i in 1:ncol(GT_mat)) {
    mixes[2, i] <- length(which(GT_mat[, i] %in% alt_calls))
  }
  for (i in 1:ncol(GT_mat)) {
    mixes[3, i] <- mixes[1, i] + mixes[2, i]
  }
  for (i in 1:ncol(GT_mat)) {
    mixes[4, i] <- (mixes[1, i] / mixes[3, i]) * 100
  }
  
  outfile[, 3] <- mixes[1, ]
  outfile[, 4] <- mixes[3, ]
  outfile[, 5] <- round(mixes[4, ],2)
  outfile[, 2] <- 'Non-mix'
  outfile[, 6] <- 1
  
  #################### ESTIMATE PROPORTIONS OF MIXED SAMPLES #######################
  
  mixnames <- outfile$SampleName[which(outfile[, 5] > 1.5 & outfile[, 3] > 10)]
  mix_GT <- as.data.frame(GT_mat[, which(outfile[, 5] > 1.5 & outfile[, 3] > 10)])
  mix_VCF <- as.data.frame(vcf[, which(outfile[, 5] > 1.5 & outfile[, 3] > 10) + format])
  positions <- vcf[, 2]
  
  if (length(mixnames) > 0) {
    BICvalues <- data.frame(Sample = mixnames, G2 = NA_real_, G4 = NA_real_, G6 = NA_real_)
    
    results <- foreach(i=1:nrow(BICvalues)) %dopar% {
      samplemix_sites <- mix_VCF[which(mix_GT[, i] %in% mixed_calls), i]
      samplemix_AD <- sapply(str_split(samplemix_sites, ":"), "[[", AD)
      samplePos <- positions[which(mix_GT[, i] %in% mixed_calls)]
      samplemaj_prop <- numeric()
      samplemin_prop <- numeric()
      finalPos <- numeric()
      for (k in 1:length(samplemix_AD)) {
        sampleAD <- as.numeric(unlist(str_split(samplemix_AD[k], ",")))
        sampleAD <- sampleAD[sampleAD != 0]
        sampleAD <- sampleAD[order(sampleAD, decreasing = TRUE)]
        if (length(sampleAD) == 2) {
          samplemaj_prop <- c(samplemaj_prop, sampleAD[1] / sum(sampleAD))
          samplemin_prop <- c(samplemin_prop, sampleAD[2] / sum(sampleAD))
          finalPos <- c(finalPos, samplePos[k])
        }
      }
      
      distances <- diff(finalPos)
      group_indices <- c(1, cumsum(distances >= SNPwindow) + 1)
      samplemaj_prop <- tapply(samplemaj_prop, group_indices, median)
      samplemin_prop <- tapply(samplemin_prop, group_indices, median)
      b <- c(samplemaj_prop, samplemin_prop)
      
      bic_values <- c(NA_real_, NA_real_, NA_real_)
      mix_status <- 'Non-mix'
      no_strains <- 1
      major_strain_proportion <- NA_real_
      
      if (length(b) > LowCov) {
        a <- mclustBIC(b, G = c(2, 4, 6), verbose = FALSE)[, 2]
        if (length(a) == 3) {
          bic_values <- a
        }
        
        d <- Mclust(b, G = 2, verbose = FALSE)
        if (!is.na(bic_values[1]) && bic_values[1] >= 20) {
          mix_status <- 'Mix'
          no_strains <- 2
          major_strain_proportion <- d$parameters$mean[order(d$parameters$mean, decreasing = TRUE)][1]
        } else if (!is.na(bic_values[3]) && bic_values[3] >= 20) {
          d <- Mclust(b, G = 6, verbose = FALSE)
          mix_status <- 'Mix'
          no_strains <- 3
          means <- d$parameters$mean[order(d$parameters$mean, decreasing = TRUE)]
          if (sum(means[5:6]) < 0.5) {
            major_strain_proportion <- means[3]
          } else {
            major_strain_proportion <- means[4]
          }
        }
      }
      print(paste("Processed sample", mixnames[i]))
      return(c(i, bic_values, mix_status, no_strains, major_strain_proportion))
    }
    
    for (res in 1:length(results)) {
      i <- as.numeric(unlist(results[res])[1])
      BICvalues[i, 2:4] <- as.numeric(unlist(results[res])[2:4])
      ind <- which(BICvalues$Sample[i] == outfile$SampleName)
      outfile[ind, 2] <- unlist(results[res])[5]
      outfile[ind, 6] <- as.numeric(unlist(results[res])[6])
      outfile[ind, 7] <- round(as.numeric(unlist(results[res])[7]),2)
    }
    
    BICvalues <- BICvalues[which(BICvalues$Sample %in% outfile$SampleName[which(outfile$`Mix.Non-mix` == "Mix")]), ]
    write.csv(BICvalues, paste0(prefix, "_BICvalues.csv"), row.names = FALSE)
    write.csv(outfile, paste0(prefix, "_MixSampleSummary.csv"), row.names = FALSE)
  } else {
    print("No mixed infection")
  }
}

