##########################################################################################################################
###################### MixInfect2 - Identify mixed samples from VCF file #################################################
##########################################################################################################################

#' @param VCFfile VCF file - must have GT and AD fields
#' @param prefix Output file prefix
#' @param maskFile CSV file with regions to mask, with start position in column 1 and end position in column 2
#' @param useFilter Use the 'FILTER' column in VCF file to filter SNPs
#' @param minQual Minimum per loci quality
#' @param LowCov Minimum read depth at site to call either a cSNP or hSNP allele frequency
#' @param popFreq_threshold Maximum threshold for hSNP to be present in all samples (set as 0 to keep all sites)
#' @param SNPwindow Take the median of hSNP allele frequencies within this distance on the genome
#' @return Two CSV files, one with summary of mixed samples, one with BIC values for mixed samples
#' @export

MixInfect2<-function(VCFfile, prefix = "output", maskFile = NULL, useFilter = TRUE, 
                     minQual = 20, LowCov = 10, popFreq_threshold = 0.1, SNPwindow = 100){
  
  if (!require("mclust")){install.packages("mclust")}
  if (!require("stringr")){install.packages("stringr")}
  library(mclust)
  library(stringr)

  options(stringsAsFactors = F)
  
  # Read VCF file
  vcf<-read.table(VCFfile)
  header_input<-as.matrix(read.table(VCFfile,comment.char=" ",sep="\n"))
  end_head<-which(grepl("#CHROM",header_input)==TRUE)
  header<-as.data.frame(header_input[1:end_head-1])
  names<-unlist(strsplit(header_input[end_head],"\t"))
  colnames(vcf)<-names
  format<-which(names=="FORMAT")
  head_start<-names[1:format]
  names<-names[(format+1):length(names)]
  rm(header_input)
  
  # Remove filtered variants
  if (useFilter){
    vcf<-vcf[vcf[,which(head_start=="FILTER")]=="PASS",]
  }
  
  # Remove low quality variants
  vcf<-vcf[vcf[,which(head_start=="QUAL")]>minQual,]
  
  # Remove spanning variants
  vcf<-vcf[grep("*",vcf[,which(head_start=="ALT")],fixed = T,invert = T),]
  
  ## Remove SNPs in repeat regions
  if (!is.null(maskFile)){
    maskFile<-read.csv(maskFile)
    res1=vector()
    for (i in 1:nrow(maskFile)){
      res1<-c(res1,c(maskFile[i,1]:maskFile[i,2]))
    }
    vcf<-vcf[which(!vcf[,2] %in% res1),]
  }
  
  ### alternative calls
  mixed_calls<-c("0/1","0/2","0/3","1/2","1/3","2/3",
                 "0|1","0|2","0|3","1|2","1|3","2|3")
  alt_calls<-c("1/1","2/2","3/3","1|1","2|2","3|3")
  
  #### Determine AD field and create matrix of AD and GT
  AD<-which(unlist(str_split(vcf[1,format], ":"))=='AD')
  GT<-which(unlist(str_split(vcf[1,format], ":"))=='GT')
  DP<-which(unlist(str_split(vcf[1,format], ":"))=='DP')
  
  #### Make new matrices of separated GT, DP, and AD fields

  GT_mat<-matrix(ncol = length((format+1):ncol(vcf)),nrow = nrow(vcf))
  AD_mat<-matrix(ncol = length((format+1):ncol(vcf)),nrow = nrow(vcf))
  DP_mat<-matrix(ncol = length((format+1):ncol(vcf)),nrow = nrow(vcf))
  for (i in 1:ncol(GT_mat)) {
    GT_mat[,i] <- sapply(str_split(vcf[, i+format], ":"), "[[", GT)
    newDP <- str_split(vcf[, i+format], ":")
    newDP <- lapply(newDP, function(x) if (length(x) < 3) c(x, 0) else x)
    DP_mat[,i] <- sapply(newDP, "[[", DP)
    AD_mat[,i] <- sapply(str_split(vcf[, i+format], ":"), "[[", AD)
  }
  
  DP_mat[which(DP_mat==".")]<-"0"
  GT_mat[which(as.numeric(DP_mat) < LowCov)]<-"?"
  
  ### Create output file
  outfile<-as.data.frame(matrix(NA,nrow=length(names),ncol=7))
  colnames(outfile)<-c("SampleName","Mix.Non-mix","hSNPs","Total.SNPs",
                       "Proportion.hSNPs_totalSNPs","No.strains","Major.strain.proportion")
  outfile[,1]<-names
  
  ## Which sites have mixed reads (not just 0/1 etc.)
  for (col in 1:ncol(AD_mat)){
    ADmix<-str_split(AD_mat[,col], ",")
    for (m in 1:length(ADmix)){
      AD_site<-as.numeric(unlist(ADmix[m],","))
      AD_site<-AD_site[!AD_site==0]
      if (length(AD_site)>1){
        AD_site<-AD_site[order(AD_site,decreasing = T)]
        if (AD_site[2]>=LowCov){
          GT_mat[m,col]<-"0/1"
        }
      }
    }
  }
  
  #### Mask  sites with mixed frequency over popFreq_threshold
  propMix<-numeric()
  for (i in 1:nrow(GT_mat)){
    propMix[i]<-length(which(GT_mat[i,] %in% mixed_calls))/ncol(GT_mat)
  }
  propMix<-which(propMix>popFreq_threshold)
  if (length(propMix) > 0){
    GT_mat<-GT_mat[-propMix,]
    vcf<-vcf[-propMix,]
  }
  
  ## Keep loci with an alternative or mixed call
  keep<-apply(as.data.frame(GT_mat), 1, function(row) {
    any(row %in% c(mixed_calls,alt_calls))
  })
  GT_mat<-GT_mat[keep,]
  vcf<-vcf[keep,]
  
  #### No of het SNPs and total and proportions
  mixes<-matrix(0,ncol=ncol(GT_mat),nrow=4)
  for (i in 1:ncol(GT_mat)){
    mixes[1,i]<-length(which(GT_mat[,i] %in% mixed_calls))
  }
  for (i in 1:ncol(GT_mat)){
    mixes[2,i]<-length(which(GT_mat[,i] %in% alt_calls))
  }
  for (i in 1:ncol(GT_mat)){
    mixes[3,i]<-mixes[1,i]+mixes[2,i]
  }
  for (i in 1:ncol(GT_mat)){
    mixes[4,i]<-(mixes[1,i]/mixes[3,i])*100
  }
  
  outfile[,3]<-mixes[1,]
  outfile[,4]<-mixes[3,]
  outfile[,5]<-mixes[4,]
  outfile[,2]<-'Non-mix'
  outfile[,6]<-1
  
  #################### ESTIMATE PROPORTIONS OF MIXED SAMPLES #######################
  
  mixnames<-outfile$SampleName[which(outfile[,5]>1.5 & outfile[,3]>10)]
  mix_GT<-as.data.frame(GT_mat[,which(outfile[,5]>1.5 & outfile[,3]>10)])
  mix_VCF<-as.data.frame(vcf[,which(outfile[,5]>1.5 & outfile[,3]>10)+format])
  positions<-vcf[,2]
  
  if (length(mixnames)>0){
    ####### Run Guassian Mclust and idenitify 2 or 3 mixes
    BICvalues<-data.frame(Sample=mixnames,G2=0,G4=0,G6=0)
    
    for (i in 1:nrow(BICvalues)){
      # find mixed sites
      samplemix_sites<-mix_VCF[which(mix_GT[,i] %in% mixed_calls),i]
      samplemix_AD<-sapply(str_split(samplemix_sites, ":"), "[[", AD)
      samplePos<-positions[which(mix_GT[,i] %in% mixed_calls)]
      samplemaj_prop<-numeric()
      samplemin_prop<-numeric()
      finalPos<-numeric()
      for (k in 1:length(samplemix_AD)){
        sampleAD<-as.numeric(unlist(str_split(samplemix_AD[k], ",")))
        sampleAD<-sampleAD[sampleAD!=0]
        sampleAD<-sampleAD[order(sampleAD,decreasing = T)] 
        if (length(sampleAD)==2){
          samplemaj_prop<-c(samplemaj_prop,sampleAD[1]/sum(sampleAD))
          samplemin_prop<-c(samplemin_prop,sampleAD[2]/sum(sampleAD))
          finalPos<-c(finalPos,samplePos[k])
        } 
      }
      
      # Find SNPs within SNPwindow of each other and take median
      distances <- diff(finalPos)
      group_indices <- c(1, cumsum(distances >= SNPwindow) + 1)
      samplemaj_prop <- tapply(samplemaj_prop, group_indices, median)
      samplemin_prop <- tapply(samplemin_prop, group_indices, median)
      b<-c(samplemin_prop,samplemaj_prop)
      
      ## Calculate BICs
      if (length(b)>LowCov){
        a<-mclustBIC(b,G=c(2,4,6),verbose = F)[,2]
        if (length(a)==3){
          BICvalues[i,2:4]<-a
        } else { BICvalues[i,2:4]<-c(a,rep(NA,3-length(a)))
        }
      }
      
      ind<-which(BICvalues$Sample[i]==outfile$SampleName)
      d<-Mclust(b,G=2,verbose = F)
      if (BICvalues[i,2]>=20 && is.na(BICvalues[i,2])==FALSE){
        outfile[ind,2]<-'Mix'
        outfile[ind,6]<-2
        outfile[ind,7]<-d$parameters$mean[order(d$parameters$mean,decreasing = T)][1]
      } else if (BICvalues[i,4]>=20 && is.na(BICvalues[i,4])==FALSE){
        d<-Mclust(b,G=6,verbose = F)
        outfile[ind,2]<-'Mix'
        outfile[ind,6]<-3
        means<-d$parameters$mean[order(d$parameters$mean,decreasing = T)]
        if (sum(means[5:6])<0.5){
          outfile[ind,7]<-d$parameters$mean[order(d$parameters$mean,decreasing = T)][3]
        } else {
          outfile[ind,7]<-d$parameters$mean[order(d$parameters$mean,decreasing = T)][4]
        }
      }
    }
    
    # Output files
    BICvalues<-BICvalues[which(BICvalues$Sample %in% 
                                 outfile$SampleName[which(outfile$`Mix.Non-mix`=="Mix")]),]
    write.csv(BICvalues,paste0(prefix,"_BICvalues.csv"),row.names = F)
    write.csv(outfile,paste0(prefix,"_MixSampleSummary.csv"),row.names = F)
  } else {
    print("No mixed infection")
  }
}
                    
