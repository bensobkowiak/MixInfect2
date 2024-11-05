findRelatedStrains<-function(non_mixed_strainDF, nonMixDistances, 
                             ADlist, nuclist,maxDistance){
  df<-non_mixed_strainDF # dataframe of mixed sites from all pure strains
  AD<-ADlist # list of read depths at mixed sites
  Nuc<-nuclist # list of nucleotides at mixed sites
  nonmix_dist<-nonMixDistances # distances to pure strains in non-mixed sites
  amb<-c("N","?","-") # ambiguous calls
  maxDistance<-maxDistance
  
  # distance in mixed sites if either allele is considered
  distance<-numeric(ncol(df)) 
  ambsInStrains<-numeric(ncol(df))
  for (m in 1:ncol(df)){
    ambsInStrains[m]<-length(which(df[,m] %in% amb))
    distance[m]<-(nrow(df)-length(which(unlist(Map(function(x, y) all(x %in% y), df[,m], Nuc))==T)))-ambsInStrains[m]
  }
  
  newdist<-distance+nonmix_dist   # full distance to pure strains from all sites
  distance_DF<-data.frame(Name=colnames(df),nonMixDistance=nonmix_dist,
                          FullDistance=newdist,MedianProps=0,PropHigherProps=0,nb.Amb=ambsInStrains) # Data frame of all distances
  for (m in 1:ncol(df)){
    result <- sapply(seq_along(df[,m]), function(p) {
      match(df[p,m], Nuc[[p]])
    })
    AD_SNP <- sapply(seq_along(result), function(i) AD[[i]][result[i]])
    AD_total<- sapply(AD, sum)
    
    AD_props<-AD_SNP[!is.na(AD_SNP)]/AD_total[!is.na(AD_SNP)]
    distance_DF$MedianProps[m]<-median(AD_props)
    distance_DF$PropHigherProps[m]<-length(which(AD_props>0.5))/length(AD_props)
  }
  
  if (nrow(distance_DF)>=1){
    return(list(Summary=distance_DF))
  } else {
    return(list())
  }
}


hammingDistance <- function(sequence1, sequence2) {
  sequence1<-sequence1
  sequence2<-sequence2
  amb<-c("N","?","-")
  missing<-unique(c(which(sequence1 %in% amb),which(sequence2 %in% amb)))
  if (length(missing)>0){
    sequence1<-sequence1[-missing]
    sequence2<-sequence2[-missing]
  }
  distance<-length(which(sequence1 != sequence2))
  return(distance)
}

makeconsensus<- function(row){
  amb<-c("N","?","-") # ambiguous calls
  unique_in_row<-unique(row[2:length(row)])
  unique_in_row<-unique_in_row[!unique_in_row %in% amb]
  if (length(unique_in_row)==1) {
    return(unique_in_row)
  } else {
    if (!as.character(row[1]) %in% amb){
      return(as.character(row[1]))
    } else {
      return("N")
    }
  }
}

if (!require("stringr")){install.packages("stringr")}
if (!require("seqinr")){install.packages("seqinr")}
if (!require("optparse")) { install.packages("optparse") }
require(seqinr)
require(stringr)
library(optparse)
options(stringsAsFactors = F)

reconstructConstituents <- function(VCFfile, prefix = "output", MixInfect2Result, 
                                    maskFile = NULL, minQual = 20, maxDistance = 5000,
                                    LowCov = 10, popFreq_threshold = 1, minDepth = 5,
                                    mixProp = 0.9, closestStrain = T, n_threads = 4) {
  
  amb<-c("N","?","-") # ambiguous sites
  ## Read VCF
  allvcf<-read.table(VCFfile) # read in VCF
  names<- as.matrix(read.table(VCFfile,comment.char=" ",sep="\n"))
  end_head<-which(grepl("#CHROM",names)==TRUE)
  names<-unlist(strsplit(names[end_head],"\t"))
  colnames(allvcf)<-names
  head_start<-names[1:9]
  format <- which(names == 'FORMAT')
  samplenames<-names[10:length(names)] # get all sample names in VCF
  
  ## Remove indels
  ind <- lapply(1:nrow(allvcf), function(i) {
    length(unlist(strsplit(allvcf[i, 4], ""))) > 1 || length(unlist(strsplit(unlist(strsplit(allvcf[i, 5], ","))[1], ""))) > 1
  })
  allvcf <- allvcf[which(ind==FALSE), ]
  
  ## Read MixInfect2 result
  mixedSummary<-read.csv(MixInfect2Result)
  mixednames<-mixedSummary[which(mixedSummary[,2]=="Mix"),1]
  nonmixed<-mixedSummary[which(mixedSummary[,2]=="Non-mix"),1]
  
  ## Filter variants using QC
  allvcf<-allvcf[allvcf[,which(head_start=="FILTER")]=="PASS",] # Remove non-pass variants
  allvcf<-allvcf[allvcf[,which(head_start=="QUAL")]>minQual,] # Remove low quality variants
  allvcf<-allvcf[grep("*",allvcf[,which(head_start=="ALT")],fixed = T,invert = T),] # Remove spanning variants
  
  # Remove SNPs in repeats/PPE etc., QC pass, and overlap variants
  regions2remove<-read.csv(maskFile)
  pos<-numeric()
  for (i in 1:nrow(regions2remove)){
    pos<-c(pos,regions2remove$START[i]:regions2remove$STOP[i])
  }
  remove<-which(allvcf[,2] %in% pos)
  if (length(remove)>0){
    allvcf<-allvcf[-remove,]
  }
  
  #### Determine AD field and create matrix of AD and GT
  AD<-which(unlist(str_split(allvcf[1,format], ":"))=='AD')
  GT<-which(unlist(str_split(allvcf[1,format], ":"))=='GT')
  DP<-which(unlist(str_split(allvcf[1,format], ":"))=='DP')
  
  #### Make new matrices of separated GT and AD fields
  GT_mat<-matrix(sapply(str_split(allvcf[,format], ":"), "[[", GT))
  DP_mat<-matrix(sapply(str_split(allvcf[,format], ":"), "[[", DP))
  AD_mat<-matrix(sapply(str_split(allvcf[,format], ":"), "[[", AD))
  for (i in (format+1):ncol(allvcf)){
    GT_mat<-cbind(GT_mat,sapply(str_split(allvcf[,i], ":"), "[[", GT))
    newDP<-str_split(allvcf[,i], ":")
    newDP <- lapply(newDP, function(x) if (length(x) < 3) c(x, 0) else x)
    DP_mat<-cbind(DP_mat,sapply(newDP, "[[", DP))
    AD_mat<-cbind(AD_mat,sapply(str_split(allvcf[,i], ":"), "[[", AD))
  }
  GT_mat<-GT_mat[,-1, drop=F]
  AD_mat<-AD_mat[,-1, drop=F]
  DP_mat<-DP_mat[,-1, drop=F]
  DP_mat[which(DP_mat==".")]<-"0"
  GT_mat[which(as.numeric(DP_mat)<LowCov)]<-"?"
  
  mixed_calls<-c("0/1","0/2","0/3","1/2","1/3","2/3","0|1","0|2","0|3","1|2","1|3","2|3")
  alt_calls<-c("1/1","2/2","3/3","1|1","2|2","3|3")
  
  ## Which sites have mixed reads (not just 0/1 etc.)
  for (col in 1:ncol(AD_mat)){
    ADmix<-str_split(AD_mat[,col], ",")
    for (m in 1:length(ADmix)){
      AD_site<-as.numeric(unlist(ADmix[m],","))
      AD_site<-AD_site[!AD_site==0]
      if (length(AD_site)>1){
        AD_site<-AD_site[order(AD_site,decreasing = T)]
        if (AD_site[2]>=minDepth){
          GT_mat[m,col]<-"0/1"
        }
      }
    }
  }
  
  #### Ignore  sites with mixed frequency over popFreq_threshold
  propMix<-numeric()
  for (i in 1:nrow(GT_mat)){
    propMix[i]<-length(which(GT_mat[i,] %in% mixed_calls))/ncol(GT_mat)
  }
  if (any(propMix>popFreq_threshold)){
    GT_mat<-GT_mat[-which(propMix>popFreq_threshold),]
    allvcf<-allvcf[-which(propMix>popFreq_threshold),]
  }
  
  ## remove sites without alternative or mixed call
  remove<-apply(as.data.frame(GT_mat), 1, function(row) {
    any(row %in% c(mixed_calls,alt_calls))
  })
  allvcf<-allvcf[remove,]
  
  # Final vcfs of mixed and non-mixed strains
  nonmixvcf<-allvcf[,c(1:9,which(colnames(allvcf) %in% nonmixed))]
  mixvcf<-allvcf[,c(1:9,which(colnames(allvcf) %in% mixednames))]
  ref<-nonmixvcf[,4]
  alt<-str_split(nonmixvcf[,5], ",")
  rm(allvcf)
  
  if (length(nonmixed)>0){
    ## Make consensus sequence of all non-mixed strains
    finalnuc_nonmix<-as.data.frame(matrix(nrow = nrow(nonmixvcf),ncol=ncol(nonmixvcf)-8))
    finalnuc_nonmix[,1]<-nonmixvcf$POS
    colnames(finalnuc_nonmix)<-c("Position",colnames(nonmixvcf[10:ncol(nonmixvcf)]))
    
    
    for (i in 2:ncol(finalnuc_nonmix)){
      col<-i+8
      finalseq_maj<-ref
      geno<-sapply(str_split(nonmixvcf[,col], ":"), "[[", 1)
      firstref<-which(geno=="1/1" | geno=="1|1")
      secondref<-which(geno=="2/2" | geno=="2|2")
      thirdref<-which(geno=="3/3" | geno=="3|3")
      if (length(firstref)>0){
        finalseq_maj[firstref]<-sapply(alt[firstref], "[[", 1)
      }
      if (length(secondref)>0){
        finalseq_maj[secondref]<-sapply(alt[secondref], "[[", 2)
      }
      if (length(thirdref)>0){
        finalseq_maj[thirdref]<-sapply(alt[thirdref], "[[", 3)
      }
      
      mix_sites<-c(which(geno=="0/1" | geno=="0|1" | geno=="1/2" | geno=="1|2" |
                           geno=="1/3" | geno=="1|3" | geno=="2/3" | geno=="2|3"))
      
      for (j in mix_sites){
        genomix<-geno[j]
        nucs<-unlist(str_split(genomix,"|" ))
        nucs<-nucs[c(2,4)]
        AD<-sapply(str_split(nonmixvcf[j,col], ":"), "[[", 2)
        AD<-as.numeric(unlist(str_split(AD,",")))
        AD_nucs<-c(0,0)
        if (nucs[1]=="0"){
          nucs[1]<-ref[j]
          AD_nucs[1]<-AD[1]
        } else if (nucs[1]=="1") {
          nucs[1]<-sapply(alt[j], "[[", 1)
          AD_nucs[1]<-AD[2]
        } else if (nucs[1]=="2") {
          nucs[1]<-sapply(alt[j], "[[", 2)
          AD_nucs[1]<-AD[3]
        }
        if (nucs[2]=="1"){
          nucs[2]<-sapply(alt[j], "[[", 1)
          AD_nucs[2]<-AD[2]
        } else if (nucs[2]=="2") {
          nucs[2]<-sapply(alt[j], "[[", 2)
          AD_nucs[2]<-AD[3]
        } else if (nucs[1]=="3") {
          nucs[2]<-sapply(alt[j], "[[", 3)
          AD_nucs[2]<-AD[4]
        } 
        if (AD_nucs[1]==AD_nucs[2]){
          majnuc<-"N"
        } else if (max(AD_nucs)>=sum(AD)*mixProp) {
          majnuc<-nucs[which(AD_nucs==max(AD_nucs))]
        } else {
          majnuc<-"N"
        }
        finalseq_maj[j]<-majnuc
      }
      
      DP<-sapply(strsplit(nonmixvcf[,col],split = ":"), function(x) x[3])
      finalseq_maj[which(DP<10)]<-"?"
      finalnuc_nonmix[,i]<-finalseq_maj ## Final sequences of non-mixed strains
    }
  }
  
  
  ###### Get sequences of mixed strains
  if (length(mixednames)>0){
    finalnuc_mix<-as.data.frame(matrix(nrow = nrow(mixvcf),ncol=ncol(mixvcf)-8))
    finalnuc_mix[,1]<-mixvcf$POS
    finalnuc_min<-as.data.frame(mixvcf$POS)
    colnames(finalnuc_mix)<-c("Position",colnames(mixvcf[10:ncol(mixvcf)]))
    names_min<-"Position"
    if (length(nonmixed)>0){
      finalnucs<-finalnuc_nonmix
    } else {
      finalnucs<-data.frame(matrix(ncol = 1, nrow = nrow(mixvcf)))
      finalnucs[,1]<-mixvcf$POS
      colnames(finalnucs)[1]<-"Position"
    }
    
    ## Closest strain output 
    if (closestStrain){
      closestStrainSummary<-data.frame(Name=mixednames,Closest.major.strains=NA,Closest.major.strains.distance=NA,
                                       Closest.minor.strains=NA,Closest.minor.strains.distance=NA)
    }
    
    # Reconstruct major and minor strains 
    for (i in 2:ncol(finalnuc_mix)){
      mixname<-colnames(finalnuc_mix)[i]
      col<-i+8
      finalseq_maj<-ref
      finalseq_min<-ref
      geno<-sapply(str_split(mixvcf[,col], ":"), "[[", 1)
      
      ## Which sites have mixed reads (not just 0/1 etc)
      ADmix<-sapply(str_split(mixvcf[,col], ":"), "[[", 2)
      newmixsites<-numeric()
      for (m in 1:length(ADmix)){
        AD_site<-as.numeric(unlist(strsplit(ADmix[m],",")))
        AD_site_new<-AD_site[!AD_site==0]
        if (length(AD_site_new)>1){
          AD_site_new<-AD_site_new[order(AD_site_new,decreasing = T)]
          if (AD_site_new[2]>=minDepth | geno[m] %in% mixed_calls){
            newmixsites<-c(newmixsites,m)
            geno[m]<-paste0(which(AD_site %in% AD_site_new[1:2])-1,collapse = "/")
          } else {
            geno[m]<-paste0(c(which(AD_site %in% AD_site_new[1])-1,
                              which(AD_site %in% AD_site_new[1])-1),collapse = "/")
          }
        }
      }
      
      firstref<-which(geno=="1/1" | geno=="1|1")
      firstref<-firstref[!firstref %in% newmixsites]
      secondref<-which(geno=="2/2" | geno=="2|2")
      secondref<-secondref[!secondref %in% newmixsites]
      thirdref<-which(geno=="3/3" | geno=="3|3")
      thirdref<-thirdref[!thirdref %in% newmixsites]
      if (length(firstref)>0){
        finalseq_maj[firstref]<-sapply(alt[firstref], "[[", 1)
        finalseq_min[firstref]<-sapply(alt[firstref], "[[", 1)
      }
      if (length(secondref)>0){
        finalseq_maj[secondref]<-sapply(alt[secondref], "[[", 2)
        finalseq_min[secondref]<-sapply(alt[secondref], "[[", 2)
      }
      if (length(thirdref)>0){
        finalseq_maj[thirdref]<-sapply(alt[thirdref], "[[", 3)
        finalseq_min[thirdref]<-sapply(alt[thirdref], "[[", 3)
      } ## these lists have made a consensus sequence at non-mixed sites for both major and minor strains
      
      # Mask any sites with read depth < LowCov
      DP<-sapply(strsplit(mixvcf[,col],split = ":"), function(x) x[3])
      DP_low<-which(DP < LowCov)
      finalseq_maj[which(DP<10)]<-"?"
      finalseq_min[which(DP<10)]<-"?"
      
      if (length(nonmixed)>0){
        ## Distance to non-mixed strains in homozygous sites
        compare<-cbind(finalseq_maj,finalnuc_nonmix[,2:ncol(finalnuc_nonmix)])
        compare<-compare[-newmixsites,]
        if (length(nonmixed)==1){
          nonMixDistances<-hammingDistance(as.character(compare[,1]), as.character(compare[,2]))
        } else {
          nonMixDistances <- sapply(compare[,2:ncol(compare), drop=F], 
                                    function(x) hammingDistance(as.character(compare[,1]), as.character(x)))
        }
      }
      nuclist<-list()
      ADlist<-list()
      remove<-numeric()
      if (length(newmixsites)>0){ # After QC, check there are heterozygous sites remaining
        for (j in 1:length(newmixsites)){ 
          genomix<-geno[newmixsites[j]] 
          nucs<-unlist(str_split(genomix,"|" ))
          nucs<-nucs[c(2,4)]
          AD<-sapply(str_split(mixvcf[newmixsites[j],col], ":"), "[[", 2)
          AD<-as.numeric(unlist(str_split(AD,",")))
          AD_nucs<-c(0,0)
          
          if (nucs[1]=="0"){
            nucs[1]<-ref[newmixsites[j]]
            AD_nucs[1]<-AD[1]
          } else if (nucs[1]=="1") {
            nucs[1]<-sapply(alt[newmixsites[j]], "[[", 1)
            AD_nucs[1]<-AD[2]
          } else if (nucs[1]=="2") {
            nucs[1]<-sapply(alt[newmixsites[j]], "[[", 2)
            AD_nucs[1]<-AD[3]
          }
          if (nucs[2]=="1"){
            nucs[2]<-sapply(alt[newmixsites[j]], "[[", 1)
            AD_nucs[2]<-AD[2]
          } else if (nucs[2]=="2") {
            nucs[2]<-sapply(alt[newmixsites[j]], "[[", 2)
            AD_nucs[2]<-AD[3]
          } else if (nucs[2]=="3") {
            nucs[2]<-sapply(alt[newmixsites[j]], "[[", 3)
            AD_nucs[2]<-AD[4]
          } 
          
          nuclist[[j]]<-nucs # List of possible nucleotides at heterozygous sites
          ADlist[[j]]<-AD_nucs # List of possible read depths at heterozygous sites
        }
        
        ## make consensus from read frequencies only
        consensusMajor<-finalseq_maj
        consensusMinor<-finalseq_min
        for (j in 1:length(newmixsites)){
          AD_nucs<-unlist(ADlist[j])
          nucs<-unlist(nuclist[j])
          if (length(which(AD_nucs==min(AD_nucs)))>1){
            consensusMajor[newmixsites[j]]<-"N"
            consensusMinor[newmixsites[j]]<-"N"
          } else {
            if (AD_nucs[which(AD_nucs==max(AD_nucs))] >=minDepth){
              consensusMajor[newmixsites[j]]<-nucs[which(AD_nucs==max(AD_nucs))]
            } else {
              consensusMajor[newmixsites[j]]<-"N"
            }
            if (AD_nucs[which(AD_nucs==min(AD_nucs))] >= minDepth){
              consensusMinor[newmixsites[j]]<-nucs[which(AD_nucs==min(AD_nucs))]
            } else {
              consensusMinor[newmixsites[j]]<-"N"
            }
          }
        }
        finalnucs<-cbind(finalnucs,data.frame(consensusMajor,consensusMinor))
        colnames(finalnucs)[c((ncol(finalnucs)-1):ncol(finalnucs))]<-
          c(paste0(mixname,"_major_consensus"),paste0(mixname,"_minor_consensus"))
        
        
        ## Make closest strain constituents if TRUE
        if (closestStrain){
          if (length(nonmixed)==0){
            print("No non-mixed strains, closest strain method not possible")
          } else {

            # Find closest strain(s) in non-mixed strains using all sites 
            non_mixed_strainDF<-finalnuc_nonmix[newmixsites,2:ncol(finalnuc_nonmix), drop=F]
            closeststrains<-findRelatedStrains(non_mixed_strainDF,nonMixDistances,
                                               ADlist,nuclist,maxDistance)
            
            ## Pick closest Strains of higher and lower read frequencies
            Summary<-closeststrains$Summary
            ClosestSumHigher<-Summary[which(Summary$PropHigherProps>0.5),]
            if (nrow(ClosestSumHigher)>0){
              ClosestSumHigher<-ClosestSumHigher[which(ClosestSumHigher$FullDistance==min(ClosestSumHigher$FullDistance)),]
              ClosestSumHigher<-ClosestSumHigher[which(ClosestSumHigher$nb.Amb==min(ClosestSumHigher$nb.Amb)),]
              if (nrow(ClosestSumHigher)>1){
                closestStrainsMajor<-finalnuc_nonmix[,which(colnames(finalnuc_nonmix) %in% ClosestSumHigher$Name)]
                consensusMajor <- apply(closestStrainsMajor, 1, makeconsensus)
              } else {
                consensusMajor<-as.character(finalnuc_nonmix[,which(colnames(finalnuc_nonmix) %in% ClosestSumHigher$Name)])
              }
              closestStrainSummary$Closest.major.strains[i-1]<-paste0(ClosestSumHigher$Name,collapse = ",")
              closestStrainSummary$Closest.major.strains.distance[i-1]<-paste0(ClosestSumHigher$FullDistance,collapse = ",")
              finalnucs<-cbind(finalnucs,consensusMajor)
              colnames(finalnucs)[ncol(finalnucs)]<-paste0(mixname,"_major_closest")
            }
            
            ClosestSumLower<-Summary[which(Summary$PropHigherProps<0.5),]
            if (nrow(ClosestSumLower)>0){
              ClosestSumLower<-ClosestSumLower[which(ClosestSumLower$FullDistance==min(ClosestSumLower$FullDistance)),]
              ClosestSumLower<-ClosestSumLower[which(ClosestSumLower$nb.Amb==min(ClosestSumLower$nb.Amb)),]
              if (nrow(ClosestSumLower)>1){
                closestStrainsMinor<-finalnuc_nonmix[,which(colnames(finalnuc_nonmix) %in% ClosestSumLower$Name)]
                consensusMinor <- apply(closestStrainsMinor, 1, makeconsensus)
              } else {
                consensusMinor<-as.character(finalnuc_nonmix[,which(colnames(finalnuc_nonmix) %in% ClosestSumLower$Name)])
              }
              closestStrainSummary$Closest.minor.strains[i-1]<-paste0(ClosestSumLower$Name,collapse = ",")
              closestStrainSummary$Closest.minor.strains.distance[i-1]<-paste0(ClosestSumLower$FullDistance,collapse = ",")
              finalnucs<-cbind(finalnucs,consensusMinor)
              colnames(finalnucs)[ncol(finalnucs)]<-paste0(mixname,"_minor_closest")
            }
          }
        }
      } else {
        print(paste0(mixname,"is not a mix after additional QC steps"))
      }
    }
    
    ## Make FASTA of consensus and closest constituent strains
    unique_counts <- apply(finalnucs[,2:ncol(finalnucs)], 1, function(row) length(unique(row[!row %in% amb])))
    finalnucs<-finalnucs[which(!unique_counts<2),]
    output_fast<-t(finalnucs[,2:ncol(finalnucs)])
    forfastaref<-as.list(apply(output_fast, 1, paste, collapse=""))
    write.fasta(forfastaref,colnames(finalnucs)[2:ncol(finalnucs)],paste0(prefix,"_constituents.fasta"),open="w")
    write.table(data.frame(SNP=1:nrow(finalnucs),Position=finalnucs$Position),paste0(prefix,"_SNPindex.txt"),quote = F,sep = "\t",row.names = F)
    if (closestStrain){
      write.csv(closestStrainSummary,paste0(prefix,"_closestStrainSummary.csv"),row.names = F)
    }
  } else {
    print("No mixed infection from MixInfect2")
  }
}

option_list <- list(
  make_option(c("-v","--VCFfile"), type = "character", help = "Input VCF file (same as for MixInfect2)", metavar = "character"),
  make_option(c("-o", "--prefix"), type = "character", default = "output", help = "Prefix for output files", metavar = "character"),
  make_option(c("-r","--MixInfect2Result"), type = "character", help = "output CSV file generated using MixInfect2", metavar = "character"),
  make_option(c("-f","--maskFile"), type = "character", default = NULL, help = "CSV file with start and stop coordinates for regions to remove variants from PE/PPE, repeat regions etc. (recommended)", metavar = "character"),
  make_option(c("-q","--minQual"), type = "integer", default = 20, help = "Minimum per loci quality ", metavar = "integer"),
  make_option(c("-l","--LowCov"), type = "numeric", default = 10, help = "Minimum read depth at site to call either a cSNP or hSNP allele frequency", metavar = "numeric"),
  make_option(c("-c","--closestStrain"), type = "logical", default = TRUE, help = "Reconstruct constituent strains using closest strain method", metavar = "logical"),
  make_option(c("-x","--maxDistance"), type = "integer", default = 5000, help = "Maximum distance to closest non-mixed strain (If -c is TRUE)", metavar = "integer"),
  make_option(c("-p","--popFreq_threshold"), type = "numeric", default = 1, help = "Remove hSNPs found in greater than this proportion of sequences in VCF (set as 1 for single sample VCF)", metavar = "numeric"),
  make_option(c("-d","--minDepth"), type = "integer", default = 10, help = "Minimum read depth of minor frequency allele for a mixed call", metavar = "integer"),
  make_option(c("-m","--mixProp"), type = "numeric", default = 0.9, help = "Minimum read depth at site to call either a cSNP or hSNP allele frequency ", metavar = "numeric"),
  make_option(c("-t","--n_threads"), type = "integer", default = 4, help = "Number of threads to use", metavar = "integer")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if input files are provided
if (is.null(opt$VCFfile) | is.null(opt$MixInfect2Result)) {
  print_help(opt_parser)
  stop("Requires both a VCF and the CSV file from results of MixInfect2.", call. = FALSE)
}

# Run the function with parsed options
reconstructConstituents(opt$VCFfile, opt$prefix, opt$MixInfect2Result, opt$maskFile, 
                        opt$minQual, opt$LowCov, opt$closestStrain, opt$maxDistance, opt$popFreq_threshold, opt$minDepth, 
                        opt$mixProp, opt$n_threads)

