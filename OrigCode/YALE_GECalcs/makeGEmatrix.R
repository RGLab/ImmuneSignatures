#' ---
#' title: "Filtering low-quality probes from Yale cohort"
#' author: "Damian Fermin"
#' date: "March 31, 2015"
#' output: pdf_document
#' ---

library(Biobase);
library(lumi);
library(limma); ## for normalizeBetweenArrays()
library(ggplot2);

rm(list=ls());

minAvgIntensity = 0; ## use this to remove low intensity probes
minDetectionFraction = 0.75; ## 0-1 repesenting the percentage of d0 samples detecting a probe.


sdy63=c( ## year 1
  dataF="/home/dfermin/data/hipc/gea/sdy63.yale/illumina.data/SMOH_H12_SAMPLE_PROBE_PROFILE_NO_NORM_072011_year1.txt", 
  targetF="/home/dfermin/data/hipc/gea/sdy63.yale/illumina.data/Targets_PBMC_year1.txt",
  immPortF="/home/dfermin/data/hipc/gea/sdy63.yale/devel/sdy63_yaleID_to_ImmportID_mapping.tsv",
  outF="~/Dropbox/Yale/gea/remake.microarrays/sdy63.GEMatrix2.txt",
  RdataF="~/Dropbox/Yale/gea/remake.microarrays/sdy63_raw_data.RData");


sdy404=c( ## year 2
  dataF="/home/dfermin/data/hipc/gea/sdy404.yale/illumina.data/SMOH_H12_SAMPLE_PROBE_PROFILE_NONORM_091212_year2.txt", 
  targetF="/home/dfermin/data/hipc/gea/sdy404.yale/illumina.data/Targets_09132012_year2.txt",
  immPortF="/home/dfermin/data/hipc/gea/sdy404.yale/devel/sdy404_yaleID_to_ImmportID_mapping.tsv",
  outF="~/Dropbox/Yale/gea/remake.microarrays/sdy404.GEMatrix2.txt",
  RdataF="~/Dropbox/Yale/gea/remake.microarrays/sdy404_raw_data.RData");


sdy400=c( ## year 3
  dataF="/home/dfermin/data/hipc/gea/sdy400.yale/illumina.data/SMOH_H12_SAMPLE_PROBE_PROFILE_NO_NORM_091213_year3.txt", 
  targetF="/home/dfermin/data/hipc/gea/sdy400.yale/illumina.data/SMOH_H12_TARGET_year3.txt",
  immPortF="/home/dfermin/data/hipc/gea/sdy400.yale/devel/sdy400_yaleID_to_ImmportID_mapping.tsv",
  outF="~/Dropbox/Yale/gea/remake.microarrays/sdy400.GEMatrix2.txt",
  RdataF="~/Dropbox/Yale/gea/remake.microarrays/sdy400_raw_data.RData");


yale <- list(sdy63, sdy404, sdy400)


for(i in 1:3) {
  ##==============================================================================
  ##  Begin loop over yale list
  ##==============================================================================
  
  dataF = yale[[i]]["dataF"];
  targetF = yale[[i]]["targetF"];
  immPortF = yale[[i]]["immPortF"];
  outF = yale[[i]]["outF"];
  RdataF = yale[[i]]["RdataF"];
  
  if(!file.exists(RdataF)) {
    
    cat("\nReading in raw data files...\n");
    ptm <- proc.time();
    targetList <- read.delim(targetF, as.is=TRUE);
    
    if(length(grep("sdy63", outF)) > 0) { #' subject IDs were off in original data for year 1
      targetList$Subject <- targetList$Subject + 91000;
    }
    
    if(length(grep("sdy400", outF)) > 0) { #' column labels for year 3 were wrong
      colnames(targetList)[c(3,4,7)] <- c("Sample", "Subject", "Date");
    }
    
    cat("Reading in Illumina data...\n");
    x.lumi <- lumiR(dataF, checkDupId=FALSE); ## read in raw data
    
    y2import <- read.delim(immPortF, as.is=T, row.names=1); ## immPort mappings
    
    targetList$immportID = NA;
    for(yaleID in unique(targetList$Subject)) {
      yid = paste("SJ_",yaleID,sep="");
      immportID <- y2import[yid, ];
      idx <- grep(yaleID, targetList$Subject);
      targetList$immportID[idx] = immportID;
    }
    
    cat("Saving objects to '", RdataF, "'...\n");
    save(targetList, x.lumi, file=RdataF);
    deltaT <- proc.time() - ptm;
    cat("\nElapsed time: ", round(deltaT[3],2), " seconds\n");
  } else { 
    load(RdataF); 
  }
    
  ################################################################################
  
  ## get from the 'lumi object' d0 information
  mat <- exprs(x.lumi);
  ph <- pData(x.lumi);
  d0 <- targetList[ targetList$Date == 0, ]
  mat <- mat[, match(d0$Sample, colnames(mat))];
  ph <- ph[ match(d0$Sample, ph$sampleID), ];
  
  ## detection p-values are in the assayData slot.
  ## Get out the detection values only for samples from day 0
  detect <- assayData(x.lumi)[["detection"]];
  detect <- detect[, match(d0$Sample, colnames(detect))];
  
  ## we only want to know if the probe was detected on day 0
  tmp <- t(apply(detect, 1, function(x) as.numeric(x < 0.01)))
  colnames(tmp) <- colnames(detect);
  detectSum <- apply(tmp, 1, sum)
  
  
  ## probe mappings
  probesInfo <- data.frame(
    probes=fData(x.lumi)$PROBE_ID,
    gs=fData(x.lumi)$SYMBOL,
    detectCount=detectSum, ## number of samples probe had p-value below 0.01 on d0 
    stringsAsFactors=FALSE);
  
  expr <- log2(exprs(x.lumi)); ## expression matrix in log2 scale
  cat("\nQuantile normalizing expression data\n");
  expr <- normalizeBetweenArrays(as.matrix(expr), method="quantile"); ## normalize
  
  ## get new labels for the expression matrix
  for(i in seq(nrow(targetList))) {
    newColLab = paste(targetList$immportID[i],"_d", targetList$Date[i],sep="");
    j <- grep(targetList$Sample[i], colnames(expr));
    if(length(j) > 0) colnames(expr)[j] <- newColLab;
  }
  
  
  ## filter expr and probesInfo to only contain data for day0
  g <- grep("_d0", colnames(expr));
  expr <- expr[,g];
  targetList <- targetList[ targetList$Date == 0, ];
  
  ## keep all probes detected in at least 'k' samples at a p-value < 0.01
  N <- round((ncol(expr) * minDetectionFraction))  
  probesInfo <- probesInfo[ probesInfo$detectCount >= N,];
  expr <- expr[ rownames(probesInfo),]
  
  ## keep probes with an average expression across samples > minAvgIntensity
  rs <- apply(expr, 1, mean);
  expr <- expr[ rs >= minAvgIntensity, ]
  probesInfo <- probesInfo[ rownames(expr), ];
  
  ################################################################################

  out <- data.frame(
    probeID=probesInfo$probes,
    geneSymbol=probesInfo$gs,
    round(expr,6),
    stringsAsFactors=FALSE);
  
  write.table(out, file=outF, row.names=F, sep="\t", quote=F);
  
  rm(out, expr, targetList, probesInfo);
 
##==============================================================================
##  End loop over yale list
##==============================================================================
}




  