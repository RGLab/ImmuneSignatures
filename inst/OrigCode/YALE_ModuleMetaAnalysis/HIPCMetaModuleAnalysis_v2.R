#########################################################################################################
## The script performed gene module meta-analysis on discovery cohorts to detect signiture gene
## modules significantly defferent between high responders and low responders to flu vaccination 
## in HIPC signature project. The identifiedy gene modules were then verified by the validation cohort.          
#########################################################################################################
##
## Author: Hailong Meng
## Updated: 2016-11-18
## ? 2016 Yale University. All rights reserved.
#########################################################################################################

rm(list=ls())
library("qusage")

#########################################################################################################
##
## 1. Setup parameters
## The working directory structure is:
##
## -Working direcory
##    --data
##        --SDY212.rds
##        --SDY63.rds
##        --SDY404.rds
##        --SDY400.rds
##        --CHI-nih.rds
##        --SDY67-batch2.rds
##        --BTM_for_GSEA_20131008.GMT
##        
##    --output
##    --HIPCMetaModuleAnalysis.R
##    --getSDY.R
##    
## data is the directory hosting all R objects and gene module file (BTM_for_GSEA_20131008.GMT)
## output is the directory containing all analysis results/plots from meta gene module analysis
#########################################################################################################
#setwd("C:/Hailong/Projects/HIPC/ImmuneSpace/") ## set up working directory

## Run gene module -meta-analysis for either young or older cohort
cohort = 'young'
# cohort ='old'
FDR.cutoff = 0.5 ## FDR cutoff for gene module analysis to call significant gene modules
pvalue.cutoff = 0.01 ## pvalue cutoff for gene module analysis to call significant gene modules

#########################################################################################################
## Discovery cohorts: SDY212, SDY63, SDY404, SDY400
## Validation Cohorts: SDY80 (Young) and SDY67 (Older)
#########################################################################################################
source('getSDY_swaps.R') ## getSDY is the script from Renaud to read R object files and make sure it is in working directory
getSDY <- .getSDY('data') ## all R objects should be in folder of "data" under working direcotory
discoverySDY = c('SDY212','SDY63','SDY404','SDY400')
if(cohort == 'young'){
  validation.sdy = 'CHI-nih' 
}else if(cohort == 'old'){
  validation.sdy ='SDY67-batch2'
}else{
  stop("The cohort name is not correct!")
}
endPoint = 'fc_res_max_d30' ## select definition of end point metrics


## read in gene modules
geneSetDB = readLines("data/BTM_for_GSEA_20131008.gmt") 
geneSetDB = strsplit(geneSetDB,"\t")                                            ##convert from vector of strings to a list
names(geneSetDB) = sapply(geneSetDB,"[",1)                               ##move the names column as the names of the list
links = sapply(geneSetDB,"[",2) 
geneSetDB = lapply(geneSetDB, "[",-1:-2)                                      ##remove name and description columns
geneSetDB= lapply(geneSetDB, function(x){x[which(x!="")]})           ##remove empty strings


#########################################################################################################
##
## 2. Identify gene modules significantly different between high responders and low responders from
## discovery cohorts
#########################################################################################################
quSageObjList=list() ## save gene module analysis result for each SDY
index.sdy = 0
## Perform gene module analysis for each discovery data set
for(sdy in discoverySDY){
  index.sdy =index.sdy+1 
  print(paste0("Processing ", sdy, " --- "))
  dataset_name=paste(sdy, cohort, endPoint, sep="_" )
  
  # read in data from R object
  eset <- getSDY(sdy, cohort, endPoint, adjusted=FALSE, baselineOnly=TRUE, swapSDY404 = TRUE)
  if(sdy=='SDY212' && cohort == 'young'){
    ## two subjects from SDY212 has the same subject ID "SUB134307", but Stanford don't know what happened to this
    ## so, we removed those two entries in data analysis
    print("Remove subject 'SUB134307' from data analysis")
    SUB134307_index=which(pData(eset)['SubjectID']=="SUB134307")
    eset=eset[,-SUB134307_index]
  }
  print(summary(pData(eset)[c('Age.class', 'Response', 'Condition')]) )
  gene_symbols  = as.character(fData(eset)$geneSymbol)
  labels <- as.character(pData(eset)[,endPoint])
  eset = exprs(eset)
  
  ## Combine by gene symbols and select probe with highest average gene expression
  nodup.rank <- apply(eset,1,function(x) sum(2^x)/length(x))  ##rank probesets by average expression
  nodup.genes <- gene_symbols
  for (p in 1:dim(eset)[1]) {
    dups <- which(nodup.genes[p]==nodup.genes)    				##find the duplicate probesets
    if (length(dups)>1) {
          nodup.genes[dups[-which.max(nodup.rank[dups])]] <- NA				##remove the 'best' probeset from the duplicated list and mark the duplicates with NA
    }										
  }
  eset.nodup <- eset[!is.na(nodup.genes),]					##remove the duplicated genes from results
  gene_symbols <- gene_symbols[!is.na(nodup.genes)]
  rownames(eset.nodup)  <- gene_symbols
  ## remove rows has empty gene names
  if(any(rownames(eset.nodup)=="")){
    eset.nodup.final <- eset.nodup[-which(rownames(eset.nodup)==""),]
  }
  ## remove samples without definition of response
  if(any(is.na(labels))){
    eset.nodup.final=eset.nodup.final[,-which(is.na(labels))]
    labels = labels[-which(is.na(labels))]
  }
  
  ## run gene module analysis. 2 is highResponder and 0 is lowResponder
  qs.results = qusage(eset.nodup.final, labels, "2-0", geneSetDB)

  ## get p values of each gene module from gene module analysis
  #p.vals = pdf.pVal(qs.results)

  ## save current gene module analysis for meta analysis later on
  quSageObjList[[index.sdy]] = qs.results 
  
  ## if would like to output gene module analysis result for each SDY
  #write.table(qsTable(qs.results, number=length(geneSetDB), sort.by = "none"), paste(dataset_name,".txt", sep=""), quote=FALSE,sep = "\t", row.names = FALSE)
  
  ## if would like to save R objects of gene module analysis for each SDY
  #rds.filename=paste("output/",dataset_name, ".rds",sep="")
  #saveRDS(qs.results, rds.filename)
}

## gene module meta analysis
combinePDFsResult = combinePDFs(quSageObjList,n.points=2^14)  

## p values for gene module meta analysis
combined.p = pdf.pVal(combinePDFsResult)
combined.q = p.adjust(combined.p, method="BH")

## get gene module activity
combined.PDF = combinePDFsResult$path.PDF
pathway.activity = c()
for (i in 1:length(combined.p)){
  x.coordinates = getXcoords(combinePDFsResult,i)
  tmp = x.coordinates[which(combined.PDF[,i]==max(combined.PDF[,i],na.rm=TRUE))]
  pathway.activity=c(pathway.activity, tmp)
}

out_matrix=cbind(Pvalue=combined.p, FDR=combined.q, pathwayActivity=pathway.activity)
rownames(out_matrix) = colnames(combinePDFsResult$path.PDF)
## output module meta-analysis results to the output folder
write.table(out_matrix, file=paste0("output/metaGeneModuleAnalysis_DiscoveryCohort_",cohort,".txt"), quote=F, sep = "\t", row.names=T, col.names=NA)

## output PDF figures for gene analysis results
index_sig = intersect(which(combined.q < FDR.cutoff), which(combined.p < pvalue.cutoff)) 
pdf(paste0("output/significantModules_DiscoveryCohort_",cohort,".pdf"), width=6, height=6)
for(i in index_sig){
  plot(combinePDFsResult,path.index=i)
  legend("topleft", legend=c(discoverySDY,"metaAnalysis"),
         lty=1, col=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","black"))  
  text(0.2, 1, paste("P value =", format(combined.p[i], digits=2),sep=""))
}
dev.off()

#########################################################################################################
##
## 3. Validate gene modules significantly different between high responders and low responders using
## validation cohort
#########################################################################################################

if(cohort == 'young'){
  print(paste0("Validate signiture gene modules for Young cohort SDY80 --- "))
  dataset_name=paste('SDY80', endPoint, sep="_" )
  eset <- getSDY(validation.sdy, 'all', endPoint, adjusted=FALSE, baselineOnly=TRUE)
  d0_index <- which(pData(eset)['Condition']=='d0')
  
  # get young samples only
  eset_d0<-eset[,d0_index]
  young_index <- which(pData(eset_d0)['Age']<36)
  eset_young<-eset_d0[,young_index]   
  eset<-eset_young
  print(summary(pData(eset)[c('Age.class', 'Response', 'Condition')]) )
  gene_symbols  = as.character(fData(eset)$gs)
  rownames(eset) = gene_symbols
  labels <- as.character(pData(eset)[,endPoint])
  
}else{
  print(paste0("Validate signiture gene modules for Old cohort SDY67_batch2 --- "))
  dataset_name=paste('SDY67_young', endPoint, sep="_" )
  eset <- getSDY('SDY67-batch2', cohort, endPoint, adjusted=FALSE, baselineOnly=TRUE)  
  print(summary(pData(eset)[c('Age.class', 'Response', 'Condition')]) )  
  gene_symbols  = as.character(fData(eset)$geneSymbol)
  labels <- as.character(pData(eset)[,endPoint])
  eset = exprs(eset)  
}

## Combine by gene symbols and select probe with highest average gene expression
nodup.rank <- apply(eset,1,function(x) sum(2^x)/length(x))  ##rank probesets by average expression
nodup.genes <- gene_symbols
for (p in 1:dim(eset)[1]) {
  dups <- which(nodup.genes[p]==nodup.genes)      			##find the duplicate probesets
  if (length(dups)>1) {
    nodup.genes[dups[-which.max(nodup.rank[dups])]] <- NA				##remove the 'best' probeset from the duplicated list and mark the duplicates with NA
  }										
}
eset.nodup <- eset[!is.na(nodup.genes),]					##remove the duplicated genes from results
gene_symbols <- gene_symbols[!is.na(nodup.genes)]
rownames(eset.nodup)  <- gene_symbols
eset.nodup.final <- eset.nodup
if(any(rownames(eset.nodup)=="")){
  eset.nodup.final <- eset.nodup[-which(rownames(eset.nodup)==""),]
}
if(any(is.na(labels))){
  eset.nodup.final=eset.nodup.final[,-which(is.na(labels))]
  labels = labels[-which(is.na(labels))]
}

## validate signature gene modules
qs.results = qusage(eset.nodup.final, labels, "2-0", geneSetDB)
p.vals = pdf.pVal(qs.results)
pathway.activity=qs.results$path.mean

selected.pathways =names(geneSetDB)[index_sig]
pvalue = p.vals[index_sig]
pathway.activity.selected = pathway.activity[index_sig]
qvalue = p.adjust(pvalue, method="BH")  

out_matrix=cbind(pvalue, qvalue,pathway.activity.selected)
rownames(out_matrix) = selected.pathways
write.table(out_matrix, file=paste0("output/metaGeneModuleAnalysis_ValidationCohort_",cohort,".txt"), quote=F, sep = "\t", row.names=T, col.names=NA)

pdf(paste0("output/significantModules_ValidationCohort_",cohort,".pdf"), width=6, height=6)
for(i in index_sig){ 
  plot(qs.results,path.index=i ,col="black",xlim=c(-1,1),ylim=c(0,10),xlab="Gene Module Activity",
         main=names(geneSetDB)[i] )  
  text(0.2, 1,paste("p value=", round(p.vals[i],digits=3),sep=" "))
  abline(v=0, lty=2)   
  }
dev.off()

