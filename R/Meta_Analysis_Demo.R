#########################################################################################################
## The script performed gene module meta-analysis on discovery cohorts to detect signiture gene
## modules significantly defferent between high responders and low responders to flu vaccination
## in HIPC signature project. The identifiedy gene modules were then verified by the validation cohort.
#########################################################################################################
##
## Original Author: Hailong Meng
## Updated last by Hailong: 2016-11-18
## ? 2016 Yale University. All rights reserved.
## Refactoring Author: Evan henrich (Fred Hutch), January 2017
## Note: All locations of substantive edits are commented with "EH NOTE:"
#########################################################################################################

# EH NOTE:
# Removed comments regarding directory structure and converted script into a function to allow
# arguments of cohort, FDR.cutoff, pvalue.cutoff, and endPoint so that users may test
# different possibilities

#---------HELPER METHODS--------------------------------------------------

adjust_eset <- function(eset, cohort, response){

  ## Limit to baseline
  eset <- eset[,grep("^((d0)|(dneg))", eset$Condition)]

  ## change from "lockedEnvironment" to "list" to modify
  storageMode(assayData(eset)) <- "list"

  ## massage phenotypic data and select which assayData to use (adjusted or not)
  hai_var <- grep("^((fc)|(d0))_", varLabels(eset), value = TRUE)
  assayData(eset) <- list(exprs=assayData(eset)[["exprs"]]) # b/c NOT adjusted

  if( cohort == 'young' ){
    eset <- eset[, eset$Age.class %in% 'Young']
    phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
                                                 c(hai_var, paste0("old_", hai_var)))]
    varLabels(eset) <- gsub("^young_", '', varLabels(eset))
  }else if( group == 'old' ){
    eset <- eset[, eset$Age.class %in% 'Old']
    phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
                                                 c(hai_var, paste0("young_", hai_var)))]
    varLabels(eset) <- gsub("^old_", '', varLabels(eset))
  }

  ## change back to "lockedEnvironment" to prevent chagnes
  storageMode(assayData(eset)) <- "lockedEnvironment"

  ## define response variable: use extreme classes
  resp <- eset[[response]]
  lev <- levels(resp)
  lev <- c(head(lev, 1L), tail(lev, 1L))
  eset$Response <- factor(c('Non-responder', 'Responder')[match(resp, lev)])

  ## result
  return(eset)
}

## Combine by gene symbols and select probe with highest average gene expression
run_qs <- function(eset, gene_symbols, labels, validation, geneSetDB){

  ## Combine by gene symbols and select probe with highest average gene expression
  ## rank probesets by average expression
  nodup.rank <- apply(eset, 1 , function(x){ sum(2^x) / length(x) })
  nodup.genes <- gene_symbols

  for (p in 1:dim(eset)[1]) {
    dups <- which(nodup.genes[p] == nodup.genes)      			##find the duplicate probesets
    if (length(dups) > 1) {
      # remove the 'best' probeset from the duplicated list and mark the duplicates with NA
      nodup.genes[dups[-which.max(nodup.rank[dups])]] <- NA
    }
  }

  ##remove the duplicated genes from results
  eset.nodup <- eset[!is.na(nodup.genes),]
  rownames(eset.nodup)  <- gene_symbols <- gene_symbols[!is.na(nodup.genes)]

  if(validation){
    eset.nodup.final <- eset.nodup
    } # EH NOTE: only difference b/t discover and val

  if(any(rownames(eset.nodup) == "")){
    eset.nodup.final <- eset.nodup[-which(rownames(eset.nodup)==""),]
  }else{
    eset.nodup.final <- eset.nodup # done for case of Yale studies with manifest annotation
  }

  if(any(is.na(labels))){
    eset.nodup.final <- eset.nodup.final[,-which(is.na(labels))]
    labels <-  labels[-which(is.na(labels))]
  }

  sink("aux")
  qs.results <-  qusage(eset.nodup.final, labels, "2-0", geneSetDB)
  sink(NULL)

  return(qs.results)
}


#########################################################################################################
## Discovery cohorts: SDY212, SDY63, SDY404, SDY400
## Validation Cohorts: SDY80 (Young) and SDY67 (Older)
#########################################################################################################
#' Function to perform meta analysis for HIPC ImmuneSignatures Project
#'
#' @param geneSetDB table defining gene sets
#' @param rds_data_dir Directory holding eset objects as rds files
#' @param cohort Study cohort, young or old
#' @param FDR.cutoff Cutoff for q-values in selecting significant gene sets
#' @param pvalue.cutoff cutoff for p-values in selecting significant gene sets
#' @param endPoint HAI table column used for categorizing response
#' @export
meta_analysis <- function(geneSetDB,
                          eset_list,
                          cohort,
                          FDR.cutoff = 0.5,
                          pvalue.cutoff = 0.01,
                          endPoint = 'fc_res_max_d30'){

  library(qusage) # load whole library b/c qusage::qusage() needs other functions
  result_dfs <- list() # holds output tables for use with markdown

  discoverySDY = c('SDY212','SDY63','SDY404','SDY400')
  if(cohort == 'young'){
    validation.sdy <- 'CHI-nih'
  }else if(cohort == 'old'){
    validation.sdy <- 'SDY67-batch2'
  }else{
    stop("The cohort name is not correct!")
  }

  # Parse Gene Module
  geneSetDB <- strsplit(geneSetDB,"\t")            ##convert from vector of strings to a list
  names(geneSetDB) <- sapply(geneSetDB,"[",1)      ##move the names column as the names of the list
  links <- sapply(geneSetDB,"[",2)
  geneSetDB <- lapply(geneSetDB, "[",-1:-2)        ##remove name and description columns
  geneSetDB <- lapply(geneSetDB, function(x){x[which(x!="")]})      ##remove empty strings


  #########################################################################################################
  ##
  ## 2. Identify gene modules significantly different between high responders and low responders from
  ## discovery cohorts
  #########################################################################################################

  quSageObjList <- list() ## save gene module analysis result for each SDY
  index.sdy <- 0

  for(sdy in discoverySDY){
    index.sdy = index.sdy + 1

    eset <- adjust_eset(eset_list[[sdy]], cohort, response = endPoint)

    if(sdy =='SDY212' && cohort == 'young'){
      ## two subjects from SDY212 has the same subject ID "SUB134307", but
      ## Stanford don't know what happened to this
      ## so, we removed those two entries in data analysis
      SUB134307_index <-  which(pData(eset)['SubjectID'] == "SUB134307")
      eset <- eset[,-SUB134307_index]
    }

    gene_symbols  <- as.character(fData(eset)$geneSymbol)
    labels <- as.character(pData(eset)[,endPoint])

    ## save current gene module analysis for meta analysis later on
    quSageObjList[[index.sdy]] <- run_qs(exprs(eset),
                                         gene_symbols,
                                         labels,
                                         validation = F,
                                         geneSetDB)
  }

  ## gene module meta analysis
  combinePDFsResult = combinePDFs(quSageObjList, n.points = 2^14)

  ## p values for gene module meta analysis
  combined.p = pdf.pVal(combinePDFsResult)
  combined.q = p.adjust(combined.p, method="BH")
  index_sig <- intersect(which(combined.q < FDR.cutoff), which(combined.p < pvalue.cutoff))

  ## get gene module activity
  combined.PDF = combinePDFsResult$path.PDF
  pathway.activity = c()
  for (i in 1:length(combined.p)){
    x.coordinates = getXcoords(combinePDFsResult,i)
    tmp = x.coordinates[which(combined.PDF[,i] == max(combined.PDF[,i],na.rm=TRUE))]
    pathway.activity=c(pathway.activity, tmp)
  }

  out_matrix <- cbind(Pvalue = combined.p,
                      FDR = combined.q,
                      pathwayActivity = pathway.activity)
  rownames(out_matrix) = colnames(combinePDFsResult$path.PDF)

  cat(paste0("DISCOVERY GROUP - SIGNIFICANT PATHWAY FIGURES"))

  result_dfs$dsc <- as.data.frame(out_matrix)

  for(i in index_sig){
    plot(combinePDFsResult, path.index = i)
    legend("topleft", legend=c(discoverySDY,"metaAnalysis"),
           lty=1, col=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","black"))
    text(0.2, 1, paste("P value =", format(combined.p[i], digits=2),sep=""))
  }



  #########################################################################################################
  ##
  ## 3. Validate significantly different gene modules between high responders and low responders
  ##
  #########################################################################################################

  eset <- adjust_eset(eset_list[[sdy]], cohort, response = endPoint) # previously this used "all" if cohort was young

  if(cohort == 'young'){
    eset <- eset[, which(pData(eset)['Condition'] == 'd0' & pData(eset)['Age'] < 36)]
    rownames(eset) <- as.character(fData(eset)$geneSymbol)
  }else{
    eset <- exprs(eset)
  }

  labels <- as.character(pData(eset)[,endPoint])
  gene_symbols  <- as.character(fData(eset)$geneSymbol)

  qs.results <- run_qs(eset, gene_symbols, labels, validation = T, geneSetDB)

  pvalue <- pdf.pVal(qs.results)[index_sig]
  qvalue <- p.adjust(pvalue, method = "BH")
  pathway.activity.selected <- qs.results$path.mean[index_sig]

  out_matrix <- cbind(pvalue, qvalue, pathway.activity.selected)
  rownames(out_matrix) <- names(geneSetDB)[index_sig] # equal to selected.pathways

  cat(paste0("VALIDATION STUDY - SIGNFICANT PATHWAY FIGURES"))

  result_dfs$val <- as.data.frame(out_matrix)

  # plot graphs
  for(i in index_sig){
    plot(qs.results,
         path.index = i,
         col = "black",
         xlim = c(-1,1),
         ylim = c(0,10),
         xlab = "Gene Module Activity",
         main = names(geneSetDB)[i]
    )
    text(0.2,
         1,
         paste("p value =",
               round(pdf.pVal(qs.results)[i], digits = 3),
               sep = " "))
    abline(v = 0, lty = 2)
  }

  return(result_dfs)
}


