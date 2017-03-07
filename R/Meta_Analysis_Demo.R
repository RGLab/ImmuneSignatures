#########################################################################################################
## The script performed gene module meta-analysis on discovery cohorts to detect signiture gene
## modules significantly defferent between high responders and low responders to flu vaccination
## in HIPC signature project. The identifiedy gene modules were then verified by the validation cohort.
#########################################################################################################
##
## Original Author: Hailong Meng
## Updated last by Hailong: 2016-11-18
## ? 2016 Yale University. All rights reserved.
## Refactoring Author: Evan henrich (Fred Hutch), March 2017
#########################################################################################################

# EH NOTE:
# This version is heavily refactored and much of the eset adjustment is done in pipeline_demo.Rmd.

#---------HELPER METHODS--------------------------------------------------

adjust_eset <- function(eset, cohort){

  ## change from "lockedEnvironment" to "list" to modify assayData
  storageMode(assayData(eset)) <- "list"

  ## Select phenoData for correct cohort and remove other values (e.g. all / old if young)
  hai_var <- grep("^((fc)|(d0))_", varLabels(eset), value = TRUE)

  if( cohort == 'young'){
    oppo <- "old"
    cap_cohort <- "Young"
  }else{
    oppo <- "young"
    cap_cohort <- "Old"
  }

  eset <- eset[, eset$Age.class %in% cohort] # set in adjust_hai
  phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
                                               c(hai_var, paste0(oppo, "_", hai_var)))]
  varLabels(eset) <- gsub(paste0("^", cohort, "_"), '', varLabels(eset))

  ## change back to "lockedEnvironment" to prevent chagnes
  storageMode(assayData(eset)) <- "lockedEnvironment"

  ## result
  return(eset)
}

#***********MAIN METHOD***************************************************************
#' Function to perform meta analysis for HIPC ImmuneSignatures Project
#'
#' @param geneSetDB table defining gene sets
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

  result_dfs <- list() # holds output tables for use with markdown

  discoverySDY = c('SDY212','SDY63','SDY404','SDY400')
  ifelse(cohort == 'young', validation.sdy <- "SDY80", validation.sdy <- "SDY67")

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

    eset <- adjust_eset(eset_list[[sdy]], cohort)
    rownames(eset) <- as.character(fData(eset)$geneSymbol)
    labels <- as.character(pData(eset)[,endPoint])
    quSageObjList[[index.sdy]] <- qusage(eset, labels, "2-0", geneSetDB)
  }

  ## gene module meta analysis
  combinePDFsResult = combinePDFs(quSageObjList, n.points = 2^14)

  ## p values for gene module meta analysis
  combined.p <- pdf.pVal(combinePDFsResult)
  combined.q <- p.adjust(combined.p, method="BH")
  index_sig <- intersect(which(combined.q < FDR.cutoff), which(combined.p < pvalue.cutoff))

  ## get gene module activity
  combined.PDF <- combinePDFsResult$path.PDF
  pathway.activity <- c()

  for (i in 1:length(combined.p)){
    x.coordinates <- getXcoords(combinePDFsResult,i)
    tmp <- x.coordinates[which(combined.PDF[,i] == max(combined.PDF[,i],na.rm=TRUE))]
    pathway.activity <- c(pathway.activity, tmp)
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

  eset <- adjust_eset(eset_list[[validation.sdy]], cohort) # previously this used "all" if cohort was young
  labels <- as.character(pData(eset)[,endPoint])
  rownames(eset)  <- as.character(fData(eset)$geneSymbol)

  qs.results <- qusage(eset, labels, "2-0", geneSetDB)

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

#DEP

# if( cohort == 'young' ){
#   eset <- eset[, eset$Age.class %in% 'Young']
#   phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
#                                                c(hai_var, paste0("old_", hai_var)))]
#   varLabels(eset) <- gsub("^young_", '', varLabels(eset))
# }else if( group == 'old' ){
#   eset <- eset[, eset$Age.class %in% 'Old']
#   phenoData(eset) <- phenoData(eset)[, setdiff(varLabels(eset),
#                                                c(hai_var, paste0("young_", hai_var)))]
#   varLabels(eset) <- gsub("^old_", '', varLabels(eset))
# }


# ## define response variable: use extreme classes
# resp <- eset[[response]]
# lev <- levels(resp)
# lev <- c(head(lev, 1L), tail(lev, 1L))
# eset$Response <- factor(c('Non-responder', 'Responder')[match(resp, lev)])

# ## Combine by gene symbols and select probe with highest average gene expression
# run_qs <- function(eset, gene_symbols, labels, validation, geneSetDB){
#
#   ## Combine by gene symbols and select probe with highest average gene expression
#   ## rank probesets by average expression
#   nodup.rank <- apply(eset, 1 , function(x){ sum(2^x) / length(x) })
#   nodup.genes <- gene_symbols
#
#   for (p in 1:dim(eset)[1]) {
#     dups <- which(nodup.genes[p] == nodup.genes)      			##find the duplicate probesets
#     if (length(dups) > 1) {
#       # remove the 'best' probeset from the duplicated list and mark the duplicates with NA
#       nodup.genes[dups[-which.max(nodup.rank[dups])]] <- NA
#     }
#   }
#
#   ##remove the duplicated genes from results
#   eset.nodup <- eset[!is.na(nodup.genes),]
#   rownames(eset.nodup)  <- gene_symbols <- gene_symbols[!is.na(nodup.genes)]
#
#   if(validation){
#     eset.nodup.final <- eset.nodup
#   } # EH NOTE: only difference b/t discover and val
#
#   if(any(rownames(eset.nodup) == "")){
#     eset.nodup.final <- eset.nodup[-which(rownames(eset.nodup)==""),]
#   }else{
#     eset.nodup.final <- eset.nodup # done for case of Yale studies with manifest annotation
#   }
#
#   if(any(is.na(labels))){
#     eset.nodup.final <- eset.nodup.final[,-which(is.na(labels))]
#     labels <-  labels[-which(is.na(labels))]
#   }
#
#   sink("aux")
#   qs.results <-  qusage(eset.nodup.final, labels, "2-0", geneSetDB)
#   sink(NULL)
#
#   return(qs.results)
# }
# assayData(eset) <- list(exprs=assayData(eset)[["exprs"]]) # b/c NOT adjusted

#
# if(cohort == 'young'){
#   eset <- eset[, which(pData(eset)['Age'] < 36)]
# }
