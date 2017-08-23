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
# This version is heavily refactored and much of the eset adjustment is done in ImmSig_Pipeline_Demo.Rmd.

#---------HELPER METHODS--------------------------------------------------
# Must sink to stop many print statements in qusage() from showing in report
run_qusage <- function(adj_eset, sdy, endPoint, gene_set){
  em <- data.frame(apply(exprs(adj_eset), 2, function(x){ as.numeric(x) } ))
  rownames(em) <- rownames(exprs(adj_eset))
  em$geneSymbol <- fData(adj_eset)[["gene_symbol"]]
  adj_em <- adjust_GE(em)
  labels <- pData(adj_eset)[, endPoint]

  sink("tmp")
  result <- qusage(adj_em, labels, "2-0", gene_set)
  sink(NULL)

  return(result)
}

#***********MAIN METHOD***************************************************************
#' Function to perform meta analysis for HIPC ImmuneSignatures Project
#'
#' @param eset_list list of expressionSet objects holding GE, HAI data
#' @param cohort Study cohort, young or old
#' @export

meta_analysis <- function(adj_eset_list, cohort, gene_set){

  # original manuscript params
  FDR.cutoff <- 0.508
  pvalue.cutoff <- 0.0105
  endPoint <- "fc_res_max_d30"

  results <- list() # holds output tables for use with markdown

  discoverySDY = c('SDY212','SDY63','SDY404','SDY400')
  validation.sdy <- ifelse(cohort == 'young', "SDY80", "SDY67")

  #########################################################################################################
  ##
  ## 2. Identify gene modules significantly different between high responders and low responders from
  ## discovery cohorts
  #########################################################################################################

  quSageObjList <- list() ## save gene module analysis result for each SDY
  idx <- 0

  for(sdy in discoverySDY){
    idx <- idx + 1
    quSageObjList[[idx]] <- results[[sdy]] <- run_qusage(adj_eset = adj_eset_list[[sdy]],
                                                         sdy = sdy,
                                                         endPoint = endPoint,
                                                         gene_set = gene_set)
  }

  ## gene module meta analysis
  combinePDFsResult <- results[["combinedPDF"]] <- combinePDFs(quSageObjList, n.points = 2^14)

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

  out_matrix <- cbind(Pvalue = round(combined.p, digits = 4),
                      FDR = round(combined.q, digits = 3),
                      pathwayActivity = round(pathway.activity, digits = 3))

  rownames(out_matrix) = colnames(combinePDFsResult$path.PDF)

  results$dsc <- as.data.frame(out_matrix)


  #########################################################################################################
  ##
  ## 3. Validate significantly different gene modules between high responders and low responders
  ##
  #########################################################################################################

  qs.results <- results[["valPDF"]] <- run_qusage(adj_eset = adj_eset_list[[validation.sdy]],
                                                        sdy = validation.sdy,
                                                        endPoint = endPoint,
                                                        gene_set = gene_set)

  pvalue <- pdf.pVal(qs.results)[index_sig]
  qvalue <- p.adjust(pvalue, method = "BH")
  pathway.activity.selected <- qs.results$path.mean[index_sig]

  out_matrix <- cbind(Pvalue = round(pvalue, digits = 4),
                      FDR = round(qvalue, digits = 3),
                      pathwayActivity = round(pathway.activity.selected, digits = 3))
  
  rownames(out_matrix) <- names(gene_set)[index_sig] # equal to selected.pathways

  results$val <- as.data.frame(out_matrix)

  return(results)
}
