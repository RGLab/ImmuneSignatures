# Evan Henrich
# January 2017
# Fred Hutchinson Cancer Research Center
# Gottardo Lab
# ehenrich@fredhutch.org

#--------GEN NOTES-------------------------------------------------------
# This script selects the highest gene expression probe within a set of probes
# defined as the same geneSymbol.  Returns the Expression matrix with rownames
# now as geneSymbols instead of probe IDs.  This function is from the
# original work by the collaborators.

#' Function to select only the probes with highest expression for each gene and
#' remove the rest from the expression matrix
#'
#' @param eset ExpressionSet object that has expression matrix in AssayData()
#' @param gene_symbols Character vector that matches rownames of exprs(eset)
#' @export
#'
adjust_GE <- function(eset){

  # Ensure that expression matrix includes only genes from annotation dataframe
  expr_mx <- as.data.frame(exprs(eset))
  expr_mx <- expr_mx[ order(rownames(expr_mx)), ]
  fData <- fData(eset)
  fData <- fData[ which(fData[,1] %in% rownames(expr_mx)), ]
  fData <- fData[ order(fData[,1]), ]

  # Determine average expression for each probe; must come before gene symbol cleaning work
  expr_mx$prb_avg <- apply(expr_mx, 1 , function(x){ sum(2^x) / length(x) })

  # Remove duplicate rows, rows with NA values, and those with blank / NA gene symbols
  expr_mx$gs <- fData[,2]
  expr_mx <- expr_mx[ !duplicated(expr_mx), ] # SDY80 has exact dup rows
  expr_mx <- expr_mx[ which(!is.na(expr_mx$gs)), ] # SDY212 has NA vals for gene_symbols
  expr_mx <- expr_mx[ which(expr_mx$gs != ""), ] # SDY63 / SDY404 have "" for gs with SDY63_anno_tbl
  NA_rows <- which(rowSums(is.na(expr_mx)) > 0)
  if(length(NA_rows) > 0 ){ expr_mx <- expr_mx[ -NA_rows, ] }

  # Select and keep only the probe with maximum average expression per gene symbol
  expr_mx <- expr_mx %>%
    group_by(gs) %>%
    filter(prb_avg == max(prb_avg)) %>%
    ungroup()

  # Ensure that rownames are gene symbols
  expr_mx <- data.frame(expr_mx)
  rownames(expr_mx) <- expr_mx$gs
  expr_mx <- expr_mx[ , -which(colnames(expr_mx) %in% c("prb_avg", "gs"))]

  return(expr_mx)
}

