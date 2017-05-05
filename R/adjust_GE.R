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
#' @export
#'
adjust_GE <- function(expr_mx){

  # Ensure that expression matrix includes only genes from annotation dataframe
  gs <- expr_mx$geneSymbol
  expr_mx <- expr_mx[ , -which(colnames(expr_mx) == "geneSymbol")]

  # Determine average expression for each probe; must come before gene symbol cleaning work
  expr_mx$prb_avg <- apply(expr_mx, 1 , function(x){ sum(2^x) / length(x) })

  # Remove duplicate rows, rows with NA values, and those with blank / NA gene symbols
  # expr_mx$gs <- fData[,2]  # uncomment after flat files are changed
  expr_mx$gs <- gs
  a <- !duplicated(expr_mx) # SDY80 has exact dup rows
  b <- !is.na(expr_mx$gs) # SDY212 has NA vals for gene_symbols
  c <- !(expr_mx$gs %in% c("","NA")) # SDY63 / SDY404 have "" and "NA"
  d <- rowSums(is.na(expr_mx)) == 0 # SDY212 has some rows with NA vals in expr
  expr_mx <- expr_mx[ (a & b & c & d) , ]

  # Select and keep only the probe with maximum average expression per gene symbol
  expr_mx <- expr_mx %>%
    group_by(gs) %>%
    filter(prb_avg == max(prb_avg)) %>%
    ungroup()

  # Ensure that rownames are gene symbols
  expr_mx <- data.frame(expr_mx, stringsAsFactors = F)
  expr_mx <- expr_mx[ !duplicated(expr_mx), ] # Dups possible after grouping in SDY67 ... not sure why
  rownames(expr_mx) <- expr_mx$gs
  expr_mx <- expr_mx[ , -which(colnames(expr_mx) %in% c("prb_avg", "gs"))]
}

