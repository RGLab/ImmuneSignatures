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

adjust_GE <- function(eset, gene_symbols){

  expr_mx <- as.data.frame(exprs(eset))

  ## Combine by gene symbols and select probe with highest average gene expression
  ## rank probexpr_mxs by average expression
  nodup.rank <- apply(expr_mx, 1 , function(x){ sum(2^x) / length(x) })
  nodup.genes <- gene_symbols

  ## find the duplicate probes and use only highest expression, mark others NA
  for ( p in 1:dim(expr_mx)[1] ) {
    dups <- which(nodup.genes[p] == nodup.genes)
    if ( length(dups) > 1) {
      nodup.genes[dups[-which.max(nodup.rank[dups])]] <- NA
    }
  }

  ##remove the duplicated genes from results
  expr_mx.nodup <- expr_mx[!is.na(nodup.genes),]
  rownames(expr_mx.nodup) <- gene_symbols[!is.na(nodup.genes)]

  return(expr_mx.nodup)
}

