# eset comparison script

orig_eset <-
new_eset <-

# compare sampleNames
smplnms <- identical(sampleNames(orig_eset), sampleNames(new_eset))
if(smplnms == F){
  diff_nms_o <- sampleNames(orig_eset)[ setdiff(sampleNames(orig_eset), sampleNames(new_eset))]
  diff_nms_n <- sampleNames(new_eset)[ setdiff(sampleNames(new_eset), sampleNames(orig_eset))]
  print(diff_nms_o)
  print(diff_nms_n)
}

# compare geneSymbols
genesyms <- identical(fData(orig_eset)$geneSymbol, fData(new_eset)$geneSymbol)
if(genesyms == F){
  print("GENE SYMS NOT EQUAL")
  print(paste0("length of orig ", length(fData(orig_eset)$geneSymbol)))
  print(paste0("length of new ", length(fData(new_eset)$geneSymbol)))
}

# pData rownames ... must match sampleNames
o_smpl <- sampleNames(orig_eset)
o_row <- rownames(pData(eset))
chk <- identical(o_smpl, o_row)
