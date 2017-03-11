# Evan Henrich
# January 2017
# Fred Hutchinson Cancer Research Center
# Gottardo Lab
# ehenrich@fredhutch.org

#--------GEN NOTES-------------------------------------------------------
# This script maps the sampleNames from an exprs(eset) matrix to the
# corresponding ImmuneSpace subject IDs using schema tables from ImmuneSpace
# or mapping tables sent by the Yale collaborators in HIPC.

#-------HELPER METHOD-------------------------------------
hash_it <- function(tmp_hash, key){
  present <- has.key(key, tmp_hash)
  ifelse(present == T, val <- tmp_hash[[as.character(key)]], NA)
}

#--------MAIN METHOD--------------------------------------
#' Function to generate HAI data table from ImmuneSpace Connection
#'
#' @param keys keys
#' @param values Values
#' @param vec_2_hash Vector to be hashed
#' @export
hashmap <- function(keys, values, vec_2_hash){
  tmp_hash <- hash(keys, values)
  mapped_res <- sapply(vec_2_hash, FUN = function(x){hash_it(tmp_hash, x) })
}
