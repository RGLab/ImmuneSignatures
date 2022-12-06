# inverse normal (rank) transformation
inv_normal_transform <- function(x){
  r = rank(x, na.last = "keep")
  p = r / (sum(!is.na(r)) + 1)
  return(qnorm(p))
}
