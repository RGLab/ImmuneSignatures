discretize <- function(x, q) {
  xq = quantile(x, c(q, 1-q), na.rm=T)
  xd = ifelse(is.na(x),NA,1)
  xd[x<=xq[1]] = 0
  xd[x>=xq[2]] = 2
  return(xd)
}
