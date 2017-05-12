plot_size_bubbles <- function(x,y, col=NA, f=NA) {
  dd = data.frame(x,y, col, f) %>%
    group_by(x,y,col,f) %>%
    summarise(n=n()) %>%
    ungroup()
  
  if (all(is.na(dd$col))) {
    p = ggplot(dd, aes(x,y))
  } else {
    p = ggplot(dd, aes(x,y, colour=col))
  }
  
  if (!all(is.na(dd$f))) {
    p = p + facet_wrap(~f, nrow=1)
  }
  p + geom_point(aes(size=n))
}
