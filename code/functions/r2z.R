r2z <- function(x){
  x <- 0.5 * log((1+x)/(1-x)); # Fisher-z transform
  return(x)
}