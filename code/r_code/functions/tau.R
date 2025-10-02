# tissue specificity tau metric index
#' @param x: a numeric vector (often gene expression across tissues, cell types, ...)
tau <- function(x){
  m = max(x)
  
  a = x / m
  
  N = length(x)
  
  S = sum(1 - a)
  
  y = S / (N - 1)
  
  return(y)
  
}