#' Find threshold value separating the first two modes of a distribution
#' inspired from: https://stackoverflow.com/a/77115370
#' hard_threshold parameter checks if the value foundexcludes more than 75% of
#' the data. If so, just keep the .75 quantile as result
antimode <- function(x, hard_threshold = FALSE, q = 0.75) {
  options(warn=-1)
  require(multimode)
  bw <- bw.nrd0(x) # bandwith
  n <- nmodes(x, bw = bw) # Calculate number of modes
  loc <- locmodes(x, mod0 = n, display = FALSE) # Calculate location of these modes
  am <- loc$locations[seq(from = 2, by = 2, length.out = n-1)] #antimode: even values
  y <- am[am > median(x)][1] # first antimode above the median values
  
  if(is.na(y)) y = quantile(x,q)
  
  if(hard_threshold == TRUE){
    if(y > quantile(x,q)) y <- quantile(x,q) 
  }
  
  return(y)
  
}