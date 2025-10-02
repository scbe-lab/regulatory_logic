#' Returns a T/F vector of same length as x, removing 2nd, 3rd, n-th instances
#' of repeated values (i.e. keeping the 1st instance of any repeated value).
#' IF USED IN ROWNAMES/COLNAMES OF MATRICES, ONLY USE WHEN THE ROW/COLUMN VALUES
#' ARE 100% THE SAME (i.e. which of the two does not matter)
remove_duplicates <- function(x,seed = 4343){
  if(any(is.na(x))) stop("E: there are NAs in input, please check.")
  
  d <- duplicated(x)
  
  y <- !d
  
  for(i in unique(x[d])){ y[y == i][1] <- TRUE }
  
  return(y)
  
}