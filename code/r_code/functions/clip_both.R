
clip_both <- function(x, Q = 0.95){
  stopifnot( is.numeric(Q) & Q >= 0 & Q <= 1 )
  y = 
    clip_q(
      clip_q(x, q = Q,method = "up"),
      q = 1-Q,method = "down"
      )
  
  return(y)

}