#' @param l = list with overlapping sets of events
make_yesno <-
  function(l){
    N = length(l)
    
    events <- unique(unname(unlist(l)))
    n <- length(events)
    
    m <- matrix(0,nrow = n, ncol = N, dimnames = list(events, names(l)))
    
    for(i in 1:N){
      m[,i] <- as.integer(sapply(rownames(m),function(x){x %in% l[[i]]}))
    }
    
    m # a matrix of presence/absence of each event in each set
  }