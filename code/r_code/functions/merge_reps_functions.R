FUN_by_rep <-
  function(x, idents, FUNCTION){
    
    f <- function(a,FUN){
      FUNCTION( x[ names(x) %in% names(idents[idents==a]) ] )
    }
    
    setNames(sapply(unique(idents), FUN = f),unique(idents))
  }

merge_reps <-
  function(x,idents,FUNCTION){
    
    t(apply(
      X = x, MARGIN = 1, FUN = FUN_by_rep,
      idents = idents,
      FUNCTION = FUNCTION 
    ))
    
  }
