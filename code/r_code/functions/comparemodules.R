comparemodules <- function(ma, mb, f){
  fam_of = function(x) {f$gfam[f$id == x]} # lookup function to grab the family of a given gene
  gimme_fams = function(x) { # x is a set of gene names. is this function slow?
    l = c()
    for (i in x){
      l = c(l, fam_of(i))
    }
    l = unique(l)
    return(l)
  }
  gene_POP <- length(unique(f$gfam)) # population size
  
  #' Create a matrix to store stats
  fon = data.frame(
    module_a = "none", 
    module_b = "none", 
    success_in_samples = 0, 
    sample_size = 0, 
    success_in_pop = 0, 
    gene_POP = gene_POP, 
    hypgeom_pval = numeric(1), 
    hypgeom_log = numeric(1), 
    binom_pval =  numeric(1), 
    gfams_common = "none",
    gfams_excl_a = "none",
    gfams_excl_b = "none"
  )
  
  #' Transform into lists for practicality
  ma_list <- split(ma$id, ma$module) 
  mb_list <- split(mb$id, mb$module)
  
  #' Create a matrix to store results (-logpval or similar)
  PVS <- data.frame() # create matrix to store result pval
  PVS_binom <- data.frame()
  
  
  #' Compute the upper tail of the hypergeometric 
  #' distribution (survival function) for each pair of modules
  #'  as a metric of enrichment #comment from panos @ skarmetalab
  for (i in 1:length(ma_list)){ #need to dramatically optimise speed
    
    a_modulei_name <- names(ma_list[i])
    a_modulei_genes <- ma_list[[i]]
    a_fams <- gimme_fams(a_modulei_genes) #is this slow?
    
    for (j in 1:length(mb_list)){
      
      b_modulej_name <- names(mb_list[j])
      b_modulej_genes <- mb_list[[j]]
      b_fams <- gimme_fams(b_modulej_genes) #is this slow?
      
      sample_size <- length(a_fams)
      success_in_pop <- length(b_fams)
      common_fams <- paste(
        a_fams[which(a_fams %in% b_fams)], collapse = ", "
      )
      exclusive_fams_a <- paste(
        a_fams[which(!(a_fams %in% b_fams))], collapse = ", "
      )
      exclusive_fams_b <- paste(
        b_fams[which(!(b_fams %in% a_fams))], collapse = ", "
      )
      success_in_sample <- length(which(a_fams %in% b_fams))
      
      if (success_in_sample > 0) {
        hypg <- phyper( # HYPGEOMTEST
          # from stackoverflow/questions/8382806/hypergeometric-test-phyper
          q = success_in_sample - 1, # no. of success balls drawn from urn
          m = success_in_pop, # no. of success balls in the urn
          n = gene_POP-success_in_pop, # no. of non-success in the urn
          k = sample_size, # no. of balls drawn from urn
          lower.tail = FALSE
        )
        binom <- binom.test(
          x = success_in_sample, 
          n = sample_size, 
          p = success_in_pop / gene_POP, 
        )$p.value
        
        fon <- rbind(
          fon, 
          c(
            a_modulei_name, b_modulej_name, 
            success_in_sample, sample_size, 
            success_in_pop, gene_POP, as.numeric(hypg), 
            -log(hypg), as.numeric(binom), common_fams, # add here which are the names of the gfams enriched.
            exclusive_fams_a, # exclusive fams in a
            exclusive_fams_b # exclusive fams in b
          )
        )
      } else {
        hypg <- 1
        binom <- 1
      }
      
      PVS[i, j] <- hypg
      PVS_binom[i, j] <- binom
      
    }
  }
  
  # Tidy up of data 
  rownames(PVS) <- names(ma_list)
  colnames(PVS) <- names(mb_list)
  rownames(PVS_binom) <- names(ma_list)
  colnames(PVS_binom) <- names(mb_list)
  fon <- fon[!(fon$module_a == "none"), ]
  fon$hypgeom_pval <- as.numeric(fon$hypgeom_pval)
  fon$hypgeom_log <- as.numeric(fon$hypgeom_log)
  fon$binom_pval <- as.numeric(fon$binom_pval)
  
  loghypg <- -log(PVS)
  loghypg[loghypg > 30] <- 30
  
  logbinom <- -log(PVS_binom)
  logbinom[logbinom > 30] <- 30
  
  
  res <- list(
    stats = fon, 
    hypgeom = PVS, 
    binom = PVS_binom, 
    loghypg = loghypg, 
    logbinom = logbinom
  )
  
  return(res)
}
