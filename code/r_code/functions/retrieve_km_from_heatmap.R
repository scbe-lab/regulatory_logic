# from: https://github.com/jokergoo/ComplexHeatmap/issues/136#issuecomment-1452693518
retrieve_km_from_heatmap <-
  function(HM, seed = 1234){
    # if(any(dim(M) != dim(HM@matrix))) stop("E: heatmap and matrix provided are not equal dim. Are they the same?")
    if(!(class(HM)[1] %in% c("Heatmap","HeatmapList"))) stop("E: heatmap class not recognised. Check the class of object provided")
    
    set.seed(seed)
    rcl.list <- row_order(HM)  #Extract clusters (output is a list)
    
    if(class(HM) == "HeatmapList"){
      rn <- rownames(HM@ht_list$expression@matrix)
    } else{
      rn <- rownames(HM@matrix)
    }
    
    clu_df <-
      stack(lapply(rcl.list,function(x){rn[x]}))
    
    colnames(clu_df) <- c("id","cluster")
    clu_df$cluster <- paste0("cluster",as.character(clu_df$cluster))
    
    clu_df$cluster <- factor(clu_df$cluster, levels = sort(unique(clu_df$cluster)))
    
    # clu_df <- clu_df[order(clu_df$id,clu_df$cluster),]
    
    return(clu_df)
  }
