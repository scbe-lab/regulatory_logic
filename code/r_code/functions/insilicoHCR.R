insilicoHCR <- function(
    scdata, feat_a, feat_b,
    col_a = c("#000000","#ff0000"),
    col_b = c("#000000","#00ff00"),
    clip = FALSE,q = 0.9999,
    col_fun = "sum_cols", col_method = "discrete"
    ){
  
  if( !(feat_a %in% rownames(scdata)) | !(feat_b %in% rownames(scdata)) ) {
    stop("error: at least one of the features specified is not in this dataset")
  }
  
  if(!(col_fun %in% c("average_cols", "sum_cols"))){
    stop("error: function for merging colors not recognised.")
  }
  
  if (!(col_method %in% c("gradient","discrete"))){
    stop("error: method of colouring cells not recognised.")
  }
  
  # values feature a
  a <- scdata@assays$RNA@counts[feat_a,]
  
  if(clip == TRUE){
    q_a <- quantile(a[a>0],q)
    a[a>q_a] <- q_a
  }
  
  a <- relativise(a)
  
  # color map and colouring feature a
  cmap_a <- colorRamp2(breaks = c(min(a),max(a)), col = col_a, space = "RGB")
  a_c <- cmap_a(a)
  
  # values feature b
  b <- scdata@assays$RNA@counts[feat_b,]
  
  if(clip == TRUE){
    q_b <- quantile(b[b>0],q)
    b[b>q_b] <- q_b
  }
  
  b <- relativise(b)
  
  # color map and colouring feature a
  cmap_b <- colorRamp2(breaks = c(min(b),max(b)), col = col_b, space = "RGB")
  b_c <- cmap_b(b)
  
  # merge colours
  col_fun <- get(col_fun)
  
  d <- data.frame(a_c = a_c, b_c = b_c)
  col_merge <- apply(d,1, col_fun)
  
  order_dots <- order(a+b)
  
  d <-
    cbind(
      d,
      data.frame(
        cell = colnames(scdata),
        a = a,
        b = b,
        comb_expr = a+b,
        x = scdata@reductions$umap@cell.embeddings[,1],
        y = scdata@reductions$umap@cell.embeddings[,2],
        gradient = col_merge
      )
    )
  
  d_$discrete <- "black"
  d_$discrete[d_$a > 0 ] <- alpha(col_a[2],.5)
  d_$discrete[d_$b > 0 ] <- alpha(col_b[2],.5)
  d_$discrete[d_$a > 0 & d_$b > 0 ] <- sum_cols(c(col_a[2],col_b[2]))
  
  d_$plotting_col <- d_[,col_method]
  
  y <- 
    d_%>%
    arrange(comb_expr) %>%
    ggplot(mapping = aes(x = x, y = y, col = plotting_col, fill = plotting_col))+
    geom_point(size = 1)+
    scale_color_identity()+
    theme_bw()+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    ggtitle(paste0("a:",feat_a," + ","b:",feat_b))
  
  return(y)
}
