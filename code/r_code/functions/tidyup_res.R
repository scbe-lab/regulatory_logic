tidyup_res <- function(d,filter_by = "padj", p_thres = 0.1, lfc_thres = 0, as_data_frame = FALSE){
  if(!(filter_by %in% c("pvalue","padj"))) stop("E: filter_by not recognised. Only 'pvalue' and 'padj' allowed.")
  d <- d[complete.cases(d),]
  f <- d[[filter_by]] < p_thres & abs(d$log2FoldChange) > lfc_thres
  d$DEG <- ifelse(f, TRUE, FALSE)
  d$updown <- "none"
  d$updown[f & d$log2FoldChange < 0 ] <- "down"
  d$updown[f & d$log2FoldChange > 0 ] <- "up"
  
  if(as_data_frame) d <- as.data.frame(d)
  
  return(d)
}