feat <- function(x, scdata){
  FeaturePlot(
    scdata,
    features = x,
    cols = c("#DAE7F2","#531ccb"),
    order = TRUE,
    pt.size = 1
  )+ggtitle(x)+NoAxes()
}
