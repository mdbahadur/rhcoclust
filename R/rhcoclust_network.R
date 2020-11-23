#' @export
rhcoclust_network <- function(CoClustObj, scale.threshold = 10)
{
  # Pre-process row and column names
  Row.names <- paste0("R", CoClustObj$NG_Cocls) # R indicates row
  Col.names <- paste0("C", CoClustObj$NC_Cocls) # C indicates column
  len <- 1 : length(Col.names)
  Combine.Row.Col <- sapply(len, function(len) c(Row.names[len], Col.names[len]))

  # Create networks
  # library(igraph)
  graph.obj <- graph(Combine.Row.Col, directed = F)
  color.edge <- rep("black", length(CoClustObj$Coclust_MeanMat$GC_CoMeanR))
  color.edge[which(CoClustObj$Coclust_MeanMat$GC_CoMeanR > CoClustObj$UpContLimit)] <- "red"
  color.edge[which(CoClustObj$Coclust_MeanMat$GC_CoMeanR < CoClustObj$LowrContLimit)] <- "blue"

  # Visualization of clustering network plot
  plot(graph.obj,
       #edge.color=rep(c("red","pink"),5),           # Edge color
       edge.color = color.edge,
       #vertex.color = rgb(0.8,0.4,0.3,0.8),
       #edge.width=(CoClustObj$Coclust_MeanMat$GC_CoMeanR/max(CoClustObj$Coclust_MeanMat$GC_CoMeanR))+0.5,
       edge.width = CoClustObj$Coclust_MeanMat$GC_CoMeanR / scale.threshold,# Edge width, defaults scale.threshold 10
       edge.arrow.size = 1,                           # Arrow size, defaults to 1
       edge.arrow.width = 1,                          # Arrow width, defaults to 1
       edge.lty = c("solid"),                           # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
       #edge.curved=c(rep(0,5), rep(1,5))            # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
       main="Clustering Network Plot"
  )

}
