# This is the function for rhcoclust interaction network (internet)
#' @importFrom utils combn
#' @importFrom igraph graph
#' @importFrom igraph as_edgelist graph_from_edgelist layout_with_sugiyama V
#' @export
rhcoclust_internet <- function(data, CoClustObj, plot.internet = FALSE, cex.nodes = 0.7, edge.width = 1)
{
#data <- simu_data[1:30,1:30]
#data <- toxygates_data
# Apply rhcoclust to identify significant co-cluster of samples and their regulatory features
# CoClustObj <- rhcoclust(data, rk=3, ck=3, method.dist= "manhattan", method.hclust = "ward.D")
# plot_rhcoclust (CoClustObj, plot.coclust = TRUE, plot.ccim = FALSE)
# plot_rhcoclust (CoClustObj, plot.coclust = FALSE, plot.ccim = TRUE)

# co-cluster
RowClust <- paste0("R", 1 : length(CoClustObj$rowclust)) # for row
ColClust <- paste0("C", 1 : length(CoClustObj$colclust)) # for col

RowCoClust <- c()
len.1 <- 1 : length(ColClust)
for (i in 1 : length(RowClust)) {
  #SortIndx<-sort(1:ncol(CoClustObj$CoClsDtMat),decreasing = T)
  #CoClustObj$CoClsDtMat1 <- CoClustObj$CoClsDtMat[SortIndx,SortIndx]
  RowCoClust[[i]] <- sapply(len.1, function(len.1) paste0(RowClust[i], ColClust[len.1]))
}

# col
len.2 <- 1 : length(RowClust)
ColCoClust <- c()
for (i in 1 : length(ColClust)) {
  #SortIndx<-sort(1:ncol(CoClustObj$CoClsDtMat),decreasing = T)
  #CoClustObj$CoClsDtMat1 <- CoClustObj$CoClsDtMat[SortIndx,SortIndx]
  ColCoClust[[i]] <- sapply(len.2, function(len.2) paste0(ColClust[i], RowClust[len.2]))
}

# row
Combine_row <- c()
for (i in 1 : length(CoClustObj$rowclust)) {
  len_rowclust <- 1 : length(CoClustObj$rowclust[[i]])
  Combine_row[[i]] <- sapply(len_rowclust, function(len_rowclust) c(CoClustObj$rowclust[[i]][len_rowclust], RowCoClust[[i]]))
}
# cbind row
resultRow <- do.call("cbind", Combine_row)

CombRow <- c()
for (ii in 1 : ncol(resultRow)) {
  CombRow[[ii]] <- combn(resultRow[, ii], 2)[, 1 : length(CoClustObj$colclust)]
}
# graph for row
graph_row <- graph(unlist(CombRow), directed = T)


# Col
Combine_col <- c()
for (j in 1:length(CoClustObj$colclust)) {
  len_colclust <- 1 : length(CoClustObj$colclust[[j]])
  Combine_col[[j]] <- sapply(len_colclust,function(len_colclust) c(CoClustObj$colclust[[j]][len_colclust], ColCoClust[[j]]))
}

# cbind col
resultCol <- do.call("cbind", Combine_col)

CombCol <- c()
for (jj in 1:ncol(resultCol)) {
  CombCol[[jj]] <- combn(resultCol[, jj], 2)[, 1 : length(CoClustObj$rowclust)]
}
# graph for col
graph_col <- graph(unlist(CombCol), directed = T)


# calculate edges weight
# len_rowclust <- 1:length(CoClustObj$rowclust)
len_colclust <- 1 : length(CoClustObj$colclust)
data_weight <- c()
# data extracted
for (k in 1:length(CoClustObj$rowclust)) {
  data_weight[[k]] <- sapply(len_colclust, function(len_colclust) CoClustObj$CoClsDtMat[match(CoClustObj$rowclust[[k]], rownames(CoClustObj$CoClsDtMat)), match(CoClustObj$colclust[[len_colclust]], colnames(CoClustObj$CoClsDtMat))])
}

# mean for genes (rows)
mean_rows <- list()
len_colclust <- 1 : length(CoClustObj$colclust)
for (kk in 1:length(CoClustObj$rowclust)) {
  mean_rows[[kk]] <- sapply(len_colclust, function(len_colclust) apply(data_weight[[kk]][[len_colclust]], 1, mean)) # mean for genes (rows)
}
# rbind rows
Rbind_row <- do.call(rbind, mean_rows)

# mean for compounds (cols)
mean_cols <- list()
len_colclust <- 1 : length(CoClustObj$colclust)
for (l in 1 : length(CoClustObj$rowclust)) {
  mean_cols[[l]] <- sapply(len_colclust, function(len_colclust) apply(data_weight[[l]][[len_colclust]], 2, mean)) # mean for compounds (rows)
}

len_mean_cols <- 1 : length(mean_cols)
Rbind_col <- sapply(len_mean_cols, function(len_mean_cols) unlist(mean_cols[[len_mean_cols]]))

# combined edges and edges weight for row and col
# row
Edges_row <- as_edgelist(graph_row)
Edges_row_weight <- as.vector(t(Rbind_row))

# col
Edges_col <- as_edgelist(graph_col)
Edges_col[, 2] <- reversestring(Edges_col[,2],2)

Edges_col_weight <- as.vector(t(Rbind_col))

# all edges rbind
All_edges <- rbind(Edges_col, Edges_row)
# all weight combined
All_weight <- c(Edges_col_weight, Edges_row_weight)

# edge colors based on weight
color_edges <- rep("black", length(All_weight))
color_edges[which(All_weight > CoClustObj$UpContLimit)] <- "red"
color_edges[which(All_weight < CoClustObj$LowrContLimit)] <- "blue"

# combined all as a data frame
Graph_data <- data.frame(from=All_edges[, 1],
                         to=All_edges[, 2],
                         weight=All_weight,
                         color=color_edges)
UnRC <- unique(Graph_data$to)
# high value variables
HighValueVar <- sapply(UnRC, function(UnRC) as.character(Graph_data$from[which(Graph_data$to==UnRC & Graph_data$weight > CoClustObj$UpContLimit)]))
# low value variables
LowValueVar <- sapply(UnRC, function(UnRC) as.character(Graph_data$from[which(Graph_data$to==UnRC & Graph_data$weight < CoClustObj$LowrContLimit)]))

# graph objects
graph_obj <- graph_from_edgelist(
  as.matrix(Graph_data[, c('from', 'to')]))

# layout
layer <- rep(2, length(V(graph_obj)$name))
# layer for row (genes)
layer[match(noquote(rownames(data)), V(graph_obj)$name)] = 1
# layer for col
layer[match(noquote(colnames(data)),V(graph_obj)$name)] = 3
# layout
layout <- layout_with_sugiyama(graph_obj,layers=layer,maxiter = 1000,hgap = 1, vgap = 1)

layout$layout[which(layout$layout[, 2]==1), 1] <- seq(1, 3500, length.out=length(which(layout$layout[, 2]==1)))
layout$layout[which(layout$layout[, 2]==2), 1] <- seq(1, 3500, length.out=length(which(layout$layout[, 2]==2)))
layout$layout[which(layout$layout[, 2]==3), 1] <- seq(1, 3500, length.out=length(which(layout$layout[, 2]==3)))

if (plot.internet == TRUE)
{

#par(mar=c(4,4,4,4))
plot (graph_obj,
      #asp = 1,
      edge.color = color_edges,
      edge.width = edge.width,
      edge.width = 1,
      #vertex.color=V(graph.1)$color,
      #label.color=c(CoClustObj$colorsC,CoClustObj$colorsG),
      #color=c(CoClustObj$colorsC,CoClustObj$colorsG),
      #edge.curved=0,
      layout = cbind(layer,layout$layout[,1]),
      #layout=layer,
      vertex.shape = c("none","circle","none")[layer],
      vertex.size = c(20,13,25)[layer], #need to control by user
      #vertex.size2=40,
      #vertex.size2=strheight("I") * 2 * 100,
      #vertex.label=10,
      vertex.label.cex = cex.nodes,
      vertex.label.dist = 0,
      vertex.label.degree = 2,
      edge.arrow.size = 0,
      margin = -0.3)
  # Add a legend
  #legend("topleft",
  #       legend = c("Up-regulated", "No-regulated","Down-regulated"),
  #       col = c("red","black","blue"),
  #       bty = "n",
  #       pch = c(20,20,20),
  #       pt.cex = 2,
  #       cex = 1.2,
  #       x.intersp = 0.2,
  #       y.intersp = 0.4,
  #       text.col = "black",
  #       horiz = F ,
  #       inset = c(-0.2, 0))
}
return(list(HighValueVar = HighValueVar,
       LowValueVar = LowValueVar))
}

