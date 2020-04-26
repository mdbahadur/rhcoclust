#' @export
plot_rhcoclust <- function(CoClustObj){

# Plot results for gene (Row) and compound (column) co-cluster graph
par(mar=c(4,7,1,3))
par(mfrow=c(1,2))
# The reorganized transformed data matrix to generate co-cluster graph.
CoClustData <- CoClustObj$CoClsDtMat

# column and their ranked co-cluster mean in the second cluster.
GC_cls_MeanMat <- CoClustObj$NGC_cls_MeanMat

# Shape of points to generate individual control chart.
PcmQC <- CoClustObj$pchmark

# Colors to generate individual control chart.
ColorQC <- CoClustObj$color

# Colors of genes/row entity clusters to generate co-cluster graph
colors.genes <- CoClustObj$colorsG

# Colors of DCCs/column entity clusters to generate co-cluster graph
colors.dcc <- CoClustObj$colorsC

# Central Line of individual control chart to generate graph of control chart and to
# identify significant co-clusters.
CntrLine_QC <- CoClustObj$CentralLine

# Upper Control Limit to generate graph of control chart and to identify significant
# co-clusters.
UCL_QC <- CoClustObj$UpContLimit

# Lower Control Limit to generate graph of control chart and to identify significant
# co-clusters.
LCL_QC <- CoClustObj$LowrContLimit

# Display a color image
image(CoClustData,
      col = colorRampPalette(c("white","gray","black"),space = "rgb")(300),
      axes = FALSE,
      main="Co-cluster graph")

# Show row names for the image
mtext(text = rownames(CoClustData),
      side = 1,
      line = 0.3,
      at = seq(0,1,1/(nrow(CoClustData)-1)),
      las = 2,
      cex = 0.6,
      col = colors.genes)

# Show column names for the image
mtext(text = colnames(CoClustData),
      side = 2,
      line = 0.3,
      at = seq(0,1,1/(ncol(CoClustData)-1)),
      las = 2,
      cex = 0.6,
      col = colors.dcc)

# Show a legend strip for the color scale
image.plot(CoClustData,
           col = colorRampPalette(c("white","gray","black"),space = "rgb")(300),
           legend.only = TRUE,
           add = TRUE,
           xpd = TRUE,
           legend.shrink=1,
           horizontal = FALSE)

# Plot graph of QCC for identification of biomarker co-cluster
plot(x = GC_cls_MeanMat[,2],
     xlab = "Combination of Row and Column Cluster",
     ylab = "Co-cluster Average",
     main = "Graph for QCC",
     xaxt = 'n',
     pch = PcmQC,
     cex = 1.1,
     col = ColorQC,
     ylim = c(min(LCL_QC,min(GC_cls_MeanMat[,2]-5)),max(GC_cls_MeanMat[,2]+2)))

# Add straight lines to a plot
abline(h = c(UCL_QC,
             CntrLine_QC,
             LCL_QC),
       lty = 3,
       lwd = 3,
       col = "green")
# Show row names to a plot
mtext(text = GC_cls_MeanMat[,1],
      side = 1,
      line = 0.3,
      at = seq(1:length(GC_cls_MeanMat[,1])),
      las = 2,
      font = 2,
      cex = 1,
      col = "black")

# Add UCL,CL,LCL to a plot
text(c(length(GC_cls_MeanMat[,1])-0.5,length(GC_cls_MeanMat[,1])-0.5,
       length(GC_cls_MeanMat[,1])-0.5),
     c(UCL_QC, CntrLine_QC, LCL_QC),
     c("UCL","CL","LCL"),
     font = 4,
     cex = 1.5,
     col = "blue")
}
