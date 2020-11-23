# This is the function for generating plot for rhcoclust objects
#' @export
plot_rhcoclust <- function(CoClustObj, plot.cocluster = FALSE, plot.ccim = FALSE){

# Plot results for gene (row) and compound (column) co-cluster graph

# Transformed data matrix (rearranged by RHCOC algorithm) to generate co-cluster graph.
CoClustData <- CoClustObj$CoClsDtMat

# Row and column cluster and their co-cluster mean.
Coclust_MeanMat <- CoClustObj$Coclust_MeanMat

# Shape of points to generate control chart for individual measurement (CCIM).
PcmQC <- CoClustObj$pchmark

# Colors to generate control chart for individual measurement.
ColorQC <- CoClustObj$color

# Colors of row entity clusters to generate co-cluster graph.
colors.genes <- CoClustObj$colorsG

# Colors of column entity clusters to generate co-cluster graph
colors.dcc <- CoClustObj$colorsC

# Central line of control chart for individual measurement to generate graph and to
# identify significant co-clusters.
CntrLine_QC <- CoClustObj$CentralLine

# Upper control limit to generate graph of control chart for individual measurement and to identify significant
# co-clusters.
UCL_QC <- CoClustObj$UpContLimit

# Lower control limit to generate graph of control chart for individual measurement and to identify significant
# co-clusters.
LCL_QC <- CoClustObj$LowrContLimit

# Display a color image for co-cluster
if (plot.cocluster == TRUE)
{
image(CoClustData,
        col = colorRampPalette(c("white", "gray", "black"),space = "rgb")(300),
        axes = FALSE)
# Show row names for the image
mtext(text = rownames(CoClustData),
        side = 1,
        line = 0.3,
        at = seq(0, 1, 1 / (nrow(CoClustData) - 1)),
        las = 2,
        cex = 0.6,
        col = colors.genes)
# Show column names for the image
mtext(text = colnames(CoClustData),
        side = 2,
        line = 0.4,
        at = seq(0, 1, 1 / (ncol(CoClustData) - 1)),
        las = 2,
        cex = 0.4,
        col = colors.dcc)
# x-axis
mtext("Row clusters",
        side = 1,
        line = 3.5,
        cex = 1.5)
# y-axis
mtext("Column clusters",
        side = 2,
        line = 5.5,
        cex = 1.5)
# Title
title("Co-cluster graph",
        adj = 0.5,
        line = 1,
        cex.main = 1.5)
# Show a legend strip for the color scale
image.plot(CoClustData,
             col = colorRampPalette(c("white", "gray", "black"), space = "rgb")(300),
             legend.only = TRUE,
             add = TRUE,
             xpd = TRUE,
             legend.shrink = 1,
             legend.mar = 1.5,
             horizontal = FALSE)
}
# Plot graph of control chart for individual measurement for identification of biomarker co-clusters.
if (plot.ccim == TRUE)
{
#par(mar=c(6,10,3,6)) # Modify if needed
# mar order: bottom, left, top, and right
plot(x = Coclust_MeanMat[,2],
     xlab = "",
     ylab = "",
     xaxt = 'n',
     pch = PcmQC,
     cex = 1.1,
     col = ColorQC,
     ylim = c(min(LCL_QC, min(Coclust_MeanMat[, 2] - 5)), max(Coclust_MeanMat[, 2] + 2)))

# Add straight lines to a plot
abline(h = c(UCL_QC, CntrLine_QC, LCL_QC),
       lty = 3,
       lwd = 3,
       col = "green")
# Show row names to a plot
mtext(text = Coclust_MeanMat[,1],
      side = 1,
      line = 0.3,
      at = seq(1 : length(Coclust_MeanMat[, 1])),
      las = 2,
      font = 2,
      cex = 1,
      col = "black")
# x-axis
mtext("Combination of Row and Column Cluster",
      side = 1,
      line = 2.5,
      cex = 1.2)
# y-axis
mtext("Co-cluster Average",
      side = 2,
      line = 2.5,
      cex = 1.2)
# Title
title("Graph for CCIM",
      adj = 0.5,
      line = 1)

# Add UCL,CL,LCL to a plot
text(c(length(Coclust_MeanMat[, 1]) - 0.5, length(Coclust_MeanMat[, 1]) - 0.5,
         length(Coclust_MeanMat[, 1]) - 0.5),
     c(UCL_QC, CntrLine_QC, LCL_QC),
     c("UCL","CL","LCL"),
     font = 4,
     cex = 1.5,
     col = "blue")
}
}
