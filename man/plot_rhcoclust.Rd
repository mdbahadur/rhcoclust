\name{plot_rhcoclust}
\alias{plot_rhcoclust}

\title{
Plot of the 'rhcoclust' objects
}
\description{
This function used for two plots from output of rhcoclust (i) plot results for gene (row) and compound (column) co-cluster graph, and (ii) plot graph of QCC for identification of biomarker co-cluster.
}
\usage{
plot_rhcoclust(CoClustObj, plot.cocluster = FALSE, plot.ccim = FALSE,
cex.xaxis = 0.7, cex.yaxis = 0.5)
}

\arguments{
  \item{CoClustObj}{
Output objects from rhcoclust
}
\item{plot.cocluster}{
To set no plotting as the default for cocluster.
}
\item{plot.ccim}{
To set no plotting as the default for icc.
}
\item{cex.xaxis}{
A numerical value giving to control/annotation text size in x-axis. Default is 0.7.
}
\item{cex.yaxis}{
A numerical value giving to control/annotation text size in y-axis. Default is 0.5.
}


}

\value{
Plots
}

\author{
Md. Bahadur Badsha <mbbadshar@gmail.com>
}

\seealso{
\link{rhcoclust} for generating a graph objects for clustering network
}

\examples{
# Load necessary library
library(rhcoclust)
library(fields)

# Load real data
data("toxygates_data")

# Load predefined simulated data
data("simu_data")

# Data analysis using real data by rhcoclust. Note, for analysis of predefined
# simulated data use data <- simu_data
data <- toxygates_data

# Apply rhcoclust to identify significant co-cluster of samples and their regulatory features
CoClustObj <- rhcoclust(data, rk = 3, ck = 3, method.dist = "manhattan", method.hclust = "ward.D")

# Plot co-cluster
# Please use par(mar=c(6, 10, 3, 6))  or modify if needed for best fit of the graph
# mar order: bottom, left, top, and right
plot_rhcoclust (CoClustObj, plot.coclust = TRUE, plot.ccim = FALSE,
cex.xaxis = 0.7, cex.yaxis = 0.5)

# Plot CCIM
# Please use dev.off() to avoid the figure margin from previous plot
plot_rhcoclust (CoClustObj, plot.coclust = FALSE, plot.ccim = TRUE)
# Please add legend with change or add any parameters if needed.
legend("topleft",
       legend = c("Upper-significant", "Insignificant","Down-significant"),
       col = c("red","black","blue"),
       bty = "n",
       pch = c(20,20,20),
       pt.cex = 2,
       cex = 1.2,
       x.intersp = 0.2,
       y.intersp = 0.4,
       text.col = "black",
       horiz = FALSE ,
       inset = c(0.3, -0.08))
}
