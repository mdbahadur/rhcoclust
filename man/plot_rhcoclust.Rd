\name{plot_rhcoclust}
\alias{plot_rhcoclust}

\title{
Plot of the 'rhcoclust' objects
}
\description{
This function used for two plots from output of rhcoclust (i) plot results for gene (row) and compound (column) co-cluster graph, and (ii) plot graph of QCC for identification of biomarker co-cluster.
}
\usage{
plot_rhcoclust(CoClustObj, plot.cocluster = FALSE, plot.SCC = FALSE,
cex.xaxis = 0.7, cex.yaxis = 0.5)
}

\arguments{
  \item{CoClustObj}{
Output objects from rhcoclust
}
\item{plot.cocluster}{
To set no plotting as the default for cocluster.
}
\item{plot.SCC}{
To set no plotting as the default for SCC.
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
data("FCGE_Data_GMP")
data("FCGE_Data_PPARs")
# Load predefined real data
# Real data use: data <- FCGE_Data_PPARs
# Real data use: data <- FCGE_Data_GMP

# Load predefined simulated data
data("simu_data")

# simulated data
data <- simu_data

# Apply rhcoclust to identify significant co-cluster of samples and their regulatory features
CoClustObj <- rhcoclust(data, rk=4, ck=3, method.dist = "manhattan", method.hclust = "ward.D")
# For real data either FCGE_Data_PPARs or FCGE_Data_GMP
#CoClustObj <- rhcoclust(data, rk=3, ck=3, method.dist = "manhattan", method.hclust = "ward.D")

# Plot co-cluster
# Please use par(mar=c(6, 10, 3, 6))  or modify if needed for best fit of the graph
# mar order: bottom, left, top, and right
plot_rhcoclust (CoClustObj, plot.coclust = TRUE, plot.SCC = FALSE,
cex.xaxis = 0.7, cex.yaxis = 0.5)

# Plot SCC
# Please use dev.off() to avoid the figure margin from previous plot
plot_rhcoclust (CoClustObj, plot.coclust = FALSE, plot.SCC = TRUE)
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
