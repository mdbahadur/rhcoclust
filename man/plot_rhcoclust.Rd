\name{plot_rhcoclust}
\alias{plot_rhcoclust}

\title{
Plot of the 'rhcoclust' objects
}
\description{
This function used for two plots from output of rhcoclust (i) plot results for gene (row) and compound (column) co-cluster graph, and (ii) plot graph of QCC for identification of biomarker co-cluster.
}
\usage{
plot_rhcoclust(CoClustObj, plot.cocluster = FALSE, plot.ccim = FALSE)
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
# par(mar=c(6, 10, 3, 6)) # Modify if needed
# mar order: bottom, left, top, and right
plot_rhcoclust (CoClustObj, plot.coclust = TRUE, plot.ccim = FALSE)

# Plot CCIM
# use dev.off() to avoid the figure margin from previous plot
plot_rhcoclust (CoClustObj, plot.coclust = FALSE, plot.ccim = TRUE)

}
