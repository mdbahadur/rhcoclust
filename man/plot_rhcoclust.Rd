\name{plot_rhcoclust}
\alias{plot_rhcoclust}

\title{
Plot of the 'rhcoclust' objects
}
\description{
This function used for two plots from output of rhcoclust (i) plot results for gene (row) and compound (column) co-cluster graph, and (ii) plot graph of QCC for identification of biomarker co-cluster.
}
\usage{
plot_rhcoclust(CoClustObj)
}

\arguments{
  \item{CoClustObj}{
Output objects from rhcoclust
}
}

\value{
Plots
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
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
# simulated data use GCmat <- simu_data
GCmat <- toxygates_data

# Apply rhcoclust to identify significant co-cluster of samples and their regulatory features
CoClustObj <- rhcoclust(GCmat, rk=3, ck=3, method.dist="manhattan",method.hclust="ward.D")

# Plot the results
plot_rhcoclust (CoClustObj)
}
