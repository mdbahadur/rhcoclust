\name{simulate_data}
\alias{simulate_data}
\title{
Simulate data for robust hirarachical clustering to identify significant co-cluster
}
\description{
We generate fold change gene expression (FCGE) data according to the characteristics of toxicogenomic data.
}
\usage{
simulate_data(no.gene, no.dcc)
}
\arguments{
  \item{no.gene}{
Number of genes in the simulated data.
}
  \item{no.dcc}{
Numbe of doses of chemical compounds (dcc) in the simulated data.
}
}
\details{
There are four gene groups and three DCCs groups in the simulated dataset. The gene group 1, 2, 3 and 4 consists the genes G1-G10, G11-G20, G21-G30 and G31-G50 respectively and DCCs group 1, 2 and 3 consists the DCCs C1_High-C5_High- C1_Middle-C5_Middle, C6_High-C10_High- C6_Middle-C10_Middle and C1_Low-C12_Low- C11_Middle- C12_Middle- C11_High- C12_High respectively. Where, G stands for gene and C stands for chemical compound arranged in the row and column of the simulated data matrix respectively. The error term N(0,0.35) from normal distribution with mean 0 and variance 0.35 is added to each element of the simulated dataset. In the simulated dataset the gene group-1 is up-regulated by the DCCs group-1, gene goup-2 is up and down-regulated by the DCCs group-2 and 1 respectively. The gene group-3 is down-regulated by the DCs group-2. The gene group-4 is not regulated by any of the DCCs groups and DCCs group-3 does not influence any of the genes in the dataset.
}
\value{
A \link{list} of object that containing the following:

SimData: A simulated data matrix as generated.

SimDataRnd: Randomly distributed row and column entity of simulated data.

GCmat: Transformed simulated data.
}

\author{
Md. Bahadur Badsha <mbbadshar@gmail.com>
}

\examples{

# Load library
library(rhcoclust)

# Number of genes in the simulated data.
no.gene <- 50

# Numbe of doses of chemical compounds (dcc) in the simulated data.
no.dcc <- 12

SimulteData <- simulate_data(no.gene,no.dcc)

# A simulated data matrix as generated.
SimulteData$SimData

# A randomly distributed row and column entity of simulated data.
SimulteData$SimDataRnd

# A transformed simulated data.
SimulteData$GCmat
}

