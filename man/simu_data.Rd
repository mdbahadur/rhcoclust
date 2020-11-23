\name{simu_data}
\alias{simu_data}

\title{
A predefined simulated data for 'rhcoclust' package
}
\description{
This is the  predefined simulated data matrix that used as an example in 'rhcoclust' package.
}

\details{
The column of the simulated data matrix represents doses of chemical compounds (DCCs) and row represents genes. Each element in the dataset represents fold change gene expression data.
}

\value{
Matrix
}

\author{
Md Bahadur Badsha <mbbadshar@gmail.com>
}

\seealso{
\link{simulate_data}
}

\examples{
library(rhcoclust)

# Load data
data("simu_data")
}
