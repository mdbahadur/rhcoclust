\name{toxygates_data}
\alias{toxygates_data}

\title{
A real dataset for 'rhcoclust' package
}
\description{
This is the real data matrix that used as an example in 'rhcoclust' package.
}

\details{
The real data matrix collected from Nyström-Persson et al., 2013 (https://toxygates.nibiohn.go.jp/toxygates/#columns).The column of the data matrix represents doses of chemical compounds (DCCs) and row represents genes. Each element in the dataset represents fold change gene expression data.
}

\value{
Matrix
}

\author{
Md. Bahadur Badsha <mbbadshar@gmail.com>
}
\references{
1. Nyström-Persson et al., (2013). Toxygates: interactive toxicity analysis on a hybrid microarray and linked data platform. Bioinformatics, Volume 29, Issue 23,Pages 3080–3086.
}
\examples{
# Load library
library(rhcoclust)

# Load real data
data("toxygates_data")
}

