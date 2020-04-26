\name{rhcoclust}
\alias{rhcoclust}
\title{The function for co-clustering sample and feature to explore significant samples and their regulatory features}
\usage{
rhcoclust(data, rk, ck)
}
\arguments{
\item{data}{A data matrix containing data having the characteristics of interval and ratio level of measurement or continuous data}

\item{rk}{Number of clusters in the row entities of the data matrix}

\item{ck}{Number of clusters in the column entities of the data matrix}
}
\value{
A \link{list} of object that containing the following:

NGC_cls_MeanMat: A data frame containing combination of row and column cluster number in the first column and their ranked co-cluster mean in the second cluster.

CoClsDtMat: The reorganized transformed data matrix to generate co-cluster graph.

Co_Gcls: The gene/row entity clusters.

Co_Ccls: The column entity clusters.

colorsG: Colors of genes/row entity clusters to generate co-cluster graph.

colorsC: Colors of DCCs/column entity clusters to generate co-cluster graph.

CentralLine: Central Line of individual control chart to generate graph of control chart and to identify significant co-clusters.

UpContLimit: Upper Control Limit to generate graph of control chart and to identify significant co-clusters.

LowrContLimit: Lower Control Limit to generate graph of control chart and to identify significant co-clusters.

color: Colors to generate individual control chart.

pchmark: Shape of points to generate individual control chart.
}
\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}
\description{
Toxicogenomic studies require co-clustering to identify biomarker genes for the assessment of chemical toxicity from gene expression levels. It is also essential in the drug discovery experiments. However, gene expression datasets are often contaminated by outliers due to several steps involve in the data generating process. This package performs robust hierarchical co-clustering between row and column entities of a data matrix by reducing the influence of outlying observations. It can be used to explore biomarker genes those are divided into upregulatory and downregulatory groups by the influence of different chemical compounds groups more accurately. It can also provide the statistical significance of the identified co-clusters.
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
CoClustObj <- rhcoclust(GCmat, rk=3, ck=3)

# A data frame containing combination of row and column cluster number in the first
# column and their ranked co-cluster mean in the second cluster.
GC_cls_MeanMat <- CoClustObj$NGC_cls_MeanMat

# The reorganized transformed data matrix to generate co-cluster graph.
CoClsDtMat <- CoClustObj$CoClsDtMat

# The gene/row entity clusters.
Co_Gcls <- CoClustObj$Co_Gcls

# The column entity clusters.
Co_Ccls <- CoClustObj$Co_Ccls

# Colors of genes/row entity clusters to generate co-cluster graph
colorsG <- CoClustObj$colorsG

# Colors of DCCs/column entity clusters to generate co-cluster graph
colorsC <- CoClustObj$colorsC

# Central Line of individual control chart to generate graph of control chart and to
# identify significant co-clusters.
CntrLine_QC <- CoClustObj$CentralLine

# Upper Control Limit to generate graph of control chart and to identify significant
# co-clusters.
UCL_QC <- CoClustObj$UpContLimit

# Lower Control Limit to generate graph of control chart and to identify significant
# co-clusters.
LCL_QC <- CoClustObj$LowrContLimit

# Colors to generate individual control chart.
ColorQC <- CoClustObj$color

# Shape of points to generate individual control chart.
PcmQC <- CoClustObj$pchmark

# Plot the results
plot_rhcoclust (CoClustObj)

}