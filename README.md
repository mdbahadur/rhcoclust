# rhcoclust

rhcoclust function performs robust hierarchical co-clustering between row and column entities of a data matrix by reducing the influence of outlying observations. It can be used to explore biomarker genes those are 
divided into upregulatory and downregulatory groups by the influence of different chemical compounds groups more accurately. It can also provide the statistical significance of the identified co-clusters.  See Hasan et al., (2020).
Development of the R package rhcoclust is at https://github.com/mdbahadur/rhcoclust, and official releases are available on CRAN: https://cran.r-project.org/web/packages/rhcoclust/index.html.

References:
  
1. Mohammad Nazmol Hasan, Md. Bahadur Badsha, Md. Nurul Haque Mollah: Robust Hierarchical Co-clustering to Explore Toxicogenomic Biomarkers and Their Regulatory Doses of Chemical Compounds. bioRxiv, 2020.05.13.094946 (2020).

## Installation

### 1. Installation of the most recent version from GitHub.

First install the R package devtools available on CRAN, if it is not already installed. This package provides function `install_github()` that enables installing packages directly from github with the following command.

Invoke R and then type with the following command:
  ```
R> install.packages ("devtools")
R> library (devtools)
# install R packages that rhcoclust depends on before running the next line 
# see details below
R> install_github ("mdbahadur/rhcoclust")
```
rhcoclust depends on several R packages from CRAN and from Bioconductor.  It is likely that some of these packages are not installed on your computer.  If the R package is available on CRAN, you may use the following command line for installation (change _packagename_ to the name of the package to be installed, e.g, spam, maps, etc.) before running function `install_github`:
  ```
R> install.packages ("packagename")
```

The following Bioconductor packages also need to be installed before running function `install_github` or `install.packages`:
  ```
R> if (!requireNamespace ("BiocManager", quietly = TRUE))
  install.packages ("BiocManager")
R> BiocManager::install ('spam')
R> BiocManager::install ('maps')
```
### 2. Installation from the source of a released package.

Download the package source rhcoclust_xxx.tar.gz.  

In Terminal, navigate to the directory where the package is stored, and run the following command line:
  ```bash
$ R CMD INSTALL rhcoclust_xxx.tar.gz
```
Again, you may need to first install the Bioconductor packages that rhcoclust depends on using the instructions above.
Alternatively, you may also run the following command line in R, after changing the working directory to where rhcoclust_xxx.tar.gz is stored on your computer:
  ```
R> install.packages ("rhcoclust_xxx.tar.gz", repos = NULL, type="source")
```
### 3. Installation from CRAN.

Official releases are available on CRAN.  To install,
```
R> install.packages ("rhcoclust")
```
## Using rhcoclust
After installation, load the rhcoclust package into R:
  ```
R> library (rhcoclust)
```
Bring up the documentation of the package:
  ```
R> library (help=rhcoclust)
```


