#' @importFrom stats cutree dist hclust rnorm

#' @export
rhcoclust <- function(data, rk, ck, method.dist = "manhattan", method.hclust = "ward.D")
{
  # Data transformation (expressed in %) using logistic function
  dataExpr <- 100*( 1 / ( 1 + exp(-data)))
  dG <- dist((dataExpr), method = method.dist)
  dC <- dist(t(dataExpr), method = method.dist)
  # Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it.
  HCGene <- hclust(dG, method = method.hclust)
  HCComp <- hclust(dC, method = method.hclust)
  # rk row clusters estimated observing dendrogram
  HCclsMatG <- as.matrix(cutree(HCGene, rk))
  # ck column clusters estimated observing dendrogram
  HCclsMatC <- as.matrix(cutree(HCComp, ck))
  # -----------------------------------------------------------------------------------

  # Loop to get list of row and column clusters
  HCclsG <- list()
  HCclsC <- list()
  for (p in 1:rk)
  {
    HCclsG[[p]] <- rownames(HCclsMatG)[which(HCclsMatG[, 1] == p)]
  }
  for (q in 1:ck)
  {
    HCclsC[[q]] <- rownames(HCclsMatC)[which(HCclsMatC[, 1] == q)]
  }
  # -----------------------------------------------------------------------------------

  # Loop to get co-cluster means and corresponding row and column cluster number
  Gcls <- NULL
  Ccls <- NULL
  Co_ClsMeanExpr <- NULL
  for (g in 1 : length(HCclsG))
  {
    for (c in 1 : length(HCclsC))
    {
      Gcls <- c(Gcls, g)
      Ccls <- c(Ccls, c)
      Co_ClsMean <- mean(abs(dataExpr[HCclsG[[g]], HCclsC[[c]]]))
      Co_ClsMeanExpr <- c(Co_ClsMeanExpr, Co_ClsMean)
    }
  }
  # -----------------------------------------------------------------------------------

  # Ordering co-cluster mean and their corresponding row and column clusters
  RankExprMean <- Co_ClsMeanExpr[order(-Co_ClsMeanExpr)]
  GclsR <- Gcls[order(-Co_ClsMeanExpr)]
  CclsR <- Ccls[order(-Co_ClsMeanExpr)]
  # -----------------------------------------------------------------------------------

  # Loop to get row and column clusters according to ranked co-cluster mean
  CoclG <- unique(GclsR)
  CoclC <- unique(CclsR)
  rowclust <- list()
  colclust <- list()

  for (i in 1:length(CoclG))
  {
    rowclust[i] <- HCclsG[CoclG[i]]
  }
  for (j in 1:length(CoclC))
  {
    colclust[j] <- HCclsC[CoclC[j]]
  }
  # -----------------------------------------------------------------------------------

  # Loop to get ranked co-cluster mean and their associated row and column cluster
  NGcls <- NULL
  NCcls <- NULL
  GC_CoMean <- NULL
  for (a in 1 : length(CoclG))
  {
    for (b in 1 : length(CoclC))
    {
      NGcls <- c(NGcls, a)
      NCcls <- c(NCcls, b)
      GCM <- mean(dataExpr[rowclust[[a]], colclust[[b]]])
      GC_CoMean <- c(GC_CoMean, GCM)
    }
  }
  # -----------------------------------------------------------------------------------

  GC_CoMeanR <- GC_CoMean[order(-GC_CoMean)]
  NG_Cocls <- NGcls[order(-GC_CoMean)]
  #Row.names <- paste0("R", NG_Cocls) # R indicates row (gene)

  NC_Cocls <- NCcls[order(-GC_CoMean)]
  #Col.names <- paste0("C", NC_Cocls) # C indicates column (compound)
  #len=1:length(Col.names)
  #Combine.Row.Col <- sapply(len,function(len) c(Row.names[len],Col.names[len]))

  NGC_Cocls <- paste(NG_Cocls, NC_Cocls, sep = ",")
  Coclust_MeanMat <- data.frame(NGC_Cocls, GC_CoMeanR)

  # Loop to get colors for the row and column clusters
  x <- unlist(rowclust)
  y <- unlist(colclust)
  colorsG <- rep(0, length(x))
  colorsC <- rep(0, length(y))
  for (k in 1:length(rowclust))
  {
    ind <- which(x %in% rowclust[[k]])
    if((k %% 2) == 0)
    {
      colorsG[ind] <- "green"
    } else {
      colorsG[ind] <- "red"
    }

  }
  for (L in 1:length(colclust))
  {
    inds <- which(y %in% colclust[[L]])
    if((L %% 2) == 0)
    {
      colorsC[inds] <- "green"
    } else {
      colorsC[inds] <- "red"
    }

  }
  # -----------------------------------------------------------------------------------

  # Re-arranged data matrix based on recovered clusters
  CoClsDtMat <- dataExpr[x, y]
  # -----------------------------------------------------------------------------------

  # For calculating moving range of co-cluster means to identify biomarker co-clusters
  mr <- NULL
  for (q in 1 : (length(GC_CoMeanR) - 1))
  {
    mr1 <- abs(GC_CoMeanR[q]-GC_CoMeanR[ q + 1])
    mr <- c(mr, mr1)
  }
  # -----------------------------------------------------------------------------------

  # Limits caculation for individual control charts
  CentralLine <- mean(GC_CoMeanR)
  Std.Dv <- mean(mr) / 1.128
  UpContLimit <- mean(GC_CoMeanR) + 3*Std.Dv
  LowrContLimit <- mean(GC_CoMeanR) - 3*Std.Dv
  # -----------------------------------------------------------------------------------

  # Pchmarks and colors for the individual control charts
  ind1 <- which(GC_CoMeanR >= UpContLimit)
  ind2 <- which((GC_CoMeanR >= LowrContLimit) & (GC_CoMeanR <= UpContLimit))
  ind3 <- which(GC_CoMeanR <= LowrContLimit)

  color <- rep(0, length(GC_CoMeanR))
  color[ind1] <- "red"
  color[ind2] <- "black"
  color[ind3] <- "blue"

  pchmark <- rep(0, length(GC_CoMeanR))
  pchmark[ind1] <- 19
  pchmark[ind2] <- 19
  pchmark[ind3] <- 19
  # -----------------------------------------------------------------------------------


  # return results
  return(list(Coclust_MeanMat = Coclust_MeanMat,
              CoClsDtMat = CoClsDtMat,
              NG_Cocls = NG_Cocls,
              NC_Cocls = NC_Cocls,
              rowclust = rowclust,
              colclust = colclust,
              colorsC = colorsC,
              colorsG = colorsG,
              CentralLine = CentralLine,
              UpContLimit = UpContLimit,
              LowrContLimit = LowrContLimit,
              color = color,
              pchmark = pchmark))
}
