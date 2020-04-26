#' @importFrom stats cutree dist hclust rnorm
#' @importFrom fields image.plot
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline image mtext par plot text

#' @export
rhcoclust <- function(data, rk, ck)
{
  # Data Transformation using logistic function
  dataExpr <- 100*(1/(1+exp(-data)))
  dG <- dist((dataExpr), method = "manhattan")
  dC <- dist(t(dataExpr), method = "manhattan")

  HCGene <- hclust(dG, method = "ward.D")
  HCComp <- hclust(dC, method = "ward.D")
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
    HCclsG[[p]] <- rownames(HCclsMatG)[which(HCclsMatG[,1]==p)]
  }
  for (q in 1:ck)
  {
    HCclsC[[q]] <- rownames(HCclsMatC)[which(HCclsMatC[,1]==q)]
  }
  # -----------------------------------------------------------------------------------

  # Loop to get co-cluster means and corresponding row and column cluster number
  Gcls <- NULL
  Ccls <- NULL
  Co_ClsMeanExpr <- NULL
  for (g in 1:length(HCclsG))
  {
    for (c in 1:length(HCclsC))
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
  Co_Gcls <- list()
  Co_Ccls <- list()

  for (i in 1:length(CoclG))
  {
    Co_Gcls[i] <- HCclsG[CoclG[i]]
  }
  for (j in 1:length(CoclC))
  {
    Co_Ccls[j] <- HCclsC[CoclC[j]]
  }
  # -----------------------------------------------------------------------------------

  # Loop to get ranked co-cluster mean and their associated row and column cluster
  NGcls <- NULL
  NCcls <- NULL
  GC_CoMean <- NULL
  for (a in 1:length(CoclG))
  {
    for (b in 1:length(CoclC))
    {
      NGcls <- c(NGcls, a)
      NCcls <- c(NCcls, b)
      GCM <- mean(dataExpr[Co_Gcls[[a]],Co_Ccls[[b]]])
      GC_CoMean <- c(GC_CoMean, GCM)
    }
  }
  # -----------------------------------------------------------------------------------

  GC_CoMeanR <- GC_CoMean[order(-GC_CoMean)]
  NG_Cocls <- NGcls[order(-GC_CoMean)]
  NC_Cocls <- NCcls[order(-GC_CoMean)]
  NGC_Cocls <- paste(NG_Cocls, NC_Cocls, sep=",")
  NGC_cls_MeanMat <- data.frame(NGC_Cocls, GC_CoMeanR)

  # Loop to get colors for the row and column clusters
  x <- unlist(Co_Gcls)
  y <- unlist(Co_Ccls)
  colorsG <- rep(0,length(x))
  colorsC <- rep(0,length(y))
  for (k in 1:length(Co_Gcls))
  {
    ind <- which(x %in% Co_Gcls[[k]])
    if((k %% 2) == 0)
    {
      colorsG[ind] <- "green"
    } else {
      colorsG[ind] <- "red"
    }

  }
  for (L in 1:length(Co_Ccls))
  {
    inds <- which(y %in% Co_Ccls[[L]])
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
  for (q in 1:(length(GC_CoMeanR)-1))
  {
    mr1 <- abs(GC_CoMeanR[q]-GC_CoMeanR[q+1])
    mr <- c(mr,mr1)
  }
  # -----------------------------------------------------------------------------------

  # Limits caculation for individual control charts
  CentralLine <- mean(GC_CoMeanR)
  Std.Dv <- mean(mr)/1.128
  UpContLimit <- mean(GC_CoMeanR)+3*Std.Dv
  LowrContLimit <- mean(GC_CoMeanR)-3*Std.Dv
  # -----------------------------------------------------------------------------------

  # Pchmarks and colors for the individual control charts
  ind1 <- which(GC_CoMeanR >= UpContLimit)
  ind2 <- which((GC_CoMeanR >= LowrContLimit)&(GC_CoMeanR <= UpContLimit))
  ind3 <- which(GC_CoMeanR <= LowrContLimit)

  color <- rep(0,length(GC_CoMeanR))
  color[ind1] <- "red"
  color[ind2] <- "black"
  color[ind3] <- "red"

  pchmark <- rep(0,length(GC_CoMeanR))
  pchmark[ind1] <- 19
  pchmark[ind2] <- 19
  pchmark[ind3] <- 19
  # -----------------------------------------------------------------------------------


  # return results
  return(list(NGC_cls_MeanMat = NGC_cls_MeanMat,
              CoClsDtMat = CoClsDtMat,
              Co_Gcls = Co_Gcls,
              Co_Ccls = Co_Ccls,
              colorsC = colorsC,
              colorsG = colorsG,
              CentralLine = CentralLine,
              UpContLimit = UpContLimit,
              LowrContLimit = LowrContLimit,
              color = color,
              pchmark = pchmark ))
}
