# This is the function for generating simulated data
#' @export
simulate_data <- function (no.gene, no.dcc)
{
  x1 <- rbind(matrix(rnorm(10*10, 3, 0.35), 10, 10), matrix(rnorm(10*10, -3, 0.35), 10, 10), matrix(rnorm(30*10, 0, 0.35), 30, 10))
  x2 <- rbind(matrix(rnorm(10*10, 0, 0.35), 10, 10), matrix(rnorm(10*10, 3, 0.35), 10, 10), matrix(rnorm(10*10, -3, 0.35), 10, 10),
              matrix(rnorm(20*10, 0, 0.35), 20, 10))
  x3 <- matrix(rnorm(50*16, 0, 0.35), 50, 16)

  SimData <- cbind(x1, x2, x3)

  GeneName <- paste("G", 1 : no.gene, sep="")
  DrugName <- paste("C", 1 : no.dcc, sep="")
  DrugNameL <- paste(DrugName, "Low", sep=("_"))
  DrugNameM <- paste(DrugName, "Middle", sep=("_"))
  DrugNameH <- paste(DrugName, "High", sep=("_"))
  DnameWdose <- c(DrugNameH[1:5], DrugNameM[1:5], DrugNameH[6:10], DrugNameM[6:10], DrugNameH[11:12], DrugNameM[11:12],
                  DrugNameL)
  dimnames(SimData) <- list(GeneName, DnameWdose)
  gname <- rownames(SimData)
  cname <- colnames(SimData)

  SimDataRnd <- SimData[sample(nrow(SimData)), sample(ncol(SimData))]
  GCmat <- 100*(1 / (1 + exp(-(SimDataRnd))))

  return(list(SimData = SimData,
              SimDataRnd = SimDataRnd,
              GCmat = GCmat))
}
