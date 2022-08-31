# getkernel function
#
# @title Calculate the kernel matrix based on the Aitchison distance
#
# @param dv A numeric vector object. The distance vector for a specific OUT.
#   {nOTU}
#
# @param OTU A matrix object. The OTU table {n x nOTU}.
#
# @param h A numeric value. The h value of the correlation, c * sd {1}
#
# @keywords internal
.aitchison <- function(dv, h, OTU) {

  # get the correlation vector
  Rmc <- exp(x = -dv^2/(max(h,1e-10)))

  Rmc  <- as.vector(x = Rmc)
  
  # get the Aitchison distance matrix
  # {nsub x nsub}
  OTU <- sweep(x = OTU, MARGIN = 2L, STATS = sqrt(x = Rmc), FUN = "*")
  AD <- as.matrix(x = stats::dist(x = OTU))
  
  # get the kernel matrix
  # {nsub x nsub}
  nsub <- nrow(x = OTU)
  mat <- matrix(data = -1.0/nsub, nrow = nsub, ncol = nsub)
  diag(x = mat) <- 1.0 - {1.0/nsub}
  
  # {nsub x nsub}
  kernel <- -0.5 * {mat %*% AD^2 %*% mat}
  
  return( kernel )
}


# @title Estimate the kernel for each OTU
#
# @param cValue A numeric object. A single 'c' value.
#
# @param dInfo A matrix object. The OTU distance matrix. {nOTU x nOTU}
#
# @param OTU A matrix object. The OTU table. {n x nOTU}
#
.getKernel <- function(cValue, dInfo, OTU) {

  h <- cValue * dInfo$sd_dMatrix

  # for each OTU, an nsub x nsub matrix is returned
  # kernelObj is {nsub*nsub x nOTU}
  kernelObj <- apply(X = dInfo$dMatrix, MARGIN = 2L,
                     FUN = function(x, h, OTU) {
                       .aitchison(dv = x, h = h, OTU = OTU)
                     },
                     h = h,
                     OTU = OTU)

  return( kernelObj )
}






