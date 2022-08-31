# Code adapted from CKAT package as described in
# Zhan, X., Girirajan, S., Zhao, N., Wu, M. C., and Ghosh, D. (2016)
# A novel copy number variants kernel association test with application
# to autism spectrum disorders studies. Bioinformatics, 32, 3603--3610.
#
#' @keywords internal
#' @importFrom CompQuadForm davies liu
CKAT.bin <- function (K,  fitObj) {

  DKD <- tcrossprod(x = fitObj$sqrtw) * K
  
  tQK <- t(x = fitObj$Q) %*% DKD
  
  QtQK <- fitObj$Q %*% tQK

  PKP1 <- DKD - QtQK - t(x = QtQK) + 
          fitObj$Q %*% (tQK %*% fitObj$Q) %*% t(x = fitObj$Q)

  q1 <- as.numeric(x = fitObj$res %*% K %*% fitObj$res)

  q1 <- q1 / fitObj$adj

  ee1 <- eigen(PKP1 - q1 * fitObj$P0, symmetric = TRUE, only.values = TRUE)

  lambda1 <- ee1$values[abs(ee1$values) >= 1e-10]

  p1 <- suppressWarnings(tryCatch(expr = davies(q = 0, lambda = lambda1, acc = 1e-6)$Qq,
                                  condition = function(ew) { return(NULL) }))

  if (is.null(x = p1)) {
    p1 <- suppressWarnings(tryCatch(expr = liu(q = 0, lambda = lambda1),
                                    error = function(e) { return( NULL ) }))
    if (is.null(x = p1)) {
      p1 <- davies(q = 0, lambda = lambda1, acc = 1e-6)$Qq
    }
  }

  return( p1 )
}

#' @keywords internal
#' @importFrom CompQuadForm davies liu
CKAT.cont <- function (K, fitObj) {

  PKP <- fitObj$P0 %*% K %*% fitObj$P0

  q <- as.numeric(x = fitObj$res %*% K %*% fitObj$res / fitObj$s2)

  ee <- eigen(PKP - q * fitObj$P0, symmetric = TRUE)

  lambda0 <- ee$values[abs(ee$values) >= 1e-10]

  p1 <- suppressWarnings(tryCatch(expr = davies(q = 0, lambda = lambda0, acc = 1e-6)$Qq,
                                  condition = function(ew) { return(NULL) }))

  if (is.null(x = p1)) {
    p1 <- suppressWarnings(tryCatch(expr = liu(q = 0, lambda = lambda0),
                                    error = function(e) { return( NULL ) }))
    if (is.null(x = p1)) {
      p1 <- davies(q = 0, lambda = lambda0, acc = 1e-6)$Qq
    }
  }

  return( p1 )
}

