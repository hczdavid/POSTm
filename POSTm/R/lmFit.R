# extracted from CKAT.bin from CKAT package and adapted to suit package needs
#
# Debashis Ghosh and Xiang Zhan. "CKAT software" (2016)
# Available at: http://works.bepress.com/debashis_ghosh/75/
#
#' @importFrom stats model.matrix lm residuals
.lmFit <- function (y, X) {

  n <- length(x = y)
  
  X1 <- stats::model.matrix(~. , as.data.frame(X))

  mod <- stats::lm(y ~ X1 - 1)
  res <- residuals(object = mod)
  s2 <- sum(res^2)

  P0  <- diag(x = n)  - X1 %*% solve(a = crossprod(x = X1), b = t(x = X1))
    
  return( list("P0" = P0, "res" = res, "s2" = s2) )

}
