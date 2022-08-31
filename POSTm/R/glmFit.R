# extracted from CKAT.bin from CKAT package and adapted to suit package needs
#
# Debashis Ghosh and Xiang Zhan. "CKAT software" (2016)
# Available at: http://works.bepress.com/debashis_ghosh/75/
#
#' @importFrom stats model.matrix glm
.glmFit <- function (y, X) {
  
  n <- length(x = y)

  # obtain model matrix
  X1 <- stats::model.matrix(object = ~. , data = as.data.frame(x = X))
  
  # logistic regression
  glmfit <- stats::glm(y ~ X1-1, family = "binomial")

  mu  <- glmfit$fitted.values
  res.wk <- glmfit$residuals
  res <- y - mu

  w <- mu * {1.0 - mu}
  sqrtw <- sqrt(x = w)

  adj <- sum({sqrtw * res.wk}^2)

  DX12 <- sqrtw * X1

  qrX <- qr(x = DX12, tol = 1e-7)
  Q <- qr.Q(qr = qrX)
  Q <- Q[, 1L:qrX$rank, drop = FALSE]

  P0 <- diag(x = n) - tcrossprod(x = Q)
  
  return( list("Q" = Q, 
               "sqrtw" = sqrtw,  
               "P0" = P0,  
               "res" = res,  
               "adj" = adj) )

}
