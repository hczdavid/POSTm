# adapted from ACAT.R by Yaowu Liu available from
# https://github.com/yaowuliu/ACAT/blob/master/R/ACAT.R on 3/25/2021

#
# Aggregated Cauchy Association Test
#
# A p-value combination method using the Cauchy distribution.
#
# @param Pvals a numeric vector/matrix of p-values. When it is a matrix, each 
#     column of p-values is combined by ACAT.
#
# @return The p-value(s) of ACAT.
#
# @author Yaowu Liu
#
# @references Liu, Y., & Xie, J. (2019). 
#     Cauchy combination test: a powerful test with analytic p-value calculation
#     under arbitrary dependency structures. 
#     \emph{Journal of American Statistical Association},115(529), 393-402. 
#
# modified to suit the package needs on 3/25/2021
# primarily style preferences and the removal of ability to incorporate
# weights and ability to initiate checks
#
ACAT <- function(pvalues) {

  #### check if there is NA
  if (any(is.na(x = pvalues))) {
    stop("cannot have NAs in the p-values", call. = FALSE)
  }

  #### check if p-values are between 0 and 1
  if (any(pvalues < 0.0) && any(pvalues > 1.0)) {
    stop("p-values must be between 0 and 1", call. = FALSE)
  }

  #### check if there are pvalues that are either exactly 0 or 1.
  isZero <- any(pvalues < 1e-15)
  isOne <- any(pvalues > {1-1e-15})
  if (isZero && isOne) {
    stop("cannot have both 0 and 1 p-values", call. = FALSE)
  }

  if (isZero) warning("some p-values are exactly 0", call. = FALSE)
  if (isOne) warning("some p-values are exactly 1", call. = FALSE)

  #### check if there are very small non-zero p values and calcuate the 
  #### cauchy statistics
  isSmall <- pvalues < 1e-15
  pvalues[!isSmall] <- tan(x = {0.5 - pvalues[!isSmall]}*pi)
  pvalues[isSmall] <- 1.0 / pvalues[isSmall] / pi
  cctStat <- mean(x = pvalues)

  #### return the ACAT p value(s).
  pval <- stats::pcauchy(q = cctStat, lower.tail = FALSE)
  return( pval )
}
