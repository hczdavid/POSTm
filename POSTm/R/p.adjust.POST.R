#' Adjust P-values for Multiple Comparisons
#'
#' Given a set of p-values, returns p-values adjusted using one of
#'  several methods.
#'
#' @param p A numeric object. A vector of p-values.
#'
#' @param ... additional arguments required for additional methods
#'
#' @importFrom stats p.adjust
#' @name p.adjust
#' @export
p.adjust <- function (p, ...) {
   UseMethod("p.adjust", p)
 }

#' @export
p.adjust.default <- function(p, ...) {
  stats::p.adjust(p = p, ...)
}


#' Adjust P-values for Multiple Comparisons
#'
#' Given an object of class POST, returns the POST and SO p-values
#'   adjusted using the specified method(s).
#'
#' Extends the p.adjust() method of stats and incorporates the TSBH
#'  method available through multtest::mt.raw2adjp. 
#'
#' @param p An object of class POST. The value object returned by post().
#'
#' @param method A character vector object. One or more correction methods.
#'   must be one or more from {"holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'   "fdr", "none", "TSBH"}.
#'
#' @param n number of comparisons, must be at least ‘length(p)’; 
#'   Per ?stats::p.adjust -- only set this (to non-default) when you know 
#'   what you are doing!
#'
#' @param alpha A numeric object. A nominal type I error rate, or a vector of 
#'   error rates, used for estimating the number of true null hypotheses in 
#'   the two-stage Benjamini & Hochberg procedure ("TSBH"). Default is 0.05.
#'   NOTE: This definition does not conform to stats::p.adjust.
#' 
#' @param ... Ignored.
#' 
#' @returns A list, the contents of which will depend on the type of analysis.
#'   If a tree object was provided for the analysis, the list will contain
#'   \item{adjPOST}{A matrix object. The adjusted p-values for the POST p-values.
#'     The rows correspond to each OTU; the
#'     columns to each specified adjustment method.}
#'   \item{adjSO}{A matrix object. The adjusted p-values for the single OTU p-values.
#'     The rows correspond to each OTU; the
#'     columns to each specified adjustment method.}
#'  If no tree object was provided, only adjSO as described above will be returned.
#'
#'
#' @importFrom stats p.adjust.methods p.adjust
#' @rdname p.adjust
#' @method p.adjust POST
#' @export
#' @examples
#'
#' data("POSTmData")
#'
#' y <- as.integer(x = metadata[,"GC"] == "BV")
#' X <- metadata[,"mRace"]
#'
#' result <- post(y = y, 
#'                X = X, 
#'                OTU = otu[,1:20], 
#'                tree = otutree,
#'                cValues = seq(0,0.05,by=0.01))
#'
#'  p.adjust(p = result, method = c("BH","BY"))
#' 
p.adjust.POST <- function(p, ..., method = p.adjust.methods, 
                          n = length(x = p), alpha = 0.05) {

  if (requireNamespace("multtest", quietly = TRUE)) {
    allMethods <- c(p.adjust.methods, "TSBH")
  } else {
    allMethods <- p.adjust.methods
    if ("TSBH" %in% method) {
      warning("multtest is not available; TSBH cannot be used")
      method <- method[method %in% allMethods]
    }
  }

  if (isTRUE(all.equal(method, p.adjust.methods))) {
    method <- allMethods
  }

  if ('none' %in% method) {
    method <- method[!(method %in% 'none')]
  }

  if (ncol(x = p) == 3L) {

    adjPOST <- matrix(data = p[,1L], ncol = 1L)

    i <- 1L
    while (i <= length(x = method)) {

      if (method[i] == "TSBH") {
        adjPOST <- cbind(adjPOST, 
                         multtest::mt.rawp2adjp(rawp = p[,1L],
                                                proc = "TSBH", 
                                                alpha = alpha)$adjp[,2L])
      } else {

        adjPOST <- cbind(adjPOST, 
                         stats::p.adjust(p = p[,1L],
                                         method = method[i],
                                         n = n))

      }
      i <- i + 1L
    }
    colnames(x = adjPOST) <- c("raw", method)
    rownames(x = adjPOST) <- rownames(p)

    adjSO <- matrix(data = p[,2L], ncol = 1L)
    i <- 1L
    while (i <= length(x = method)) {

      if (method[i] == "TSBH") {
        adjSO <- cbind(adjSO, 
                       multtest::mt.rawp2adjp(rawp = p[,2L],
                                              proc = "TSBH", 
                                              alpha = alpha)$adjp[,2L])
      } else {
        adjSO <- cbind(adjSO, 
                       stats::p.adjust(p = p[,2L],
                                       method = method[i],
                                       n = n))
      }
      i <- i + 1L
    }

    colnames(x = adjSO) <- c("raw", method)
    rownames(x = adjSO) <- rownames(p)

    return( list("adjPOST" = adjPOST, "adjSO" = adjSO) )

  } else {

    adjSO <- matrix(data = p[,1L], ncol = 1L)

    i <- 1L
    while (i <= length(x = method)) {

      if (method[i] == "TSBH") {
        adjSO <- cbind(adjSO, 
                       multtest::mt.rawp2adjp(rawp = p[,1L],
                                              proc = "TSBH", 
                                              alpha = alpha)$adjp[,2L])
      } else {
        adjSO <- cbind(adjSO, 
                       stats::p.adjust(p = p[,1L],
                                       method = method[i],
                                       n = n))
      }
      i <- i + 1L
    }

    colnames(x = adjSO) <- c("raw", method)
    rownames(x = adjSO) <- rownames(p)

    return( list("adjSO" = adjSO) )
  }

}
