#' Print the Primary Results of a post() Analysis
#'
#' Print the primary results of a post() analysis.
#'
#' @param x A POST object. The value object returned by a call to post()
#'
#' @param ... Further arguments passed to or from other methods
#'
#' @param siglevel A numeric object. If < 1, shows only the OTUs for which
#'   the p-value is below the specified siglevel.
#'
#' @export
#' @name print
#' @method print POST
#'
#' @returns No return value, called to display key results.
#'
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
#' print(x = result)
#' 
print.POST <- function(x, ..., siglevel = 1.0) {

  attr(x = x, which = "tree") <- NULL
  x <- unclass(x = x)

  if (siglevel < 1.0) {

    # if POST analysis included
    if (ncol(x = x) == 3L) {

      showi <- x[,1L] <= siglevel
      if (any(showi)) {
        cat("POST: p-values <= ", siglevel, "\n")
        print(x = x[showi,c(1L:3L),drop=FALSE])
      } else {
        cat("POST: No p-values <= ", siglevel, "\n")
      }

      showi <- x[,2L] <= siglevel
      if (any(showi)) {
        cat("\nSingle OTU test: p-values <= ", siglevel, "\n")
        print(x = x[showi,2L,drop=FALSE])
      } else {
        cat("\nSingle OTU test: No p-values <= ", siglevel, "\n")
      }

    } else {
      # only single OTU test included in original analysis
      showi <- x[,1L] <= siglevel
      if (any(showi)) {
        cat("Single OTU test: p-values <= ", siglevel, "\n")
        print(x = x[showi,1L,drop=FALSE])
      } else {
        cat("Single OTU test: No p-values <= ", siglevel, "\n")
      }
    }
  } else {
    print(x = x)
  }
}
