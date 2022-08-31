#' Plot Phylogenies with Significant OTU
#'
#' Plot phylogenetic tree on the current graphical device as defined by
#'  package ape with tip labels highlighted to show the significant OTU.
#'
#' @param x A POST object. The value object returned by post(). Note that
#'   a tree of class "phylo" must have been provided to the original
#'   analysis.
#'
#' @param ... Additional arguments to be passed to plot.phylo.
#' 
#' @param siglevel A numeric object. The significance level; OTU's with 
#'   (adjusted or raw) p-values below this level will be highlighted.
#'   Default is 0.05.
#'
#' @param method A character object. Must be one of \{"holm", "hochberg", 
#'   "hommel", "bonferroni", "BH", "BY", "fdr", "TSBH", "none"\}, where
#'   'none' indicates the raw POST p-values. "TSBH" is the
#'   two-stage Benjamini & Hochberg procedure. See ?stats::p.adjust for 
#'   descriptions of all other possible values.
#'
#' @param alpha A numeric object. A nominal type I error rate, or a vector of 
#'   error rates, used for estimating the number of true null hypotheses in 
#'   the two-stage Benjamini & Hochberg procedure ("TSBH"). Default is 0.05.
#'   Ignored if method is not TSBH.
#'
#' @param subTree A logical object. If TRUE, only the OTU used in the
#'   analysis are included in the plot.
#' 
#' @method plot POST
#'
#' @returns No return value, called to produce graphical elements.
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
#' plot(x = result)
#' 
#' @name plot
#' @export
#' @importFrom graphics plot
#' @importFrom ape plot.phylo tiplabels keep.tip
#'
#' @include p.adjust.POST.R
plot.POST <- function(x, ..., 
                      siglevel = 0.05, 
                      method = "none", 
                      alpha = 0.05,
                      subTree = TRUE) {

  tree <- attr(x = x, which = "tree")

  if (ncol(x = x) == 1L) stop("no tree provided in original analysis")

  if (!is(object = tree, class2 = "phylo")) {
    stop("tree is not of an appropriate class")
  }

  if (subTree) {
    tree <- ape::keep.tip(phy = tree, tip = rownames(x = x))
  }

  args <- list(...)
  args[[ "x" ]] <- tree
  args[[ "type" ]] <- "fan"
  args[[ "use.edge.length" ]] <- TRUE
  args[[ "show.tip.label" ]] <- FALSE

  if (tolower(x = method) == 'none') {
    args[[ "main" ]] <- "Raw p-values"
    show <- which(x = x[,1L] <= siglevel)
  } else {
    if ("TSBH" == method[1L]) {
      if (!requireNamespace("multtest", quietly = TRUE)) {
        stop("multtest is not available; TSBH cannot be used", call. = FALSE)
      }
    }
    pvalues <- p.adjust(p = x, method = method[1L], alpha = alpha)
    args[[ "main" ]] <- paste(method[1L], "adjusted p-values")
    show <- which(x = pvalues$adjPOST[,2L] <= siglevel)
  }

  do.call(plot, args = args)

  if (any(show)) {

    tiplabels(tree$tip.label[tree$tip.label %in% rownames(x = x)[show]],
              which(tree$tip.label %in% rownames(x = x)[show]),
              adj = c(0,0), cex = 0.5,
              col = "red", frame = "none")

  } else {
    warning("no p-values below provided significant level")
  }
}
