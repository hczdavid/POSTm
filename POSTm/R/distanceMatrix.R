# @title Calculate all pairwise distance based on branch lengths
#
# @param tree A phylo object or a matrix. The phylogenetic tree or a
#   matrix of distance as provided by the user
#
# @param OTU A matrix object. The OTU table {n x nOTU}
#
#' @importFrom ape cophenetic.phylo
#' @importFrom stats sd
#
distanceMatrix <- function(tree, otunames) {

  # {nTips x nTips}
  if (is.null(x = tree)) {
    Dmat <- diag(x = length(x = otunames))
    rownames(x = Dmat) <- otunames
    colnames(x = Dmat) <- otunames
  } else if (is(object = tree, class2 = "phylo")) {
    Dmat <- as.matrix(x = ape::cophenetic.phylo(x = tree))
  } else {
    Dmat <- tree
  }
  
  # the choice of sd for the full distance matrix is intentional
  # it is believed that SD is fixed no matter if some OTUs are filtered
  sd_dMatrix <- sd(x = Dmat)

  # for the analysis, remove any distance provided for OTUs no in OTU table

  # {nOTU x nOTU}
  Dmat <- Dmat[otunames, otunames, drop = FALSE]

  return( list("dMatrix"    = Dmat,
               "sd_dMatrix" = sd_dMatrix) )
}

