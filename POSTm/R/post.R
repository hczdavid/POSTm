#' Phylogeney-Guided OTU-Specific Association Test for Microbiome Data
#' 
#' The POST implements a
#'   phylogeny-guided OTU-specific association test for microbiome data. This
#'   method boosts the 
#'   testing power by adaptively borrowing information from phylogenetically 
#'   close OTUs of the target OTU. Whether or not borrowing information or the 
#'   amount of information from the neighboring OTUs is data adaptive and 
#'   supervised by phylogenetic distance and the outcome variable. POST is 
#'   built on a kernel machine regression framework and inherited the 
#'   advantages including flexibly model complex microbiome effects (e.g., 
#'   effects from opposite direction), easily adjust for covariates, and 
#'   accommodate both continuous and binary outcomes. 
#'
#' It is assumed that the OTU table is defined by a 97\% sequence similarity. 
#'   Though this
#'   threshold is not (cannot be) enforced in the implementation, the
#'   recommended maximum c-value may not be appropriate for other thresholds.
#'
#' There are numerous packages available for generating phylogenetic trees.
#'   The sole purpose of the tree input is to enable the calculation of 
#'   pairwise distances between the pairs of tips using branch lengths.
#'   Because it is not feasible to support this functionality for all 
#'   packages that generate such trees, the package allows for the specification
#'   of the distance matrix as an alternative to providing a specifically
#'   formatted tree object. For objects of class "phylo", "hclust", and "phylog", 
#'   the ape package provides tools to obtain the distance matrix, and these
#'   tools are used in this implementation. For all others, the distance
#'   matrix (defined by branch lengths) must be provided by the user through 
#'   input tree.
#'
#' @param y A numerical vector. The outcome of interest. 
#'   Data can be binary or continuous.
#'
#' @param X A data.frame object, matrix object or NULL. The covariates data.
#'   If NULL, an intercept only model is assumed. Factor covariates are allowed.
#'
#' @param OTU A matrix object. The operational taxonomic units (OTU).
#'   Data can be provided as counts or as proportions.
#'   Each row indicates a single sample; each column a single OTU.
#'   NA/0 values are allowed, but their presence will trigger a shift of all
#'   data by a small internally defined value. 
#'   The matrix must include column headers providing unique
#'   identification for each OTU; these identifiers are expected to be included
#'   in the tip labels of the input tree object. Any identifiers that are not
#'   included as tip labels are removed from the analysis.
#'
#' @param tree An object of class "phylo", "hclust", "phylog", a matrix object
#'   or NULL. If NULL, only the single OTU test will be estimated.
#'   Objects of class "phylo", "hclust", and "phylog" are phylogenetic trees,
#'   the tip labels of which must include all of the
#'   identifiers used as column headers of OTU. If a matrix, a
#'   square symmetric matrix containing the pairwise distances between OTUs
#'   as defined by the branch lengths. Note that the full tree should be
#'   provided/used and should not be subset or truncated, even if OTU does
#'   not contain all tips.
#'   See details for further information.
#'
#' @param cValues A numeric vector. The c values at which p-values are to be
#'   estimated. The default is a vector of evenly spaced values between
#'   zero and the recommended maximum value for OTUs defined at
#'   97\% sequence similiary, 
#'   c_max = 0.05.
#'   If no tree is provided, cValues will be set to 0.
#'
#' @returns Returns a POST object. Analysis results are provided as
#'   a matrix, the contents of which will depend
#'   on the inputs selected for the analysis. Possible columns include
#'   \item{POST_pvalue}{A numeric object. The POST p-value.}
#'   \item{SO_pvalue}{A numeric object. The single OTU test p-value.}
#'   \item{BEST_C}{A numeric object. The c value corresponding to the minimum
#'     POST p-value.}
#'   Row names indicate the OTU to which the results pertain.
#'   Results are ordered according to the POST (if tree provided) or SO (if tree
#'   not provided) p-values.
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
#' @export
#' @include getKernel.R CKAT.R glmFit.R lmFit.R distanceMatrix.R
#' @importFrom ape as.phylo
#' @importFrom methods is
#' @references Huang, C., Callahan, B., Wu, M. C., Holloway, S. T., 
#'   Brochu, H., Lu, W., Peng, X., and Tzeng, J-Y. (2021). 
#'   Phylogeny-guided microbiome OTU-specific association test (POST). 
#'   Bioinformatics, under revision.

post <- function(y,
                 OTU,
                 tree = NULL,
                 X = NULL,
                 cValues = seq(from = 0, to = 0.05, by = 0.01)) {

  if (!is.matrix(x = OTU)) OTU  <- as.matrix(x = OTU)

  nsubj  <- nrow(x = OTU)
    
  # check if y is numeric
  if (!is.numeric(x = y)) {
    stop("response must be numeric")
  }
  
  # check the length of y, x and OTU
  if (length(y) != nsubj) {
    stop("number of samples in y does not match number of samples in OTU table")
  }
  
  # check if X is provided, if not specify intercept only model. If it is
  # provided ensure that it is a matrix object
  if (!is.null(x = X)) {
    if (is.data.frame(x = X)) {
      X <- data.matrix(frame = X)
    } else if (!is.matrix(x = X)) {
      X <- as.matrix(x = X)
    }
    if (length(x= y) != nrow(x = X)) {
      stop("number of samples in y does not match number of samples in x")
    }
    message(ncol(x = X), " covariate(s) included in outcome regression model")
  } else {
    X <- matrix(data = 1.0, nrow = nsubj, ncol = 1L)
    message("intercept only outcome regression model")
  }

  # check missing value of OTU
  if (any(is.na(x = OTU)) ) {
    OTU[is.na(x = OTU)] <- 0
    message("missing value(s) detected in OTU; NA replaced by 0")
  }
  
  # can replace the otu.type input with the following check
  # if integers are provided, then the rounded version will be equal
  # to the input version. If proportions are given, this will not be
  # true. 
  tst <- as.integer(x = round(x = OTU[,1L], digits = 0L))

  if (isTRUE(x = all.equal(target = OTU[,1L],  
                           current = tst, 
                           check.attributes = FALSE))) {
    otu.type <- "counts"
    if (any(OTU < 0.5)) OTU <- OTU + 0.5

  } else {
    otu.type <- "proportions"
    if (any(OTU < 1e-6)) OTU <- OTU + 1e-6
  }
  message("OTU type identified as ", otu.type)
  
  # can replace the y.type input with the following check
  # if integers are provided, then the rounded version will be equal
  # to the input version and there will only be two values. 
  # If continuous are given, this will not be
  # true. 
  tst <- as.integer(x = round(x = y, digits = 0L))
  if (isTRUE(x = all.equal(target = y, current = tst))) {
    if (length(x = unique(x = tst)) > 2L) {
      y.type <- "non-binary"
    } else if (length(x = unique(x = tst)) == 2L) {
      y.type <- "binary"
    } else {
      stop("could not identify response type")
    }
  } else {
    y.type <- "non-binary"
  }
  message("response type identified as ", y.type)
  
  # check if OTU name is given
  otunames <- colnames(x = OTU)
  if (is.null(x = otunames)) {
    stop("OTU column names are required and must match those provided in tree")
  }

  if (is.null(x = tree)) {

    cValues <- 0.0
    message("no tree provided; performing single OTU test")
    tiplabel <- otunames

  } else if (is(object = tree, class2 = "phylo") ||
             is(object = tree, class2 = "hclust") ||
             is(object = tree, class2 = "phylog")) {
    # input tree is a phylogenetic tree of an expected class

    # get tip labels
    tree <- ape::as.phylo(x = tree)
    tiplabel <- tree$tip.label

    message("tree input identified as a phylogenetic tree")

  } else if (is.matrix(x = tree)) {
    # tree is a distance matrix

    # must be square
    if (nrow(x = tree) != ncol(x = tree)) {
      stop("if provided as distance matrix, tree must be of dimension square")
    }

    # diagonal elements must be zero
    if (!all(diag(x = tree) < 1e-8)) {
      stop("if provided as distance matrix, diagonal elements must be zero")
    }

    # must be symmetric
    if (!isTRUE(all.equal(tree[lower.tri(tree)], t(tree)[lower.tri(tree)]))) {
      stop("if provided as distance matrix, tree must be symmetric")
    }

    # column headers must be provided
    tiplabel <- colnames(x = tree)
    if (is.null(x = tiplabel)) {
      stop("if provided as distance matrix, tree must have column names")
    }

    message("tree input identified as a distance matrix")
  } else {
    # tree is not recognized
    stop("tree must be of class 'phylo', 'hclust', 'phylog' or a matrix")
  }

  nmsIn <- otunames %in% tiplabel
        
  if (any(!nmsIn)) {

    # if there are any otunames not in tip labels remove from OTU matrix

    OTU <- OTU[,nmsIn, drop = FALSE]

    if (ncol(x = OTU) == 0L) {
      # if there is no overlap at all stop execution
      stop("column headers of OTU do not overlap with tip labels of tree")
    }

    message(sum(!nmsIn), " OTU(s) not found in tree tip labels; OTU(s) removed")
  }

  if (!is.numeric(x = cValues)) {
    stop("cValues must be numeric")
  }
  if (max(cValues) > 0.05) {
    warning("the maximum cvalue is larger than the recommended value of 0.05;",
            " this may result in more false positives.",
            " See Huang et al. (2021) discussion section for details.")
  }

  cValues <- sort(x = unique(x = c(0.0, cValues)))
        
  nOTU <- ncol(x = OTU)

  # CLR transformation
  OTU <- log(x = OTU) - 1.0/nOTU*rowSums(x = log(x = OTU))

  # get the distance matrix and the standard deviation
  dInfo <- distanceMatrix(tree = tree, otunames = otunames)

  # fit outcome of interest based on outcome type
  if (y.type == "binary") {

    fitObj = .glmFit(y = y, X = X)

  } else {

    fitObj = .lmFit(y = y, X = X)

  }

  postp <- numeric(length = nOTU)
  minc  <- numeric(length = nOTU)
  c0p   <- numeric(length = nOTU)
      
  pvCKatMat <- matrix(data = 0.0, nrow = length(x = cValues), ncol = nOTU)

  for (cv in 1L:length(x = cValues)) {

    message("cvalue: ", cValues[cv])

    # {nsubj*nsubj x nOTU}
    kernellist <- .getKernel(cValue = cValues[cv], dInfo = dInfo, OTU = OTU)

    message("OTU: ", appendLF = FALSE)

    for (m in 1L:nOTU) {

      
      if (m %% 10 == 0L) message(m, " ", appendLF = FALSE)
      if (m %% 100 == 0L) message()
      
      kern <- matrix(data = kernellist[,m], nrow = nsubj, ncol = nsubj)

      if (y.type == "binary") {
        
        pvCKAT <- CKAT.bin(K = kern, fitObj = fitObj)
        if (pvCKAT > 1.0) pvCKAT <- 1.0
        
      } else {
        
        pvCKAT <- CKAT.cont(K = kern, fitObj = fitObj)
        if (pvCKAT > 1.0) pvCKAT <- 1.0
        
      }

      pvCKatMat[cv,m] <- pvCKAT
    }
    message()

  }
      
  for (m in 1L:nOTU) {
    # if vector pvCKATv have both 0 and 1
    if (any(pvCKatMat[,m] >= {1.0-1e-15}) && any(pvCKatMat[,m] <= 1e-15)) {
      postp[m] <- 0.0
    } else {
      postp[m] <- ACAT(pvalues = pvCKatMat[,m])
    }
    minc[m] <- cValues[which.min(x = pvCKatMat[,m])]
    c0p[m] <- pvCKatMat[1L,m]
  }

  if (!is.null(x = tree)) {
    results  <- cbind(postp, c0p, minc)

    rownames(x = results) <- otunames
    colnames(x = results) <- c("POST_pvalue", "SO_pvalue", "Best_C")
  
    results  <- results[order(results[,1L]),]
    results <- round(results, 5L)

    class(results) <- "POST"

    attr(x = results, which = "tree") <- tree

  } else {
    results  <- matrix(data = c0p, ncol = 1L)

    rownames(x = results) <- otunames
    colnames(x = results) <- c("SO_pvalue")

    results  <- results[order(results[,1L]),,drop = FALSE]
    results <- round(results, 5L)

    class(results) <- "POST"

    attr(x = results, which = "tree") <- NA

  }

  return( results )
}

