#' @title All-resolutions inference (ARI) for cluster thresholding
#' @name ARICluster
#' @aliases ARICluster
#' @description Based on ARI, the adaptive thresholding algorithm is employed to answer queries.
#' @usage ARICluster(p, adj, alpha = 0.05)
#' @param p A vector of p-values.
#' @param adj A list of neighbours (unsorted orders that start from 1) for all nodes.
#' @param alpha Significance level. \code{alpha = 0.05} by default.
#' @return Returns an \code{\link{ARICluster-class}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel
#' @export

ARICluster <- function(p, adj, alpha=0.05) {
  
  # check for p
  if (min(p)<0 || max(p)>1) stop("P-values must be within [0,1].")
  # check for p & adj
  if (length(p)!=length(adj)) stop("The length of p & adj does not match!")
  
  # compute size of the multiple testing problem
  m <- length(p)

  # perform hommel to find whole-set TDP bound
  hom <- hommel::hommel(p, simes=TRUE)
  tdp <- hommel::tdp(hom)
  if (tdp==0) {
    stop("The full-set TDP bound is zero!") 
  } else {
    cat("The", 1-alpha, "confidence lower bound for the TDP concerning all nodes is", tdp, "\n")
  }
  # compute h(alpha)
  halpha      <- hommel:::findHalpha(hom@jumpalpha, alpha, m)
  simeshalpha <- hom@simesfactor[halpha+1]
  
  #
  # run adaptive thresholding algorithm
  #
  
  # sort input p-values in ascending order
  ordp  <- order(p)
  ordp  <- as.integer(ordp)   # sorted orders (starts from 1)
  # find the sorting ranks for unsorted p-values
  rankp <- integer(m)
  rankp[ordp] <- 1:m
  rankp <- as.integer(rankp)  # sorting ranks (starts from 1)
  
  # find all admissible STCs & compute TDP bounds
  reslist <- findClusters(m, adj, ordp, rankp)
  tdps    <- forestTDP(m, halpha, alpha, simeshalpha, p, ordp, reslist$SIZE, reslist$ROOT, reslist$CHILD)
  stcs    <- queryPreparation(m, reslist$ROOT, tdps, reslist$CHILD)
  
  # initialize a vector for marking found clusters
  marks <- integer(m)
  
  out <- new("ARICluster",
             m = m,
             alpha = alpha,
             p = p,
             ordp = ordp,
             stcs = stcs,
             sizes = reslist$SIZE,
             tdps = tdps,
             childs = reslist$CHILD,
             marks = marks)
  
  return(out)
}
