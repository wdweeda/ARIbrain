#' @title Answer queries with a TDP threshold
#' @name TDPquery
#' @aliases TDPquery
#' @description Given a TDP threshold, it quickly finds all maximal supra-threshold clusters, each with the TDP not smaller than the threshold.
#' @usage TDPquery(ARIBrainCluster, gamma)
#' @param ARIBrainCluster An \code{\link{ARIBrainCluster-class}} object.
#' @param gamma A TDP threshold for forming maximal clusters.
#' @examples
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#' 
#' # (1) create an ARIBrainCluster object
#' ari <- ARIBrainCluster(Pmap = pvalue_name, mask = mask_name)
#'
#' # (2) answer queries: find all maximal clusters given a TDP threshold
#' res <- TDPquery(ari, gamma = 0.7)
#' res@clusterlist  # access cluster list (sorting ranks that start from 0)
#' 
#' # (3) summarize cluster information and write image file
#' summaryCluster(ari, res, rest = TRUE)
#' writeCluster(ari, res, file = "ari.nii.gz", template = mask_name)
#'     
#' @return Returns a \code{\link{TDPquery-class}} object of cluster list.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import 
#' @export
#' 
TDPquery <- function(ARIBrainCluster, gamma) {
  
  # check for gamma
  if (missing(gamma)) stop("TDP threshold gamma should be provided.")
  
  # find all maximal STCs
  clusterlist  <- answerQuery(gamma, ARIBrainCluster@stcs,
                              ARIBrainCluster@ordp, ARIBrainCluster@sizes,
                              ARIBrainCluster@marks, ARIBrainCluster@tdps,
                              ARIBrainCluster@childs)
  n <- length(clusterlist)
  if (n==0) stop("No clusters were found for gamma=", gamma)
  # sort clusters by descending cluster size.
  # R's built-in order function defaults to counting sort for integer vectors 
  # of length >= 200 and range < 100000, and uses radix sort when the range 
  # condition is not met. We therefore estimate whether radix sort will be 
  # output-sensitive, and otherwise ensure a counting sort is used, which is 
  # output-sensitive.
  if (n>1) {
    cluster_sizes <- sapply(clusterlist, length)
    d             <- diff(range(cluster_sizes))
    maxsize       <- max(cluster_sizes)
    if (n<200 || d<100000 || n*log2(d)<=sum(cluster_sizes)) {
      clusterlist <- clusterlist[order(cluster_sizes, decreasing=TRUE)]
    } else {
      cluster_ord <- counting_sort(n, maxsize, cluster_sizes)
      clusterlist <- clusterlist[cluster_ord+1]
    }
  }
  
  out <- new("TDPquery",
             gamma = gamma,
             clusterlist = clusterlist)
  
  return(out)
}

