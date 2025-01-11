#' @title Class "ARICluster-class" for storing results of adaptive thresholding algorithm
#' @name ARICluster-class
#' @docType class
#' @aliases ARICluster-class
#' @description The \code{ARICluster-class} class is the output of a call to \code{\link{ARICluster}}. It stores the information needed for the next answering-query step, where all maximal supra-threshold clusters can be quickly found given a TDP threshold.
#' @slot m Object of class "integer". Stores the total number of nodes.
#' @slot alpha Object of class "numeric". Stores the significance level.
#' @slot p Object of class "numeric". Stores input p-values in the original order.
#' @slot stcs Object of class "integer". Stores the representative nodes for all admissible supra-threshold clusters in the ascending order of TDP bounds.
#' @slot sizes Object of class "integer". Stores the subtree sizes (cluster extents) for all nodes in the ascending order of p-values.
#' @slot tdps Object of class "numeric". Stores the TDP lower confidence bounds for all supra-threshold clusters in the ascending order of p-values.
#' @slot childs Object of class "list". Stores a list of vectors of child nodes for all nodes in the ascending order of p-values.
#' @slot marks Object of class "integer". Initializes a vector to mark nodes within found clusters.
#' @return Returns an \code{\link{ARICluster-class}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel
#' @export
#' 
setClass("ARICluster",
         slots = list(
           m = "integer",       # stores number of p-values
           alpha = "numeric",   # stores significance level
           p = "numeric",       #(m) stores input unsorted p-values
           stcs = "integer",    #(m*<=m) stores representative nodes for all admissible STCs
           sizes = "integer",   #(m) stores subtree sizes for all nodes
           tdps = "numeric",    #(m) stores TDP lower bounds for all STCs, each represented by a node
           childs = "list",     #(m) stores children for all nodes
           marks = "integer"    #(m) stores marks for nodes within found clusters
         )
)


#' @title All-resolutions inference for cluster thresholding (ARICluster)
#' @name ARICluster
#' @aliases ARICluster
#' @description Based on ARI, the adaptive thresholding algorithm is employed to answer queries.
#' @usage ARICluster(p, adj, alpha = 0.05)
#' @param p A vector of p-values.
#' @param adj A list of neighbours for all nodes.
#' @param alpha Significance level. \code{alpha = 0.05} by default.
#' @return Returns an \code{\link{ARICluster-class}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel
#' @export
#' 
ARICluster <- function(p, adj, alpha=0.05) {
  
  # check for p
  if (missing(p)) stop("The p-value vector is not defined.")
  if (length(p)==0) stop("The p-value vector is empty.")
  if (!is.numeric(p)) stop("The p-value vector should be numeric.")
  if (anyNA(p)) stop("Found missing values in p-value vector. Please remove them!")
  if (min(p)<0 || max(p)>1) stop("P-values must be within [0,1].")
  # compute size of the multiple testing problem
  m <- length(p)
  # check for adj
  if (missing(adj)) stop("The adjacency list is not defined.")
  if (length(adj)==0) stop("The adjacency list is empty.")
  if (!is.list(adj) || !all(sapply(adj, is.integer))) stop("The adjacency list should be a list of integer vectors.")
  if (!all(Reduce(range, adj, c(1,m))==c(1,m))) stop("The elements of adjacency list must be within [1,m].")
  # check for p & adj
  if (length(p)!=length(adj)) stop("The length of p & adj does not match!")
  
  # perform hommel to find whole-set TDP bound
  hom <- hommel::hommel(p, simes=TRUE)
  tdp <- hommel::tdp(hom, alpha=alpha)
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
  rankp[ordp] <- seq_len(m)
  rankp <- as.integer(rankp)  # sorting ranks (starts from 1)
  
  # find all admissible STCs & compute TDP bounds
  reslist <- findClusters(m, adj, ordp, rankp)
  tdps    <- forestTDP(m, halpha, alpha, simeshalpha, p, reslist$SIZE, reslist$ROOT, reslist$CHILD)
  stcs    <- queryPreparation(m, reslist$ROOT, tdps, reslist$CHILD)
  
  # initialize a vector for marking found clusters
  marks <- integer(m)
  
  out <- new("ARICluster",
             m = m,
             alpha = alpha,
             p = p,
             stcs = stcs,
             sizes = reslist$SIZE,
             tdps = tdps,
             childs = reslist$CHILD,
             marks = marks)
  
  return(out)
}


#' @title Answer query with a TDP threshold
#' @name TDPQuery
#' @aliases TDPQuery
#' @description \code{TDPQuery} is a generic function used to quickly find all maximal supra-threshold clusters given a TDP threshold, where each cluster has the TDP not smaller than the threshold. The function invokes particular methods that depend on the class of the first argument.
#' @usage TDPQuery(aricluster, threshold)
#' @param aricluster An \code{\link{ARICluster-class}} or \code{\link{ARIBrainCluster-class}} object.
#' @param threshold A TDP threshold for forming maximal clusters.
#' @examples
#' 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#' 
#' # (1) create an ARIBrainCluster-class object
#' aricluster <- ARIBrainCluster(Pmap = pvalue_name, mask = mask_name)
#'
#' # (2) answer query: find all maximal clusters given a TDP threshold
#' tdpclusters <- TDPQuery(aricluster, threshold = 0.7)
#' summary(tdpclusters)
#'     
#' @return Returns a \code{\link{TDPClusters}}or \code{\link{TDPBrainClusters}} object of cluster list.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import 
#' @export
#'
setGeneric("TDPQuery", function(aricluster, threshold) standardGeneric("TDPQuery"))
setMethod("TDPQuery", "ARICluster", function(aricluster, threshold) {
  
  # check for threshold
  if (missing(threshold)) stop("A TDP threshold should be provided.")
  
  # find all maximal STCs
  clusterlist <- answerQuery(threshold, aricluster@stcs,
                             aricluster@sizes, aricluster@marks,
                             aricluster@tdps, aricluster@childs)
  n <- length(clusterlist)
  if (n==0) stop("No clusters were found for TDP threshold = ", threshold)
  # Sort clusters by descending cluster size.
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
  
  out <- new("TDPClusters",
             aricluster = aricluster,
             threshold = threshold,
             clusterlist = clusterlist)
  
  return(out)
})


#' @title Find all local minima
#' @name findLMs
#' @aliases findLMs
#' @description \code{findLMs} is a generic function used to search for all local minima of p-values.
#' @usage findLMs(aricluster)
#' @param aricluster An \code{\link{ARICluster-class}} or \code{\link{ARIBrainCluster-class}} object.
#' @return Returns a vector of all local minima.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import 
#' @export
#'
setGeneric("findLMs", function(aricluster) standardGeneric("findLMs"))
setMethod("findLMs", "ARICluster", function(aricluster) {
  # find node indices of all local minima
  out <- findLMS(aricluster@childs)+1
  return(out)
})
