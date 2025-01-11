#' @title Class "TDPClusters" for outputting cluster information
#' @name TDPClusters
#' @docType class
#' @aliases TDPClusters
#' @description The class \code{TDPClusters} is the output of a call to \code{\link{TDPQuery}} or \code{\link{TDPChange}}. It stores the resulting cluster information.
#' @slot aricluster Object of class "ARICluster-class". Stores the \code{\link{ARICluster-class}} object, usually, a result of a call to \code{\link{ARICluster}}.
#' @slot threshold Object of class "numeric". Stores the TDP threshold.
#' @slot clusterlist Object of class "list". Stores a list of found clusters, each including indices of nodes within that cluster. Here, the node index starts from 0. Please use \code{tdpclusters@clusterlist[[i]]+1} to access node indices for the ith largest cluster.
#' @return Returns a \code{\link{TDPClusters}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import
#' @export
#' 
setClass("TDPClusters",
         slots = list(
           aricluster = "ARICluster",  # stores an ARICluster object
           threshold = "numeric",      # stores TDP threshold
           clusterlist = "list"        #(n) stores a list of found clusters
         )
)


#' @title Summarizing cluster information and generating summary table for found clusters.
#' @name summary.TDPClusters
#' @aliases summary.TDPClusters
#' @description \code{summary} method for class \code{\link{TDPClusters}}.
#' @usage summary(object, ..., rest = FALSE)
#' @param object Object of class "TDPClusters". Stores the \code{\link{TDPClusters}} object, usually, a result of a call to \code{\link{TDPQuery}} or \code{\link{TDPChange}}.
#' @param rest Object of class "logical". By default, \code{rest = FALSE} indicates that information on the rest nodes will not be shown in the summary table.
#' @details If the output is not assigned, the summary table will be printed on console.
#' @examples
#' 
#' table <- summary(tdpclusters)      # get summary table (without rest)
#' summary(tdpclusters, rest = TRUE)  # print summary table (with rest)
#' 
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @return Returns a summary table reporting Size, TDN (lower bound), #TrueNull (upper bound), TDP (lower bound), maximum statistic and the corresponding ID for each found cluster.
#' @import plyr
#' @export
#'
setMethod("summary", "TDPClusters", function(object, ..., rest=FALSE) {
  # compute number of clusters
  n <- length(object@clusterlist)
  # apply summaries to each cluster
  if (n>0) {
    sumtable <- plyr::laply(1:n, function(i) {
      clus_stat <- -qnorm(object@aricluster@p[object@clusterlist[[i]]+1])
      id_clus   <- which.max(clus_stat)
      
      clus_size <- length(object@clusterlist[[i]])
      clus_tdp  <- object@aricluster@tdps[object@clusterlist[[i]][clus_size]+1]
      unlist(c(Size=clus_size, 
               FalseNull=round(clus_size*clus_tdp), 
               TrueNull=round(clus_size*(1-clus_tdp)), 
               ActiveProp=clus_tdp,
               maxZ=clus_stat[id_clus],
               maxID=object@clusterlist[[i]][id_clus]+1))
    })
  } else {
    stop("No clusters were found for TDP threshold = ", object@threshold)
  }
  if (is.null(dim(sumtable))) sumtable <- t(as.matrix(sumtable))
  
  # modify output summary table by adding the "rest" information
  if (rest) {
    ord_rest  <- (1:object@aricluster@m)[-(unlist(object@clusterlist)+1)]
    rest_size <- length(ord_rest)
    if (rest_size>0) {
      rest_stat <- -qnorm(object@aricluster@p[ord_rest])
      id_rest   <- which.max(rest_stat)
      
      hom       <- hommel::hommel(object@aricluster@p, simes=TRUE)
      rest_disc <- hommel::discoveries(hom, ix=ord_rest, alpha=object@aricluster@alpha, incremental=FALSE)
      
      sumtable  <- rbind(sumtable, c(rest_size,
                                     rest_disc,
                                     rest_size-rest_disc,
                                     rest_disc/rest_size,
                                     rest_stat[id_rest],
                                     ord_rest[id_rest]))
    } else {
      sumtable <- rbind(sumtable, c(0,0,0,NA,NA,NA))
      warning("The rest of the volume does not contain any nodes.")
    }
    rownames(sumtable) <- c(paste0("Cluster", n:1), "rest")
  } else {
    rownames(sumtable) <- paste0("Cluster", n:1)
  }
  # update column names
  colnames(sumtable) <- c("Size", "TDN(lower)", "#TrueNull(upper)", "TDP(lower)", "Stat", "ID")
  
  sumtable
})


#' @title Displaying cluster summary
#' @name show.TDPClusters
#' @aliases show.TDPClusters
#' @description \code{show} method for class \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}}.
#' @param object Object of class "TDPClusters" or "TDPBrainClusters". Stores the \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}} object, usually, a result of a call to \code{\link{TDPQuery}} or \code{\link{TDPChange}}.
#' @details The summary table will be printed on console, and information on the rest of the brain will not be shown.
#' @examples
#' 
#' # four ways to display summary table (without rest):
#' tdpclusters
#' show(tdpclusters)
#' print(tdpclusters)
#' summary(tdpclusters)
#' 
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' 
setMethod("show", "TDPClusters", function(object) print(summary(object)))


#' @title Operations on cluster list
#' @name length.TDPClusters
#' @aliases length.TDPClusters
#' @description \code{length} method for class \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}}.
#' @param x Object of class "TDPClusters" or "TDPBrainClusters". Stores the \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}} object, usually, a result of a call to \code{\link{TDPQuery}} or \code{\link{TDPChange}}.
#' @details Some operations on \code{clusterlist} for class \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}}.
#' @examples
#' 
#' # get number of clusters
#' length(tdpclusters)
#' 
#' # two ways to get all clusters:
#' tdpclusters[1:length(tdpclusters)]
#' tdpclusters@clusterlist
#' 
#' # get indices of the largest cluster
#' tdpclusters[[1]]
#' 
#' # get indices of the top 2 largest clusters
#' tdpclusters[1:2]
#' 
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' 
setMethod("length", "TDPClusters", function(x) length(x@clusterlist))
setMethod("[", "TDPClusters", function(x, i, j, ..., drop=TRUE) {
  lapply(x@clusterlist[i], function(clus) clus+1)
})
setMethod("[[", "TDPClusters", function(x, i, j, ...) x@clusterlist[[i]]+1)
#setMethod("[", c("TDPClusters", "numeric", "missing", "ANY"), function(x, i, j, ..., drop=TRUE) x@clusterlist[i])
#setMethod("[[", c("TDPClusters", "numeric", "missing"), function(x, i, j, ...) x@clusterlist[[i]])


# ---------- NEWLY ADDED: CHANGE CLUSTER SIZE ---------- #

#' @title Change the size of a chosen cluster based on the TDP change inquiry
#' @name TDPChange
#' @aliases TDPChange
#' @description \code{TDPChange} is a generic function used to enlarge or shrink a cluster in terms of decreasing or increasing the TDP bound.
#' @usage TDPChange(object, v, tdpchg=0.01)
#' @param object Object of class "TDPClusters" or "TDPBrainClusters". Stores a \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}} object, usually, a result of a call to \code{\link{TDPQuery}} or \code{\link{TDPChange}}.
#' @param v Object of class "numeric". Stores a node index or a vector of 3D coordinates (one-based indexing), which is used to select a cluster in the cluster list saved in \code{object}.
#' @param tdpchg Object of class "numeric". Stores an expected change in the TDP bound, which should be within the range of \eqn{(-1,0) \cup (0,1)}. A positive value indicates the TDP bound is increased (i.e., cluster size is decreased), and a negative value indicates the TDP bound is decreased (i.e., cluster size is increased). \code{tdpchg = 0.01} by default.
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
#' # (3) change query: change the size of chosen cluster based on request
#' tdpchanges <- TDPChange(tdpclusters, v = 1, tdpchg =  0.01)  # increase TDP bound / decrease cluster size
#' tdpchanges <- TDPChange(tdpchanges,  v = 1, tdpchg =  0.01)  # further update the chosen cluster
#' summary(tdpchanges)
#' 
#' tdpchanges <- TDPChange(tdpchanges,  v = 1, tdpchg = -0.01)  # decrease TDP bound / increase cluster size
#' summary(tdpchanges)
#'     
#' @return Returns a \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}} object of cluster list.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import 
#' @export
#'
setGeneric("TDPChange", function(object, v, tdpchg) standardGeneric("TDPChange"))
setMethod("TDPChange", "TDPClusters", function(object, v, tdpchg=0.01) {
  
  # check for inputs
  if (missing(object)) stop("Missing argument 'object'")
  if (missing(v)) stop("Missing argument 'v'")
  
  # # check for tdpchg
  # maxtdp  <- object@aricluster@tdps[object@aricluster@stcs[length(object@aricluster@stcs)]+1]
  # mintdp  <- object@aricluster@tdps[object@aricluster@stcs[1]+1]
  # currtdp <- object@aricluster@tdps[object@clusterlist[[iclus+1]][length(object@clusterlist[[iclus+1]])]+1]
  # if (tdpchg<=-1 || tdpchg==0 || tdpchg>=1) stop("'tdpchg' must be non-zero & within (-1,1)")
  # if (mintdp==currtdp || maxtdp==currtdp) stop("No further changes can be attained")
  # if (tdpchg<0 && mintdp-currtdp>tdpchg) 
  #   stop("A further TDP reduction of ", abs(tdpchg), " cannot be achieved as min(TDP) = ", mintdp, " and current TDP is ", currtdp)
  # if (tdpchg>0 && maxtdp-currtdp<tdpchg) 
  #   stop("A further TDP augmentation of ", tdpchg, " cannot be achieved as max(TDP) = ", maxtdp, " and current TDP is ", currtdp)
  
  # update the chosen cluster
  clusterlist <- changeQuery(sum(object@aricluster@indexp<=v)-1, 
                             tdpchg,
                             object@aricluster@stcs,
                             object@aricluster@sizes,
                             object@aricluster@marks,
                             object@aricluster@tdps,
                             object@aricluster@childs,
                             object@clusterlist)
  
  n <- length(clusterlist)
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
  
  # if (length(unlist(clusterlist))==length(unlist(object@clusterlist)))
  #   warning("No further changes can be attained")
  
  out <- new("TDPClusters",
             aricluster = object@aricluster,
             threshold = object@threshold,
             clusterlist = clusterlist)
  
  return(out)
})
