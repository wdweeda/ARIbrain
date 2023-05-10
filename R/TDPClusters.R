#' @title Class "TDPClusters" for outputting cluster information
#' @name TDPClusters
#' @docType class
#' @aliases TDPClusters
#' @description The class \code{TDPClusters} is the output of a call to \code{\link{TDPQuery}}. It stores the resulting cluster information.
#' @slot aricluster Object of class "ARICluster-class". Stores the \code{\link{ARICluster-class}} object, usually, a result of a call to \code{\link{ARICluster}}.
#' @slot dims Object of class "integer". Stores the dimensions of the input p-values.
#' @slot threshold Object of class "numeric". Stores the TDP threshold.
#' @slot clusterlist Object of class "list". Stores a list of found clusters, each including indices of nodes within that cluster. Here, the node index starts from 0. Please use \code{aricluster@indexp[tdpclusters@clusterlist[[i]]+1]} to access voxel indices in 3D space for the ith largest cluster.
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
#' @slot object Object of class "TDPClusters". Stores the \code{\link{TDPClusters}} object, usually, a result of a call to \code{\link{TDPQuery}}.
#' @slot rest Object of class "logical". By default, \code{rest = FALSE} indicates that information on the rest of the brain will not be shown in the summary table.
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
      warning("No clusters were found for the rest of the brain.")
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
#' @slot object Object of class "TDPClusters". Stores the \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}} object, usually, a result of a call to \code{\link{TDPQuery}}.
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
#' @slot x Object of class "TDPClusters". Stores the \code{\link{TDPClusters}} or \code{\link{TDPBrainClusters}} object, usually, a result of a call to \code{\link{TDPQuery}}.
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
