#' @include TDPClusters.R
#' @title Class "TDPBrainClusters" for outputting cluster information
#' @name TDPBrainClusters
#' @docType class
#' @aliases TDPBrainClusters
#' @description The class \code{TDPBrainClusters} is the output of a call to \code{\link{TDPQuery}} or \code{\link{TDPChange}}. It stores the resulting cluster information.
#' @slot aricluster Object of class "ARIBrainCluster-class". Stores the \code{\link{ARIBrainCluster-class}} object, usually, a result of a call to \code{\link{ARIBrainCluster}}.
#' @slot threshold Object of class "numeric". Stores the TDP threshold.
#' @slot clusterlist Object of class "list". Stores a list of found clusters, each including indices of nodes within that cluster. Here, the node index starts from 0. Please use \code{aricluster@indexp[tdpclusters@clusterlist[[i]]+1]} to access voxel indices in 3D space for the ith largest cluster.
#' @return Returns a \code{\link{TDPBrainClusters}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import
#' @export
#'
setClass("TDPBrainClusters",
         contains = "TDPClusters"
) 
# setClass("TDPBrainClusters",
#          contains = "TDPClusters",
#          slots = list(
#            dims = "integer",   #(3) stores 3D image dimensions
#            indexp = "integer"  #(m) stores 3D voxel indices in the original order
#          )
# )


#' @title Summarize cluster information and generate summary table for found clusters. 
#' @name summary.TDPBrainClusters
#' @aliases summary.TDPBrainClusters
#' @description \code{summary} method for class \code{\link{TDPBrainClusters}}.
#' @usage summary(object, ..., rest = FALSE)
#' @param object Object of class "TDPBrainClusters". Stores the \code{\link{TDPBrainClusters}} object, usually, a result of a call to \code{\link{TDPQuery}} or \code{\link{TDPChange}}.
#' @param rest Object of class "logical". By default, \code{rest = FALSE} indicates that information on the rest of the brain will not be shown in the summary table.
#' @details If the output is not assigned, the summary table will be printed on console.
#' @examples
#' 
#' table <- summary(tdpclusters)      # get summary table (without rest)
#' summary(tdpclusters, rest = TRUE)  # print summary table (with rest)
#' 
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @return Returns a summary table reporting Size, TDN (lower bound), #TrueNull (upper bound), TDP (lower bound), maximum statistic and the corresponding XYZ coordinates for each found cluster.
#' @export
#'
setMethod("summary", "TDPBrainClusters", function(object, ..., rest=FALSE) {
  sumtable <- callNextMethod()
  
  # store row & column names
  names_row <- rownames(sumtable)
  names_col <- c(colnames(sumtable)[1:5],"X","Y","Z")
  # find ids of rows without NA values
  ids_row   <- which(!is.na(sumtable[,4]))
  # expand sumtable by adding xyz coordinates
  sumtable  <- cbind(sumtable, matrix(NA,dim(sumtable)[1],2))
  if (length(ids_row)>0) {
    sumtable[ids_row, 6:8] <- ids2xyz(as.integer(object@aricluster@indexp[sumtable[ids_row,6]]-1), object@aricluster@dims)
  }
  # update row & column names
  rownames(sumtable) <- names_row
  colnames(sumtable) <- names_col
  
  sumtable
})


setMethod("[", "TDPBrainClusters", function(x, i, j, ..., drop=TRUE) {
  # translate node indices in cluster list to 3D voxel indices
  lapply(x@clusterlist[i], function(clus) x@aricluster@indexp[clus+1])
})
setMethod("[[", "TDPBrainClusters", function(x, i, j, ...) {
  x@aricluster@indexp[x@clusterlist[[i]]+1]
})


# ---------- NEWLY ADDED: CHANGE CLUSTER SIZE ---------- #
setMethod("TDPChange", "TDPBrainClusters", function(object, v, tdpchg=0.01) {
  if (length(v)==2 || length(v)>3) stop("'v' must be an index or XYZ coordinates")
  if (length(v)==3) v <- (v[3]-1)*object@aricluster@dims[1]*object@aricluster@dims[2] + (v[2]-1)*object@aricluster@dims[1] + v[1]
  
  out <- new("TDPBrainClusters",
             callNextMethod())
  return(out)
})