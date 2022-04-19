#' @title Class "ARICluster" for storing results of adaptive thresholding algorithm
#' @name ARICluster-class
#' @docType class
#' @aliases ARICluster-class
#' @slot m Object of class "integer". Stores the total number of nodes.
#' @slot alpha Object of class "numeric". Stores the significance level.
#' @slot p Object of class "numeric". Stores input p-values in the original order.
#' @slot ordp Object of class "integer". Stores sorted orders for all input p-values.
#' @slot stcs Object of class "integer". Stores the representative nodes (sorting ranks that start from 0) for all admissible supra-threshold clusters in the ascending order of TDP bounds.
#' @slot sizes Object of class "integer". Stores the subtree sizes (cluster extents) for all nodes in the ascending order of p-values.
#' @slot tdps Object of class "numeric". Stores the TDP lower confidence bounds for all supra-threshold clusters in the ascending order of p-values.
#' @slot childs Object of class "list". Stores a list of vectors of child nodes (sorting ranks that start from 0) for all nodes in the ascending order of p-values.
#' @slot marks Object of class "integer". Initializes a vector to mark nodes within found clusters.
#' @description The class ARICluster is the output of a call to \code{\link{ARICluster}}. It stores the information needed for the next answering-query step, where all maximal supra-threshold clusters can be quickly found given a TDP threshold.
#' @return Returns an \code{\link{ARICluster-class}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel
#' @export

setClass("ARICluster",
         #representation(
         slots = list(
           m = "integer",       # stores number of p-values
           alpha = "numeric",   # stores significance level
           p = "numeric",       #(m) stores input unsorted p-values
           ordp = "integer",    #(m) stores sorted orders
           stcs = "integer",    #(m*) stores representative nodes for all admissible STCs (m*<=m)
           sizes = "integer",   #(m) stores subtree sizes for all nodes
           tdps = "numeric",    #(m) stores TDP lower bounds for all STCs, each represented by a node
           childs = "list",     #(m) stores children for all nodes
           marks = "integer"    #(m) stores marks for nodes within found clusters
         )
)
