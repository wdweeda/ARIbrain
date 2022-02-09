#' @title Class "ARIcluster" for storing results of adaptive thresholding algorithm
#' @name ARIcluster-class
#' @docType class
#' @aliases ARIcluster-class
#' @slot m Object of class "integer". Stores the number of in-mask voxels.
#' @slot alpha Object of class "numeric". Stores the significance level.
#' @slot indexp Object of class "integer". Stores the voxel indices of sorted p-values.
#' @slot sortp Object of class "numeric". Stores sorted p-values for all in-mask voxels.
#' @slot stcs Object of class "integer". Stores the representative node indices of all admissible supra-threshold clusters.
#' @slot sizes Object of class "integer". Stores the subtree sizes (cluster extents) for all nodes.
#' @slot tdps Object of class "numeric". Stores the TDP lower confidence bounds for all supra-threshold clusters.
#' @slot childs Object of class "list". Stores a list of vectors of child nodes for all nodes.
#' @slot marks Object of class "integer". Initializes a vector to mark nodes within found clusters.
#' @slot img_clus Object of class "array". Initializes a cluster array for outputting nifti file.
#' @slot tmp Object of class "character". Stores the reference nifti file name for writing cluster image in niti format. If empty, default image parameters will be used (see Details in \code{?RNifti::asNifti}).
#' @description The class ARIcluster is the output of a call to \code{\link{ARIcluster}}. It stores the information needed to quickly answer queries, i.e., find all maximal supra-threshold clusters given a TDP threshold.
#' @examples 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#' 
#' # (1) create an ARIcluster object
#' ariclstr <- ARIcluster(Pmap = pvalue_name, mask = mask_name)
#' show(ariclstr)
#'
#' # (2) answer queries: find all maximal clusters given a TDP threshold
#' res <- TDPquery(ariclstr, gamma = 0.7)
#' res
#'     
#' @return Returns an \code{\link{ARIcluster-class}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel, RNifti
#' @export

setClass("ARIcluster",
         representation(
           m = "integer",       # stores number of in-mask voxels
           alpha = "numeric",   # stores significance level
           indexp = "integer",  #(m) stores voxel indices of sorted p-values
           sortp = "numeric",   #(m) sorted p-values
           stcs = "integer",    #(m*) stores representative voxel index of each admissible STC (m*<=m)
           sizes = "integer",   #(m) stores subtree size, i.e. cluster extent, for each voxel
           tdps = "numeric",    #(m) stores TDP lower bounds for all STCs, each represented by a voxel
           childs = "list",     #(m) stores child nodes for each voxel
           marks = "integer",   #(m) marks voxels within found clusters
           img_clus = "array",  #(dim1*dim2*dim3) cluster array, with positive cluster labels
           tmp = "character"    # stores reference nifti file name for writing cluster image in nifti format
         )
)

setMethod("show", "ARIcluster", function(object) {
  cat("An ARIcluster object for", object@m, "in-mask voxels.\n")
  cat("Use TDPquery() to access this object.\n")
})

