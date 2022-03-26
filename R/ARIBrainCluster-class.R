#' @include ARICluster-class.R
#' @title Class "ARIBrainCluster" for storing neuroimaging results of adaptive thresholding algorithm
#' @name ARIBrainCluster-class
#' @docType class
#' @aliases ARIBrainCluster-class
#' @slot indexp Object of class "integer". Stores the 3D in-mask voxel indices in the original order.
#' @slot clusimg Object of class "array". Initializes a cluster array for writing 3D nifti file.
#' @description The \code{ARIBrainCluster} class is the output of a call to \code{\link{ARIBrainCluster}}. It inherits from the base class \code{ARICluster}, with the additional slots shown below, and stores the information needed for the next answering-query step, where all maximal supra-threshold clusters can be quickly found given a TDP threshold. 
#' @return Returns an \code{\link{ARIBrainCluster-class}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel, RNifti
#' @export

setClass("ARIBrainCluster",
         contains = "ARICluster",
         #representation(
         slots = list(
           indexp = "integer",  #(m) stores 3D voxel indices for unsorted p-values
           clusimg = "array"    #(dim1*dim2*dim3) stores cluster array of all zeros
         )
)