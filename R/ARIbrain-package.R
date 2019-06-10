#' @description  It performs All-Resolutions Inference on fMRI data. As a main feature, 
#' it estimates lower bounds for the proportion of active voxels in a set of clusters as, for example, given by a cluster-wise analysis.
#' @author all of us
#' @docType package
#' @name ARIbrain-package
#' @title All-Resolutions Inference
#' @import hommel
#' @importFrom stats cutree dist hclust qnorm
#' @examples 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package="ARIbrain")
#' cluster_name <- system.file("extdata", "cluster_th_3.2.nii.gz", package="ARIbrain")
#' zstat_name <- system.file("extdata", "zstat.nii.gz", package="ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package="ARIbrain")
#' 
#' ARI(Pmap = pvalue_name, clusters= cluster_name, 
#'     mask=mask_name, Statmap = zstat_name)
#' 
NULL