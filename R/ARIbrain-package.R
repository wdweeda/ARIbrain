#' @name ARIbrain-package
#' @title All-Resolutions Inference
#' @description  It performs All-Resolutions Inference (ARI) on fMRI data. As a main feature, it estimates lower bounds for the proportion of active voxels (true discovery proportion, TDP) in a set of clusters as, for example, given by a cluster-wise analysis. Additionally, it can quickly find maximal clusters using ARI under certain TDP thresholds. 
#' @author all of us
#' @docType package
#' @import hommel RNifti
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