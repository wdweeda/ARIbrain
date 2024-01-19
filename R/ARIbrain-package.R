#' @name ARIbrain-package
#' @aliases ARIbrain-package
#' @title All-Resolutions Inference
#' @description  It performs All-Resolutions Inference (ARI) on fMRI data. As a main feature, it estimates lower bounds for the proportion of active voxels (true discovery proportion, TDP) in a set of clusters as, for example, given by a cluster-wise analysis. Additionally, given certain TDP thresholds, ARICluster analysis can quickly identify maximal clusters using ARI and an efficient adaptive cluster thresholding algorithm.
#' @author all of us
#' @docType _PACKAGE
#' @import hommel plyr RNifti methods fastcluster Rcpp (>= 1.0.7)
#' @importFrom stats cutree dist hclust qnorm
#' @examples
#'
#' # (1) ARI analysis
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package="ARIbrain")
#' cluster_name <- system.file("extdata", "cluster_th_3.2.nii.gz", package="ARIbrain")
#' zstat_name <- system.file("extdata", "zstat.nii.gz", package="ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package="ARIbrain")
#' 
#' ARI(Pmap = pvalue_name, clusters= cluster_name, 
#'     mask=mask_name, Statmap = zstat_name)
#' 
#' # (2) ARICluster analysis
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#'
#' # create an ARIBrainCluster object
#' aricluster <- ARIBrainCluster(Pmap = pvalue_name, mask = mask_name)
#' # answer query: find all maximal clusters given a TDP threshold
#' tdpclusters <- TDPQuery(aricluster, threshold = 0.7)
#'
NULL
