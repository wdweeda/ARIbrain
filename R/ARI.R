#' @title Valid Circular Inference (ARI) for Brain Imaging
#' @description Valid Circular Inference (ARI) for Brain Imaging
#' @param Pmap 3D array of p-values or a (character) nifti file name.
#' @param clusters 3D array of cluster ids (0 when voxel does not belong to any cluster)
#'  or a (character) nifti file name.
#' @param mask 3D array of locicals (i.e. \code{TRUE}/\code{FALSE} in/out of the brain). Alternatively 
#' it may be a (character) nifti file name. If \code{mask=NULL}, it is assumed that non of the voxels have to be excluded.
#' @param alpha Significance level. \code{alpha=.05} by default.
#' @param Statmap Statistics (usually t-values) on which the summaries are based. Can be either
#' a 3D array, a (character) nifti file name or a function with argument \code{ix} used in the function to select the voxels belonging to a given cluster.
#' By default \code{Statmap = function(ix) -qnorm(Pmap[ix])} which convert the p-values in one-sided z-score.
#' @param summary_stat Choose among \code{=c("max", "center-of-mass")}.
#' @param silent \code{FALSE} by default.
#' @examples 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package="ARIbrain")
#' cluster_name <- system.file("extdata", "cluster_th_3.2.nii.gz", package="ARIbrain")
#' zstat_name <- system.file("extdata", "zstat.nii.gz", package="ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package="ARIbrain")
#' 
#' print(mask_name)
#' print(pvalue_name)
#' print(cluster_name)
#' print(zstat_name)
#' 
#' ARI(Pmap = pvalue_name, clusters= cluster_name, 
#'     mask=mask_name, Statmap = zstat_name)
#'     
#' @return A \code{matrix} reporting Size, FalseNull,  TrueNull, ActiveProp and other statistics for each cluster.
#' @export
#' @import hommel

ARI <- function(Pmap, clusters, mask=NULL, alpha=.05,Statmap=function(ix) -qnorm(Pmap[ix]),
                summary_stat=c("max", "center-of-mass"),silent=FALSE){
  
  # get, fix, check parameters inconsistencies
  Pmap = get_array(Pmap)
  clusters = get_array(clusters,map_dims=dim(Pmap))
  mask = get_array(mask,map_dims=dim(Pmap))
  if(is.function(Statmap)) {
    StatFun=Statmap
  } else {
    Statmap= get_array(Statmap,map_dims=dim(Pmap))
    StatFun <- function(ix) Statmap[ix]
    }
  
  # called=match.call()
  summary_stat=match.arg(summary_stat,c("max", "center-of-mass"))
  
  # get the indices of the mask
  mask=which(mask!=0)
  
  #perform hommel
  hom <-hommel(Pmap[mask])
  if(!silent) {temp=(summary(hom))
  cat("\n")}
  
  # define number of clusters
  clstr_id=sort(unique(as.vector(clusters[mask])),decreasing = TRUE)

  #apply summaries to each cluster (and all the rest in an extra cluster)
  out=plyr::laply(clstr_id,function(i){
    ix=clusters==i
    ix[-mask]=FALSE

    cluster_ids=which(ix,arr.ind = TRUE)
    cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
    
    unlist(c(summary_hommel_roi(hommel = hom,ix=ix[mask], alpha = alpha),
      summary_cluster(cluster_ids)[-1])
    )
  })
  if(is.null(dim(out))) out <- t(as.matrix(out))
  rownames(out)=paste("cl",sep="",clstr_id)
 
   # attr(out,"call")=called
  if(!silent) print(out)
  out
}
