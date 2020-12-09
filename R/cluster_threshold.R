#' @title cluster_threshold
#' @description Get spatially-connected clusters starting from a 3D map of logical values
#' @param map 3D map of logical values. \code{TRUE} if the voxel it to be clustered (e.g. it is supra-threshold).  
#' @param max_dist maximum distance allowed to in the same cluster. By default:
#' \code{max_dist=sqrt(3)} i.e. comprises all the voxels up to the corners souranding the target voxel. A value such as
#' \code{max_dist=sqrt(2)} excludes the corners.
#' @return a 3D map (same size of \code{map}) with integer values identifying the cluster and 0 elsewhere.
#' @examples 
#' \dontrun{
#' Tmap = RNifti::readNifti(system.file("extdata", "zstat.nii.gz", package="ARIbrain"))
#' clstr=cluster_threshold(Tmap>3.2)
#' table(clstr)
#' }
#' @export
#' @importFrom fastcluster hclust.vector

cluster_threshold <- function(map, max_dist=sqrt(3)){
  ### slower:
  # map=spmT
  # threshold=3.2
  # nmat <- expand.grid(-1:1, -1:1, -1:1)
  # nmat <- nmat[-c(1,3,7,9,14,19,21,25,27), ]
  # system.time(
  # {Suprathreshold_TF = cluster.threshold(spmT>=3.2, nmat=nmat,size.thr = .5)})
  # table(Suprathreshold_TF)
  
  #an alternative and faster way:
  map = ifelse(is.na(map), FALSE, map)
  Suprathreshold_TF=which(map,arr.ind = TRUE)
  #########
  #dd = dist(Suprathreshold_TF)
  #hc = hclust(dd, "single")
  hc =  hclust.vector(Suprathreshold_TF, "single")
  # plot(hc)
  # ct = cutree(hc,k=5)
  # pander(table(ct))
  # 
  
  ct = cutree(hc,h=max_dist)
  
  ## sort the cluster names on the basis of their size
  new_cluster_names=rank(table(ct),ties.method = "random")
  ct_new=rep(NA,length(ct))
  for(i in 1:length(new_cluster_names)){
    ct_new[ct==as.numeric(names(new_cluster_names)[i])]=new_cluster_names[i]
  }
  # table(ct,ct_new)
  # table(ct_new)
  ct=ct_new
  rm(ct_new)
  #########
  cluster_map=array(0,dim(map))
  cluster_map[map] =ct
  # pander::pander(table(cluster_map))
  # print(table(cluster_map))
  cluster_map
}
