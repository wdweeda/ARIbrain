#' @title All-resolutions inference (ARI) for cluster thresholding
#' @name ARIcluster
#' @aliases ARIcluster
#' @description Based on ARI, the adaptive thresholding algorithm is employed to answer queries.
#' @param Pmap 3D array of p-values or a (character) nifti file name.
#' @param mask 3D array of logicals (i.e., \code{TRUE}/\code{FALSE} for in/out of the brain). Alternatively, it may be a (character) nifti file name. If \code{mask = NULL}, it is assumed that none of the voxels have to be excluded.
#' @param conn Connectivity criterion, which can be chosen among 6, 18 (the default), and 26 for face, edge, and vertex connectivity.
#' @param alpha Significance level. \code{alpha = 0.05} by default.
#' @examples 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#' 
#' print(pvalue_name)
#' print(mask_name)
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

ARIcluster <- function(Pmap, 
                       mask=NULL, conn=18, alpha=0.05) {
  
  # check if connectivity criterion is correctly specified 
  if (!any(conn==c(6, 18, 26)))
    stop("3D connectivity criterion: 6 (face), 18 (edge), or 26 (vertex).")
  conn <- as.integer(conn)
  
  # get, fix, check parameters inconsistencies
  if(is.character(Pmap))
    tmp <- Pmap
  else {
    tmp <- character(0)
    warning("Reference Nifti file name (Pmap) is not provided for writing result cluster image.")
  }
  Pmap <- get_array(Pmap)
  mask <- get_array(mask, map_dims=dim(Pmap))
  
  # get voxel indices of the mask in 3D space
  indexp <- (1:length(mask))[mask!=0]
  #indexp <- which(mask!=0)
  
  # compute image dimensions & size of the multiple testing problem
  dims   <- dim(Pmap)
  m      <- length(indexp)
  # sort p-values of in-mask voxels in ascending order
  p      <- Pmap[indexp]
  ord    <- order(p)
  sortp  <- p[ord]
  # voxel indices of sorted p-values in 3D image (starts from 0)
  indexp <- indexp[ord] - 1
  indexp <- as.integer(indexp)
  # create 3D whole-brain mask of sorting ranks (maskI)
  maskI  <- array(0, dims)
  maskI[indexp+1] <- 1:m
  maskI  <- as.integer(maskI)
  
  # perform hommel
  hom <- hommel::hommel(p, simes=TRUE)
  tdp <- hommel::tdp(hom)
  if (tdp==0) 
    stop("There are no activations in the image!")
  else
    cat("The", 1-alpha, "confidence lower bound for the TDP for all in-mask voxels is", tdp, "\n")
  # compute h(alpha)
  halpha      <- findHalpha(hom@jumpalpha, alpha, m)
  simeshalpha <- hom@simesfactor[halpha+1]
  
  #
  # run adaptive thresholding algorithm
  #
  sizes  <- as.integer(rep( 1,m))  # subtree size for all in-mask voxels
  roots  <- as.integer(rep(-1,m))  # -1 for non-roots & 1 for roots
  # find all admissible STCs & compute TDP bounds
  adj    <- findAdjList(maskI, indexp, dims, m, conn)
  childs <- findClusters(adj, sizes, roots, m, conn)
  tdps   <- forestTDP(m, halpha, alpha, simeshalpha, sortp, sizes, roots, childs)
  stcs   <- queryPreparation(m, roots, tdps, childs)
  
  # initialize marks vector & cluster image: >0 for voxels in clusters
  marks    <- as.integer(rep(0, m))
  img_clus <- array(0, dims)
  
  out <- new("ARIcluster",
             m = m,
             alpha = alpha,
             indexp = indexp,
             sortp = sortp,
             stcs = stcs,
             sizes = sizes,
             tdps = tdps,
             childs = childs,
             marks = marks,
             img_clus = img_clus,
             tmp = tmp)
  
  return(out)
}
