#' @title All-resolutions inference (ARI) for cluster thresholding in neuroimaging
#' @name ARIBrainCluster
#' @aliases ARIBrainCluster
#' @description \code{\link{ARIBrainCluster}} is specially designed for brain imaging data analysis.
#' @usage ARIBrainCluster(Pmap, mask, conn = 18, alpha = 0.05)
#' @param Pmap 3D array of p-values, or a (character) nifti file name.
#' @param mask 3D numeric/logical array. Non-zero/TRUE values in the mask indicate voxels in the brain. Voxels with missing values are treated as out of the brain. Alternatively, it can be a (character) nifti file name. If \code{mask} is not specified, it is assumed that none of the voxels have to be excluded.
#' @param conn Connectivity criterion (numeric or character). Can be chosen among 6 (or "face"), 18 (or "edge"), and 26 (or "vertex") for face, edge, and vertex connectivity. The default is using edge connectivity.
#' @param alpha Significance level. \code{alpha = 0.05} by default.
#' @examples 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#' 
#' # create an ARIBrainCluster object
#' ari <- ARIBrainCluster(Pmap = pvalue_name, mask = mask_name)
#' 
#' @return Returns an \code{\link{ARIBrainCluster-class}} object, based on the base class \code{\link{ARICluster-class}}.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel, RNifti
#' @export

ARIBrainCluster <- function(Pmap, mask, conn=18, alpha=0.05) {
  
  # check for conn
  if (!is.numeric(conn) && !is.character(conn)) stop("conn should be numeric or character type.")
  if (is.numeric(conn) && !any(conn==c(6, 18, 26))) stop("3D connectivity criterion: 6 (face), 18 (edge), or 26 (vertex).")
  if (is.character(conn)) {
    conn <- match.arg(tolower(conn), c("face","edge","vertex"))
    if (conn=="face") {
      conn <- 6
    } else if (conn=="edge") {
      conn <- 18
    } else if (conn=="vertex") {
      conn <- 26
    }
  }
  conn <- as.integer(conn)
  
  # check for Pmap
  if (missing(Pmap)) stop("The p-value map is not defined.")
  if (is.character(Pmap)) Pmap <- RNifti::readNifti(Pmap)
  if (!is.numeric(Pmap)) stop("The p-value map should be numeric.")
  if (any(is.na(Pmap))) Pmap[is.na(Pmap)] <- 1
  # check for mask
  if (missing(mask)) mask <- array(TRUE, dim=dim(Pmap))
  if (is.character(mask)) mask <- RNifti::readNifti(mask)
  if (!is.logical(mask) && !is.numeric(mask)) stop("The mask must be logical or numeric.")
  if (any(is.na(mask))) mask[is.na(mask)] <- FALSE
  # check for Pmap & mask
  if (!all(dim(Pmap)==dim(mask))) stop("The dims of p-value map & mask don't match.")
  
  # get in-mask voxel indices in 3D space (starts from 1)
  indexp <- (1:length(mask))[mask!=0]  #indexp <- which(mask!=0)
  # extract p-values of in-mask voxels
  p      <- Pmap[indexp]
  # compute image dimensions & size of the multiple testing problem
  dims   <- dim(Pmap)
  m      <- length(indexp)
  # create 3D whole-brain mask of unsorted orders (starts from 1)
  maskI  <- array(0, dims)
  maskI[indexp] <- 1:m
  maskI  <- as.integer(maskI)
  
  # find the adjacency list
  adj <- findAdjList(maskI, as.integer(indexp-1), dims, m, conn)
  
  # initialize cluster image: >0 for voxels within clusters
  clusimg <- array(0, dims)
  
  # perform ARICluster
  out <- ARICluster(p, adj, alpha=alpha)
  
  # convert out to ARIBrainCluster class
  out <- as(out, "ARIBrainCluster")
  out@indexp  <- indexp
  out@clusimg <- clusimg
  
  return(out)
}
