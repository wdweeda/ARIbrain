#' @include ARICluster.R
#' @title Class "ARIBrainCluster-class" for storing neuroimaging results of adaptive thresholding algorithm
#' @name ARIBrainCluster-class
#' @docType class
#' @aliases ARIBrainCluster-class
#' @slot dims Object of class "integer". Stores the 3D image dimensions.
#' @slot indexp Object of class "integer". Stores the 3D in-mask voxel indices in the original order.
#' @description The \code{ARIBrainCluster-class} class is the output of a call to \code{\link{ARIBrainCluster}}. It inherits from the base class \code{ARICluster-class}, with the additional slots shown below, and stores the information needed for the next answering-query step, where all maximal supra-threshold clusters can be quickly found given a TDP threshold. 
#' @return Returns an \code{\link{ARIBrainCluster-class}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel, RNifti
#' @export

setClass("ARIBrainCluster",
         contains = "ARICluster",
         slots = list(
           dims = "integer",   #(3) stores 3D image dimensions
           indexp = "integer"  #(m) stores 3D voxel indices in the original order
         )
)


#' @title All-resolutions inference for cluster thresholding in neuroimaging (ARIBrainCluster)
#' @name ARIBrainCluster
#' @aliases ARIBrainCluster
#' @description \code{ARIBrainCluster} is specially designed for brain imaging data analysis.
#' @usage ARIBrainCluster(Pmap, mask, conn = 18, alpha = 0.05)
#' @param Pmap 3D array of p-values, or a (character) nifti file name.
#' @param mask 3D numeric/logical array. Non-zero/TRUE values in the mask indicate voxels in the brain. Voxels with missing values are treated as out of the brain. Alternatively, it can be a (character) nifti file name. If \code{mask} is not specified, it is assumed that none of the voxels have to be excluded.
#' @param conn Connectivity criterion (numeric or character). Can be chosen among 6 (or "face"), 18 (or "edge"), and 26 (or "vertex") for face, edge, and vertex connectivity. The default is using edge connectivity.
#' @param alpha Significance level. \code{alpha = 0.05} by default.
#' @examples 
#' 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#' 
#' # create an ARIBrainCluster-class object
#' aricluster <- ARIBrainCluster(Pmap = pvalue_name, mask = mask_name)
#' 
#' @return Returns an \code{\link{ARIBrainCluster-class}} object, based on the base class \code{\link{ARICluster-class}}.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import hommel, RNifti
#' @export
#' 
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
  # check for mask
  if (missing(mask)) mask <- array(TRUE, dim=dim(Pmap))
  if (is.character(mask)) mask <- RNifti::readNifti(mask)
  if (!is.logical(mask) && !is.numeric(mask)) stop("The mask must be logical or numeric.")
  # check for Pmap & mask
  if (!all(dim(Pmap)==dim(mask))) stop("The dims of p-value map & mask don't match.")
  # check for missing values
  if (anyNA(mask)) {
    warning("Found missing values in the mask. Those NA values will be masked out.")
    mask[is.na(mask)] <- FALSE
  }
  if (anyNA(Pmap)) {
    warning("Found missing values in p-value map. Those NA values will be masked out.")
    mask[is.na(Pmap)] <- FALSE
  }
  
  # get in-mask voxel indices in 3D space (starts from 1)
  indexp <- (1:length(mask))[mask!=0]  #indexp <- which(mask!=0)
  # extract p-values of in-mask voxels
  p      <- Pmap[indexp]
  # compute image dimensions & size of the multiple testing problem
  dims   <- dim(Pmap)
  m      <- length(indexp)
  # create 3D whole-brain mask of original orders (starts from 1)
  maskI  <- array(0, dims)
  maskI[indexp] <- seq_len(m)
  maskI  <- as.integer(maskI)
  
  # find the adjacency list
  adj <- findAdjList(maskI, as.integer(indexp-1), dims, m, conn)
  
  # # perform ARICluster
  # out <- ARICluster(p, adj, alpha=alpha)
  # 
  # # convert out to ARIBrainCluster class
  # out <- as(out, "ARIBrainCluster")
  # out@dims   <- dims
  # out@indexp <- indexp
  
  out <- new("ARIBrainCluster",
             ARICluster(p, adj, alpha=alpha),
             dims=dims,
             indexp=indexp)
  
  return(out)
}


setMethod("TDPQuery", "ARIBrainCluster", function(aricluster, threshold) {
  # tdpclusters <- callNextMethod()
  # return(tdpclusters)
  out <- new("TDPBrainClusters",
             callNextMethod())
  return(out)
})


#' @title Write out cluster image
#' @name writeClusters
#' @aliases writeClusters
#' @description \code{writeClusters} is defined for \code{\link{TDPBrainClusters}} object to write cluster image in Nifti format.
#' @usage writeCluster(tdpclusters, file, template, ...)
#' @param tdpclusters A \code{\link{TDPBrainClusters}} object, usually, a result of a call to \code{\link{TDPQuery}} or \code{\link{TDPChange}}.
#' @param file Character; The file name for outputting cluster image, where a file name that ends in .gz will be gzipped.
#' @param template A template object for writing Nifti image (please refer to \code{template} of \code{\link[RNifti]{writeNifti}} for more details).
#' @details \code{writeClusters} works similar to \code{\link[RNifti]{writeNifti}}. All its parameters but the first are relayed to \code{\link[RNifti]{writeNifti}}.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import RNifti
#' @export
#' 
writeClusters <- function(tdpclusters, file, template, ...) {
  
  # check for file
  if (missing(file)) {
    warning("File name is not defined. The default 'aribrain.nii.gz' will be used.")
    file <- "aribrain.nii.gz"
  }
  # check for template
  if (missing(template)) {
    warning("Default image parameters will be used for writing output image.")
    template <- NULL
  }
  
  # initialize cluster image: >0 for voxels within clusters
  clusimg <- array(0, tdpclusters@aricluster@dims)
  # compute number of clusters
  n <- length(tdpclusters@clusterlist)
  # set cluster labels
  for (i in 1:n) {
    clusimg[tdpclusters@aricluster@indexp[tdpclusters@clusterlist[[i]]+1]] <- n+1-i
  }
  # save output cluster image
  RNifti::writeNifti(clusimg, file=file, template=template, ...)
}
