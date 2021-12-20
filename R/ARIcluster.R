#' @title All-Resolutions Inference (ARI) for Cluster Thresholding
#' @name ARIcluster
#' @description Given a TDP threshold, it quickly finds all maximal supra-threshold clusters using ARI, each with the TDP not smaller than the threshold.
#' @param Pmap 3D array of p-values or a (character) nifti file name.
#' @param gamma A pre-specified TDP threshold for forming maximal clusters.
#' @param outname Specify nifti file name for writing out cluster array of cluster ids. By default, "ARIcluster.nii.gz" is used.
#' @param mask 3D array of logicals (i.e., \code{TRUE}/\code{FALSE} in/out of the brain). Alternatively, it may be a (character) nifti file name. If \code{mask = NULL}, it is assumed that none of the voxels have to be excluded.
#' @param Statmap Statistics (usually t-values) on which the summaries are based. Can be either a 3D array, a (character) nifti file name or a function with argument \code{ix} used in the function to select the voxels belonging to a given cluster. By default, \code{Statmap = function(ix) -qnorm(Pmap[ix])} converts the p-values to one-sided z-scores.
#' @param conn Connectivity criterion, which can be chosen among 6, 18 (by default), and 26 for face, edge, and vertex connectivity.
#' @param alpha Significance level. \code{alpha = 0.05} by default.
#' @param rest logical; information for the rest of the brain is not shown in the summary table if \code{FALSE} (by default).
#' @param silent logical; print summary table on screen if \code{FALSE} by default.
#' @examples 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' zstat_name <- system.file("extdata", "zstat.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#' 
#' print(mask_name)
#' print(pvalue_name)
#' print(zstat_name)
#' 
#' ARIcluster(Pmap = pvalue_name, mask = mask_name, Statmap = zstat_name, gamma = 0.5)
#'     
#' @return A \code{matrix} reporting Size, FalseNull, TrueNull, ActiveProp and other summary statistics for each cluster.
#' @import hommel, plyr, RNifti
#' @export

ARIcluster <- function(Pmap, gamma, outname="ARIcluster.nii.gz", 
                       mask=NULL, Statmap=function(ix) -qnorm(Pmap[ix]), 
                       conn=18, alpha=0.05, rest=FALSE, silent=FALSE) {
  
  # get, fix, check parameters inconsistencies
  Pmap_name <- Pmap
  mask_name <- mask
  Pmap <- get_array(Pmap)
  mask <- get_array(mask, map_dims=dim(Pmap))
  if(is.function(Statmap)) {
    StatFun <- Statmap
  } else {
    Statmap <- get_array(Statmap, map_dims=dim(Pmap))
    StatFun <- function(ix) Statmap[ix]
  }
  
  # check if connectivity criterion is correctly specified 
  if (!any(conn==c(6, 18, 26)))
    stop("3D connectivity criterion: 6 (face), 18 (edge), or 26 (vertex).")
  conn <- as.integer(conn)
  
  # check if TDP threshold is provided
  if (missing(gamma)) stop("Missing TDP threshold gamma.")
  
  # compute image dimensions & size of the multiple testing problem
  dims   <- dim(Pmap)
  m      <- sum(mask!=0)
  # compute sorting ranks of in-mask voxels
  p      <- Pmap[mask!=0]
  ord    <- order(p)
  sortp  <- p[ord]
  ranks  <- rep(0, m)
  ranks[ord] <- 1:m
  # voxel indices of sorted p-values in 3D image
  indexp <- which(mask!=0)
  indexp <- indexp[ord] - 1
  indexp <- as.integer(indexp)
  # create 3D whole-brain mask of sorting ranks (maskI)
  maskI  <- array(0, dims)
  maskI[mask!=0] <- ranks
  maskI  <- as.integer(maskI)
  
  # perform hommel
  hom <- hommel::hommel(p, simes=TRUE)
  tdp <- hommel::tdp(hom)
  if (tdp==0) warning("There are no activations in the image!")
  # compute h(alpha)
  halpha      <- findHalpha(hom@jumpalpha, alpha, m)
  simeshalpha <- hom@simesfactor[halpha+1]
  
  #
  # run adaptive thresholding algorithm
  #
  sizes <- as.integer(rep( 1,m))  # subtree size for all in-mask voxels
  roots <- as.integer(rep(-1,m))  # -1 for non-roots & 1 for roots
  # find all admissible STCs & compute TDP bounds
  child <- findClusters(maskI, indexp, dims, sizes, roots, m, conn)
  tdps  <- forestTDP(m, halpha, alpha, simeshalpha, sortp, sizes, roots, child)
  L     <- queryPreparation(m, roots, tdps, child)
  # generate maximal clusters under TDP threshold
  mark  <- answerQuery(m, gamma, L, sizes, tdps, child)
  
  # sort cluster ids by descending cluster size
  if (max(mark)>1) {
    clstr_id <- seq.int(max(mark), 1, -1)
    clstr_size <- sapply(clstr_id, function(i) sum(mark==i))
    clstr_id <- clstr_id[sort(clstr_size, decreasing=TRUE, index.return=TRUE)$ix]
    if (!all(clstr_id==seq.int(max(mark), 1, -1))) {
      mark_new <- integer(length(mark))
      i_new <- max(mark)
      for (i in clstr_id) {
        mark_new[mark==i] <- i_new
        i_new <- i_new-1
      }
      mark <- mark_new
    }
  }
  
  # convert results to 3D array
  clus <- rep(0, m)
  clus[ord] <- mark
  img_clus  <- array(0, dims)
  img_clus[mask!=0] <- clus
  
  # write out nifti file
  if (is.character(Pmap_name))
    RNifti::writeNifti(img_clus, outname, template=Pmap)
  else if (is.character(mask_name))
    RNifti::writeNifti(img_clus, outname, template=mask)
  else
    RNifti::writeNifti(img_clus, outname)
  
  # get the indices of the mask
  mask <- which(mask!=0)
  #mask <- as.logical(mask)
  
  # define label of clusters
  clstr_id <- sort(unique(mark), decreasing=TRUE)
  
  # apply summaries to each cluster (and all the rest in an extra cluster)
  out <- plyr::laply(clstr_id, function(i) {
    ix <- img_clus==i
    ix[-mask] <- FALSE
    #ix <- ix & mask
    
    cluster_ids <- which(ix, arr.ind=TRUE)
    cluster_ids <- cbind(cluster_ids, Stat=StatFun(ix))
    
    unlist(c(summary_hommel_roi(hommel=hom, ix=ix[mask], alpha=alpha),
             summary_cluster(cluster_ids, summary_stat="max")[-1]))
  })
  if (is.null(dim(out))) out <- t(as.matrix(out))
  
  # modify output summary table
  min_mark <- min(mark)
  max_mark <- max(mark)
  if (rest) {
    if (min_mark>0) {
      out <- rbind(out, integer(dim(out)[2]))
      rownames(out) <- c(paste0("cl", clstr_id), "rest")
    } else if (max_mark==0) {
      rownames(out) <- "rest"
    } else {
      rownames(out) <- c(paste0("cl", clstr_id[clstr_id>0]), "rest")
    }
  } else {
    if (min_mark==0 & max_mark==0) {
      out <- matrix(0, 0, dim(out)[2], dimnames=list(NULL, colnames(out)))
    } else if (min_mark==0 & max_mark==1) {
      out <- t(as.matrix(out[1,]))
      rownames(out) <- paste0("cl", clstr_id[clstr_id>0]) 
    } else if (min_mark==0 & max_mark>1) {
      out <- out[-dim(out)[1],]
      rownames(out) <- paste0("cl", clstr_id[clstr_id>0]) 
    } else {
      rownames(out) <- paste0("cl", clstr_id) 
    }
  }
  
  # attr(out, "call") = called
  if (!silent) print(out)
  
  out
}
