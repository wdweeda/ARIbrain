#' @title Answer queries with TDP thresholds.
#' @name TDPquery
#' @aliases TDPquery
#' @description Given a TDP threshold, it quickly finds all maximal supra-threshold clusters, each with the TDP not smaller than the threshold.
#' @param gamma A TDP threshold for forming maximal clusters.
#' @param writeNifti logical; \code{writeNifti = FALSE} by default. If \code{TRUE}, write cluster image in Nifti format ("ARIcluster.nii.gz").
#' @param rest logical; information for the rest of the brain is not shown in the summary table if \code{FALSE} (by default).
#' @param silent logical; print summary table on screen if \code{silent = FALSE} by default.
#' @examples
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package = "ARIbrain")
#' mask_name <- system.file("extdata", "mask.nii.gz", package = "ARIbrain")
#' 
#' # (1) create an ARIcluster object
#' ariclstr <- ARIcluster(Pmap = pvalue_name, mask = mask_name)
#' show(ariclstr)
#'
#' # (2) answer queries: find all maximal clusters given a TDP threshold
#' res <- TDPquery(ariclstr, gamma = 0.7)
#' res
#'     
#' @return Returns a \code{matrix} reporting Size, FalseNull, TrueNull, ActiveProp and other summary statistics for each found cluster.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import plyr, hommel, RNifti
#' @export
#' 
TDPquery <- function(ARIcluster, gamma, 
                     writeNifti=FALSE, rest=FALSE, silent=TRUE) {
  
  # check if TDP threshold is provided
  if (missing(gamma)) stop("Missing TDP threshold gamma.")
  
  # find all maximal STCs
  cluster_list  <- answerQuery(ARIcluster@m, gamma, 
                               ARIcluster@stcs, ARIcluster@sizes, 
                               ARIcluster@marks, ARIcluster@tdps, 
                               ARIcluster@childs)
  n <- length(cluster_list)
  if (n==0) stop("No clusters were found for gamma=", gamma)
  # sort clusters by descending cluster size.
  # R's built-in order function defaults to counting sort for integer vectors 
  # of length >= 200 and range < 100000, and uses radix sort when the range 
  # condition is not met. We therefore estimate whether radix sort will be 
  # output-sensitive, and otherwise ensure a counting sort is used, which is 
  # output-sensitive.
  if (n>1) {
    cluster_sizes <- sapply(cluster_list, length)
    d             <- diff(range(cluster_sizes))
    maxid         <- max(cluster_sizes)
    if (n<200 || d<100000 || n*log2(d)<=sum(cluster_sizes)) {
      cluster_list <- cluster_list[order(cluster_sizes, decreasing=TRUE)]
    } else {
      cluster_ord  <- counting_sort(n, maxid, cluster_sizes)
      cluster_list <- cluster_list[cluster_ord+1]
    }
  }
  
  # write out result nifti file
  if (writeNifti) {
    # set cluster labels & output cluster image
    for (i in 1:n) {
      ARIcluster@img_clus[ARIcluster@indexp[cluster_list[[i]]+1]+1] <- n+1-i
    }
    if (length(ARIcluster@tmp)>0)
      RNifti::writeNifti(ARIcluster@img_clus, "ARIcluster.nii.gz", template=RNifti::readNifti(ARIcluster@tmp)) 
    else {
      RNifti::writeNifti(ARIcluster@img_clus, "ARIcluster.nii.gz")
      warning("Default image parameters are used for writing 'ARIcluster.nii.gz'.")
    }
    # clear cluster labels
    for (i in 1:n) {
      ARIcluster@img_clus[ARIcluster@indexp[cluster_list[[i]]+1]+1] <- 0
    }
  }
  
  # apply summaries to each cluster (and all the rest in an extra cluster)
  out <- plyr::laply(1:n, function(i) {
    ix   <- ARIcluster@indexp[cluster_list[[i]]+1]
    xyz  <- ids2xyz(ix, dim(ARIcluster@img_clus))
    pval <- ARIcluster@sortp[cluster_list[[i]]+1]
    # ( x, y, z, max(Z) )
    cluster_xyzs <- cbind(dim1=xyz[,1], dim2=xyz[,2], dim3=xyz[,3], 
                          Stat=-qnorm(pval))
    
    clus_size <- length(cluster_list[[i]])
    clus_tdp  <- ARIcluster@tdps[cluster_list[[i]][clus_size]+1]
    unlist(c(Size=clus_size, 
             FalseNull=round(clus_size*clus_tdp), 
             TrueNull=round(clus_size*(1-clus_tdp)), 
             ActiveProp=clus_tdp,
             summary_cluster(cluster_xyzs, summary_stat="max")[-1]))
  })
  if (is.null(dim(out))) out <- t(as.matrix(out))
  
  # modify output summary table
  if (rest) {
    ix_rest   <- ARIcluster@indexp[-(unlist(cluster_list)+1)]
    ord_rest  <- (1:ARIcluster@m)[-(unlist(cluster_list)+1)]
    xyz_rest  <- ids2xyz(ix_rest, dim(ARIcluster@img_clus))
    pval_rest <- ARIcluster@sortp[ord_rest]
    rest_xyzs <- cbind(dim1=xyz_rest[,1], dim2=xyz_rest[,2], dim3=xyz_rest[,3],
                       Stat=-qnorm(pval_rest))
    
    hom       <- hommel::hommel(ARIcluster@sortp, simes=TRUE)
    rest_disc <- hommel::discoveries(hom, ix=ord_rest, alpha=ARIcluster@alpha, 
                                     incremental=FALSE)
    rest_size <- length(ix_rest)
    
    out <- rbind(out, c(rest_size,
                        rest_disc,
                        rest_size-rest_disc,
                        rest_disc/rest_size, 
                        summary_cluster(rest_xyzs, summary_stat="max")[-1]))
    rownames(out) <- c(paste0("cl", n:1), "rest")
  } else {
    rownames(out) <- paste0("cl", n:1) 
  }
   
  # attr(out, "call") = called
  if (!silent) print(out)
  
  out
}