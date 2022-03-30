#' @title Class "TDPquery" for outputting cluster information
#' @name TDPquery-class
#' @docType class
#' @aliases TDPquery-class
#' @slot gamma Object of class "numeric". Stores the TDP threshold.
#' @slot clusterlist Object of class "list". Stores a list of found clusters, each including sorting ranks for voxels within that cluster.
#' @description The class \code{TDPquery} is the output of a call to \code{\link{TDPquery}}. It stores the resulting cluster information.
#' @return Returns an \code{\link{TDPquery-class}} object.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import
#' @export

setClass("TDPquery",
         #representation(
         slots = list(
           gamma = "numeric",    # stores TDP threshold
           clusterlist = "list"  #(n) stores a list of found clusters
         )
)


#' @title Generate summary table for found clusters
#' @description It is defined for \code{\link{ARIBrainCluster-class}} & \code{\link{TDPquery-class}} objects to create the summary table. If the output is not assigned, the summary table will be printed on console..
#' @usage summaryCluster(obj, object, rest = FALSE)
#' @param obj An \code{\link{ARIBrainCluster-class}} object.
#' @param object A \code{\link{TDPquery-class}} object.
#' @param rest Logical; if \code{FALSE} (by default), information on the rest of the brain is not shown in the summary table.
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @return Returns a summary table reporting Size, FalseNull, TrueNull, ActiveProp and other summary statistics for each found cluster.
#' @import plyr
#' @export
#' 
setGeneric("summaryCluster", function(obj, object, rest) {
  standardGeneric("summaryCluster")
})
setMethod("summaryCluster",
          signature(obj="ARIBrainCluster", object="TDPquery"),
          function(obj, object, rest=FALSE) {
            
            # compute number of clusters
            n <- length(object@clusterlist)
            # apply summaries to each cluster
            sumtable <- plyr::laply(1:n, function(i) {
              ix   <- obj@indexp[obj@ordp[object@clusterlist[[i]]+1]]
              xyz  <- ids2xyz(ix, dim(obj@clusimg))
              pval <- obj@p[obj@ordp[object@clusterlist[[i]]+1]]
              # ( x, y, z, max(Z) )
              cluster_xyzs <- cbind(dim1=xyz[,1], dim2=xyz[,2], dim3=xyz[,3], Stat=-qnorm(pval))
              
              clus_size <- length(object@clusterlist[[i]])
              clus_tdp  <- obj@tdps[object@clusterlist[[i]][clus_size]+1]
              unlist(c(Size=clus_size, 
                       FalseNull=round(clus_size*clus_tdp), 
                       TrueNull=round(clus_size*(1-clus_tdp)), 
                       ActiveProp=clus_tdp,
                       summary_cluster(cluster_xyzs, summary_stat="max")[-1]))
            })
            if (is.null(dim(sumtable))) sumtable <- t(as.matrix(sumtable))
            
            # modify output summary table by adding the "rest" information
            if (rest) {
              ord_rest  <- obj@ordp[-(unlist(object@clusterlist)+1)]
              pval_rest <- obj@p[ord_rest]
              ix_rest   <- obj@indexp[ord_rest]
              xyz_rest  <- ids2xyz(ix_rest, dim(obj@clusimg))
              rest_xyzs <- cbind(dim1=xyz_rest[,1], dim2=xyz_rest[,2], dim3=xyz_rest[,3], Stat=-qnorm(pval_rest))
              
              hom       <- hommel::hommel(obj@p, simes=TRUE)
              rest_disc <- hommel::discoveries(hom, ix=ord_rest, alpha=obj@alpha, incremental=FALSE)
              rest_size <- length(ix_rest)
              
              sumtable <- rbind(sumtable, c(rest_size,
                                            rest_disc,
                                            rest_size-rest_disc,
                                            rest_disc/rest_size,
                                            summary_cluster(rest_xyzs, summary_stat="max")[-1]))
              rownames(sumtable) <- c(paste0("cl", n:1), "rest")
            } else {
              rownames(sumtable) <- paste0("cl", n:1)
            }
            
            sumtable
          })


#' @title Write out result cluster image
#' @description It is defined for \code{\link{ARIBrainCluster-class}} & \code{\link{TDPquery-class}} objects to write cluster image in Nifti format.
#' @usage writeCluster(obj, object, file = "aribrain.nii.gz", template = NULL)
#' @param obj An \code{\link{ARIBrainCluster-class}} object.
#' @param object A \code{\link{TDPquery-class}} object.
#' @param file Character; The file name \code{"aribrain.nii.gz"} is used by default, where a file name that ends in .gz will be gzipped.
#' @param template A template object for writing Nifti image (please refer to \code{template} of \code{\link[RNifti]{writeNifti}} for more details). 
#' @author Xu Chen, Thijmen Krebs, Wouter Weeda.
#' @import RNifti
#' @export
#' 
setGeneric("writeCluster", function(obj, object, file, template) {
  standardGeneric("writeCluster")
})
setMethod("writeCluster", 
          signature(obj="ARIBrainCluster", object="TDPquery"), 
          function(obj, object, file="aribrain.nii.gz", template=NULL) {
            
            # compute number of clusters
            n <- length(object@clusterlist)
            # set cluster labels
            for (i in 1:n) {
              obj@clusimg[obj@indexp[obj@ordp[object@clusterlist[[i]]+1]]] <- n+1-i
            }
            # save output cluster image
            if (is.null(template)) warning("Default image parameters are used for writing output image.")
            RNifti::writeNifti(obj@clusimg, file=file, template=template)
            if (file=="aribrain.nii.gz") cat("'aribrain.nii.gz' is the output cluster image.\n")
            # clear cluster labels
            for (i in 1:n) {
              obj@clusimg[obj@indexp[obj@ordp[object@clusterlist[[i]]+1]]] <- 0
            }
          })

