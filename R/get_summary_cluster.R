# @param coord_and_values \code{array} or \code{data.frame} with the following columns: the first column contains the values (usually t-values) to be summarized, the remaining columns contain the coordinates of the statistics.
# @param summary_stat \code{=c("max", "center-of-mass")}
# @return a \code{list} with size of the cluster, coordinates of the max or the center-of-mass and its value


# coordinates finding (not the index, while the cohordinates in the space ... 
# spmP = readNifti(paste(sep="","./",data_folder,"/pvalue_stat1.nii.gz"))
# > str(spmP)
# niftiImage [1:91, 1:109, 1:91] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
# - attr(*, "pixdim")= num [1:3] 2 2 2
# - attr(*, "pixunits")= chr [1:2] "mm" "s"
# - attr(*, ".nifti_image_ptr")=<externalptr> 

summary_cluster <- function(coord_and_values, summary_stat=c("max", "center-of-mass")) {
  # compute max and/or centre of gravity, see below
  name_stat = names(coord_and_values)[1]
  summary_stat = match.arg(summary_stat, c("max", "center-of-mass"))
  out = list(Size=nrow(coord_and_values))
  if (summary_stat=="max") {
    id_max = which.max(coord_and_values[,4])
    out = c(out, coord_and_values[id_max,])
  } else if (summary_stat=="center-of-mass") {
    id_mean = colMeans(coord_and_values[, -1, drop=FALSE])
    id_closest_to_baricenter = which.min(rowSums(t(t(coord_and_values[, -1, drop=FALSE])-id_mean)^2))
    out = c(out, coord_and_values[id_closest_to_baricenter,])
    names(out)[names(out)==name_stat] = summary_stat
  }
  out
}
