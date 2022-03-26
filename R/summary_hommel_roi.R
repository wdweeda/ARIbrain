# @description computes alpha-level estimate of active voxels.
#   

summary_hommel_roi <- function(hommel, ix, alpha=0.05) {
  Total = length(hommel@p[ix])
  False_Null = hommel::discoveries(hommel, alpha=alpha, ix=ix)
  True_Null = Total-False_Null
  Active_Proportion = False_Null / Total
  list(Size=Total, FalseNull=False_Null, TrueNull=True_Null, ActiveProp=Active_Proportion)
}

