# @title get_array
# @description get the array from a map parameter, and make compatibility checks.
# @param map array, character giving the nii file (address and) name.
# @param map_dims vector with the 3 dimension the map must agree. if \code{NULL} (default) no checks are made.
get_array <- function(map, map_dims=NULL){
  if(is.null(map)){
    if(is.null(map_dims)) {
      stop("The dims of map are not defined")
      } else {
        map = array(TRUE, map_dims)
        return(map)
      }
    }
  
  if(is.character(map))
    map = RNifti::readNifti(map)
  
  if(!is.null(map_dims)) 
    if(any(dim(map)!=map_dims)) 
      stop("The dims of map: ", paste(dim(map), sep="x"), 
           " don't fit the dims of map_dims: ", paste(map_dims, sep="x"))
  map
}