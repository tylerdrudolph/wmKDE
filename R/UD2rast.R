##########################################################################################
## Convertir un objet de format UD en 'RasterLayer'
##
## Suppose un syst?me de coordonn?es g?ographiques projet?es en Qu?bec Lambert NAD 1983)

#' Convert a UD object to a RasterLayer object
#'
#' @param UD input UD object; a list containing three elements named 'x1', 'x2' and 'fhat'.
#' @param dcrs projected coordinate system of the input UD.
#'
#' @return RasterLayer object retaining original UD cell values at original xy grid cell locations.
#' @export
#'
#' @examples
UD2rast = function(UD, dcrs=NULL) {
  if('crs' %in% names(UD)) dcrs <- UD$crs
  return(raster::raster(list(x=UD$x1, y=UD$x2, z=UD$fhat), crs=dcrs))
}
