#' Convert a UD object to a RasterLayer object
#'
#' @param UD input utilization distribution (UD) object; a list containing three elements named 'x1', 'x2' and 'fhat'.
#' @param sproj projected coordinate system of the input UD.
#'
#' @return spatRaster object retaining original UD cell values at original xy grid cell locations.
#' @export
#'
#' @examples
#'
UD2rast <- function(UD, sproj = NULL) {

  x <- expand.grid(x = UD$x1, y = UD$x2)
  x$z <- as.vector(UD$fhat)
  return(terra::rast(x, crs = sproj$wkt))
  
  # return(terra::flip(x = terra::rast(x = t(UD$fhat),
  #                             extent = terra::ext(min(UD$x1), max(UD$x1), min(UD$x2), max(UD$x2)),
  #                             crs = sproj$wkt)))

}
