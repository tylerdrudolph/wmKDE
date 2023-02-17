#' Convert a spatRaster object to a list representing a UD (utilization distribution) object (list of x1, x2, fhat)
#'
#' @param r spatRaster object
#' @param sproj projected CRS of input spatRaster
#'
#' @return UD object retaining cell values of original spatRaster object at original grid cell locations.
#' @export
#'
#' @examples
#'
rast2UD <- function (r, sproj = NULL) {

  x <- terra::xFromCol(r, 1:ncol(r))
  y <- terra::yFromRow(r, nrow(r):1)
  z <- t(terra::as.matrix(r)[nrow(r):1, ])
  return(list(x1 = x, x2 = y, fhat = z, sproj = sproj))

}
