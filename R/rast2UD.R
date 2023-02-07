#' Convert a RasterLayer object to a UD object
#'
#' @param r RasterLayer object
#'
#' @return UD object retaining cell values of original RasterLayer object at original grid cell locations.
#' @export
#'
#' @examples
rast2UD <- function (r) {
  x <- raster::xFromCol(r, 1:ncol(r))
  y <- raster::yFromRow(r, nrow(r):1)
  z <- t(raster::as.matrix(r)[nrow(r):1, ])
  rproj <- raster::crs(r)
  return(list(x1 = x, x2 = y, fhat = z, crs=rproj))
}
