#' Invert the values of a raster
#'
#' @param x RasterLayer object
#'
#' @return
#' @export
#'
#' @examples
raster.invert <- function (x) {
  if (!inherits(x, "RasterLayer"))
    stop("MUST BE RasterLayer OBJECT")
  rmax <- raster::cellStats(x, stat = "max", na.rm = TRUE,
                            asSample = FALSE)
  rmin <- raster::cellStats(x, stat = "min", na.rm = TRUE,
                            asSample = FALSE)
  return(((x - rmax) * -1) + rmin)
}
