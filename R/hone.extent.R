#' hone.extent
#' Generate a buffered Extent around a polygon or other spatial object from which an extent object can be extracted.
#'
#' @param spatpol Spatial polygon from which to generate the extent object
#' @param r expansion radius in projected coordinates (m)
#'
#' @return
#' @export
#'
#' @examples
hone.extent <- function(spatpol, r=1000) {
  exbox <- raster::extent(spatpol)
  out <- spatstat.geom::dilation(spatstat.geom::owin(xrange=exbox[1:2], yrange=exbox[3:4]), r=r)
  return(raster::extent(out$xrange[1], out$xrange[2], out$yrange[1], out$yrange[2]))
}
