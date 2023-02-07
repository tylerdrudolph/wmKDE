#' Given a Utilization Distribution, generate polygon(s) corresponding to specified probability contours
#'
#' @param UD object representing a Utilization Distribution. A list containing 3 elements: x (x1) and y (x2) coordinates, and corresponding density estimates (fhat)
#' @param PID character identifier to be written to SpatialPolygonsDataFrame output
#' @param proj projected coordinate system of input UD
#' @param probs probability contour(s) at which to generate the polygon
#' @param levels threshold minimum density value at which to generate the polygon contours. Useful to approximate the 100% probability contour polygon. Overrides the probs argument. NOTE`` that for the 100% isolines, use 'levels=1e-13' and not 'probs=1'.
#'
#' @return
#' @export
#'
#' @examples
UD2sp = function(UD, PID, proj = sp::CRS("+init=epsg:32198"), probs = 0.95, levels=NULL) {

  # generate polygon contours
  if(is.null(levels)) {
    temp.hr = calcHR(UD, p = probs)
    fhat.contlines = grDevices::contourLines(x=UD$x1, y=UD$x2, z=UD$fhat, levels=temp.hr$ps[[1]])
  } else {
    fhat.contlines = grDevices::contourLines(x=UD$x1, y=UD$x2, z=UD$fhat, levels=levels)
  }

  # convert to SpatialLinesDataFrame
  sldf <- maptools::ContourLines2SLDF(fhat.contlines)
  sp::proj4string(sldf) <- proj

  # convert to PolySet
  ps <- maptools::SpatialLines2PolySet(sldf)
  ps$PID = PID
  # attr(ps,"projection") <- "UTM" # to ensure correct format for the maptools functions; does not alter crs!
  # attr(ps,"zone") <- 18 # to ensure correct format for the maptools functions; does not alter crs!

  # convert PolySet to Spatial Polygons and redefine projection just in case...
  options(warn=-1)
  spatpol <- maptools::PolySet2SpatialPolygons(ps)
  sp::proj4string(spatpol) <- proj

  return(spatpol)

}
