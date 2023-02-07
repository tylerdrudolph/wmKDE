#' Calculate the percentage of relocations situated within a given polygon.
#'
#' @param xy 2-column data.frame or matrix containing projected x and y coordinates
#' @param poly  a SpatialPolygons or SpatialPolygonsDataFrame object
#' @param w optional vector of weights equal in length to nrow(xy)
#'
#' @return
#' @export
#'
#' @examples
percent.pts.in.poly = function(xy, poly, w=NULL) {
  inpoly = sp::over(x=sp::SpatialPointsDataFrame(sp::coordinates(xy[,c("x","y")]), data=xy, proj4string=sp::CRS(sp::proj4string(poly))), y=poly)
  inpoly[is.na(inpoly)] <- 0
  if(is.null(w)) {
    return(sum(inpoly) / length(inpoly))
  } else {
    if(all(inpoly==1)) return(1) else return(as.vector(weights::wpct(x=inpoly, weight=w)[2]))
  }
}
