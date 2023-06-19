#' Given a Utilization Distribution, generate polygon(s) corresponding to specified probability contours
#'
#' @param UD object representing a Utilization Distribution. A list containing 3 elements: x (x1) and y (x2) coordinates, and corresponding density estimates (fhat)
#' @param sproj projected coordinate system of input UD
#' @param probs probability contour(s) at which to generate the polygon
#' @param spType target output geometry type (sf object). Default is 'poly', otherwise 'line'.
#'
#' @return sf polygons object representing isopleth contours of input UD
#' @export
#'
UD2sf = function(UD, sproj = NULL, probs = 0.95, spType = 'poly') {

  isopleth <- plevel <- x <- NULL

  # derive isopleth contours
  p <- as.vector(calcHR(UD, p = probs, silent = T)$ps)
  if(any(probs == 1)) p[probs == 1] <- 1e-13
  fhat.contlines = lapply(p, function(i) grDevices::contourLines(x = UD$x1, y = UD$x2, z = UD$fhat, levels = i))

  # convert to sf multipolygons object
  sldf <- lapply(fhat.contlines, function(contlines, spType = spType, sproj = sproj) {
    suppressMessages(try(sf::st_as_sf(do.call('c', lapply(1:length(contlines), function(i) {
          if(spType == 'poly') {
            return(sf::st_sfc(sf::st_polygon(list(cbind(x = contlines[[i]]$x, y = contlines[[i]]$y))), crs = sproj))
          } else {
            return(sf::st_sfc(sf::st_linestring(cbind(x = contlines[[i]]$x, y = contlines[[i]]$y)), crs = sproj))
          }
    }))), silent = T))
  })

  spatpol <- do.call(rbind, lapply(1:length(sldf), function(i) {
    if(!inherits(sldf[[i]], 'try-error')) return(dplyr::bind_cols(sldf[[i]], isopleth = probs[i], plevel = p[i]))
  })) %>%
    dplyr::filter(sf::st_is_valid(.)) %>%
    dplyr::group_by(isopleth, plevel) %>%
    dplyr::summarize(.groups = 'drop') %>%
    dplyr::rename(geometry = x)

  return(spatpol)

}
