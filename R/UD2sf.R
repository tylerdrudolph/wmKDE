#' Given a Utilization Distribution, generate polygon(s) corresponding to specified probability contours
#'
#' @param UD object representing a Utilization Distribution. A list containing 3 elements: x (x1) and y (x2) coordinates, and corresponding density estimates (fhat)
#' @param sproj projected coordinate system of input UD
#' @param probs probability contour(s) at which to generate the polygon
#'
#' @return
#' @export
#'
#' @examples
UD2sf = function(UD, sproj = NULL, probs = 0.95) {

  # derive isopleth contours
  p <- as.vector(calcHR(UD, p = probs, silent = T)$ps)
  if(any(p == 1)) p[p == 1] <- 1e-13
  fhat.contlines = lapply(p, function(i) grDevices::contourLines(x = UD$x1, y = UD$x2, z = UD$fhat, levels = i))
  # if(length(fhat.contlines) == 1) fhat.contlines <- fhat.contlines[[1]]

  # convert to sf multipolygons object
  sldf <- lapply(fhat.contlines, function(contlines) {
    suppressMessages(try(sf::st_as_sf(do.call('c', lapply(1:length(contlines), function(i) {
          sf::st_sfc(sf::st_polygon(list(cbind(x = contlines[[i]]$x, y = contlines[[i]]$y))), crs = sproj)
    }))), silent = T))
  })

  spatpol <- do.call(rbind, lapply(1:length(sldf), function(i) {
    if(!inherits(sldf[[i]], 'try-error')) return(dplyr::bind_cols(sldf[[i]], isopleth = probs[i], plevel = p[i]))
  })) %>%
    dplyr::filter(sf::st_is_valid(.)) %>%
    dplyr::group_by(isopleth, plevel) %>%
    dplyr::summarize(., .groups = 'drop') %>%
    dplyr::rename(geometry = x)

  return(spatpol)

}
