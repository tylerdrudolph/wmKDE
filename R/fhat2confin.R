#' Convert UD density estimates into relative probabilities
#'
#' @param r spatRaster corresponding to UD kernel density values
#'
#' @return spatRaster correponding to kernel probability isopleths
#' @export
#'
fhat2confin <- function(r, ...) {
  
  z <- t(terra::as.matrix(r, wide = TRUE)[nrow(r):1, ])
  
  dimxy = dim(z)
  if(sum(z, na.rm = T) != 1) z = z / sum(z, na.rm = T)
  index  = 1:length(z)
  vord   = z[order(z, decreasing = TRUE)]
  indord = index[order(z, decreasing = TRUE)]
  vsu    = cumsum(vord)
  vreord = vsu[order(indord)] * 100
  
  x <- expand.grid(x = terra::xFromCol(r, 1:ncol(r)), 
                   y = terra::yFromRow(r, nrow(r):1))
  x$z = 100 - as.vector(matrix(vreord, ncol = dimxy[2])) 
       
  return(rast(x, crs = terra::crs(r)))
}

# fhat2confin <- function(z) {
#   dimxy = dim(z)
#   if(sum(z, na.rm = T) != 1) z = z / sum(z, na.rm = T)
#   index  = 1:length(z)
#   vord   = z[order(z, decreasing = TRUE)]
#   indord = index[order(z, decreasing = TRUE)]
#   vsu    = cumsum(vord)
#   vreord = vsu[order(indord)] * 100
#   return(matrix(vreord,ncol = dimxy[2]))
# }