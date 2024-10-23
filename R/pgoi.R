#' Compute the probabilistic general overlap index as per McIntyre & Barry (2018). “A Lattice-Based Smoother for Regions With Irregular Boundaries and Holes.” Journal of Computational and Graphical Statistics 27, no. 2 (2018): 360–67. https://doi.org/10.1080/10618600.2017.1375935.
#'
#' @param r A spatRaster object consisting of exactly 2 layers, the values of which are to be compared
#'
#' @return numeric pgoi
#' @import terra
#' @export
#'
pgoi <- function(r) {
  stopifnot(nlyr(r) == 2)
  ((terra::nlyr(r) - terra::global(terra::app(r, fun = 'max'), 'sum')) / (terra::nlyr(r) - 1))[1,1]
}
