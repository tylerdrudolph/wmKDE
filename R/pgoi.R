pgoi <- function(r) {
  # stopifnot(nlyr(r) == 2) 
  ((terra::nlyr(r) - terra::global(terra::app(r, fun = 'max'), 'sum')) / (terra::nlyr(r) - 1))[1,1]
}
