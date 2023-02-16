#' Cell-wise weighted means of multiple Utilization Distributions
#'
#' @param UDlist list of utilization distributions of equal dimensions.
#' @param w vector of individual UD weights for cell-wise weighted means.
#' @param sproj desired output spatial resolution of mean weighted UD in projected coordinate units (m)s.
#' @param silent logical; should messages be printed?
#' @param sproj projected coordinate system of input UDs
#' @param dres desired pixel resolution of output UD
#' @param checksum logical; should adjustments be made to ensure sure that all cell volumes sum to 1?
#'
#' @return
#' @export
#'
#' @examples
wmUD = function(udList, w, sproj = NULL, silent = FALSE, dres = NULL, checksum = TRUE) {

  if(length(w) != length(udList)) stop("poids (w) doit être de la même longueur que udList")

  if(!silent) {
    if(length(unique(w))==1) {
      cat("Calcul du kernel moyen", "\n")
    } else {
      cat("Calcul du kernel pondéré", "\n")
    }
  }

  # Ensure all UDs sum to 1 and convert to raster
  if(checksum) udList <- lapply(udList, function(i) wmKDE::reweightedUD(i, out=1))

  # Convert to raster stack
  rstack <- terra::rast(lapply(udList, function(x) wmKDE::ud2rast(x, sproj = sproj)))

  pstack = terra::weighted.mean(rstack[[w > 0]], w[w > 0])

  # resample to specified resolution
  if(!is.null(dres)) {
    if(unique(round(terra::res(pstack))) != dres) {
      rastref <- terra::rast(terra::ext(pstack), crs = sproj, res = dres)
      if(!silent) message("Resampling to ", dres, "m...")
      pstack = terra::resample(x = pstack, y = rastref, method='bilinear')
    }
  }

  # ensure values sum to 1
  if(checksum) return(wmKUD::reweightedUD(wmKUD::rast2UD(pstack), out=1)) else return(wmKUD::rast2UD(pstack))

}

