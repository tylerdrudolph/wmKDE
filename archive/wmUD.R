#' Cell-wise weighted means of multiple Utilization Distributions
#'
#' @param udList list of utilization distributions of equal dimensions.
#' @param w vector of individual UD weights for cell-wise weighted means.
#' @param sproj desired output spatial resolution of mean weighted UD in projected coordinate units (m).
#' @param silent logical; should messages be printed?
#' @param sproj projected coordinate system of input UDs
#' @param dres desired pixel resolution of output UD
#' @param checksum logical; should adjustments be made to ensure sure that all cell volumes sum to 1?
#'
#' @return weighted mean spatRaster object
#' @export
#'
wmUD = function(udStack, w, sproj = NULL, silent = FALSE, dres = NULL, checksum = TRUE) {

  if(length(w) != terra::nlyr(udStack)) stop("w must be of same length as udList")

  if(!silent) {
    if(length(unique(w))==1) {
      cat("Calculating mean kernel...", "\n")
    } else {
      cat("Calculating weighted mean kernel...", "\n")
    }
  }

  # Ensure all UDs sum to 1 and convert to raster
  if(checksum) udList <- lapply(udList, function(i) reweightUD(i, out=1))

  # Convert to raster stack
  rstack <- rast(lapply(udList, function(x) UD2rast(x, sproj = sproj)))

  # ensure weights sum to sample size
  w = length(w) * w / sum(w)
  pstack = terra::weighted.mean(rstack[[w > 0]], w[w > 0])

  # resample to specified resolution
  if(!is.null(dres)) {
    if(unique(round(res(pstack))) != dres) {
      rastref <- rast(ext(pstack), crs = sproj, res = dres)
      if(!silent) message("Resampling to ", dres, "m...")
      pstack = terra::resample(x = pstack, y = rastref, method='bilinear')
    }
  }

  # ensure values sum to 1
  if(checksum) return(reweightUD(rast2UD(pstack, sproj = sproj), out=1)) else return(rast2UD(pstack, sproj = sproj))

}

