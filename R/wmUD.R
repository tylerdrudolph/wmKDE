############################################################################################
## Produire un kernel pondéré à partir d'une liste de kernels et un vecteur de poids
#
# UDlist = liste de listes ayant des noms "x1", "x2", et "fhat"
# res = résolution désirée en mètres

#' Cell-wise weighted means of multiple Utilization Distributions
#'
#' @param UDlist list of utilization distributions of equal dimensions.
#' @param w vector of individual UD weights for cell-wise weighted means.
#' @param dres desired output spatial resolution of mean weighted UD in projected coordinate units (m)s.
#' @param silent logical; should messages be printed?
#' @param dcrs projected coordinate system of input UDs
#' @param checksum logical; should adjustments be made to ensure sure that all cell volumes sum to 1?
#'
#' @return
#' @export
#'
#' @examples
mwUD = function(UDlist, w, dres=NULL, silent=FALSE, dcrs=NULL, checksum=TRUE) {

  if(length(w) != length(UDlist)) stop("poids (w) doit être de la même longueur que UDlist")

  # Ensure all UDs sum to 1 and convert to raster
  if(checksum) UDlist = lapply(UDlist, function(i) reweightedUD(i, out=1))

  # Convert to raster stack
  UDprobz.stack = raster::stack(lapply(UDlist, UD2rast))

  # derive mean weighted UD
  if(!silent) {
    if(length(unique(w)==1)) cat("Calcul du kernel moyen", "\n") else cat("Calcul du kernel pondéré", "\n")
  }

  mwKUDprobz = raster::weighted.mean(UDprobz.stack, w)
    raster::crs(mwKUDprobz) <- dcrs

  # resample to specified resolution
  if(!is.null(dres)) {
    if(unique(round(raster::res(mwKUDprobz))) != dres) {
      rastref <- raster::raster(raster::extent(mwKUDprobz), crs=raster::crs(mwKUDprobz))
      raster::res(rastref) <- dres
      if(!silent) {
        raster::rasterOptions(overwrite=TRUE, progress="text", timer=TRUE)
        message("Échantillonnage au", dres, "m")
      }
      mwKUDprobz.resamp = raster::resample(x=mwKUDprobz, y=rastref)
      if(checksum) return(reweightedUD(rast2UD(mwKUDprobz.resamp), out=1)) else return(rast2UD(mwKUDprobz.resamp))
    }
  }

  # ensure values sum to 1
  if(checksum) return(reweightedUD(rast2UD(mwKUDprobz), out=1)) else return(rast2UD(mwKUDprobz))

}

