#' Batch deployment of multiple kernel density estimations based on levels of a specified factor.
#'
#' @param xy 2-column matrix or data.frame defining projected x and y coordinates of individual GPS relocations.
#' @param id vector defining unique sample animals from which individual GPS relocations were derived.
#' @param wts vector defining weights to be attributed to individual relocations. If the sum of relocation weights does not equal the sample size, weights are adjusted accordingly and a warning is issued.
#' @param ncores number of threads (processes) over which to dispatch individual kernel density estimators.
#' @param userGrid dimensions of the matrix over which the kernel will be estimated. Mandatory.
#' @param bwType specify bandwidth selection method. current options include plug-in ('pi' = ks::Hpi(); default argument), 'silv' (calls kernelboot::bw.silv()) and 'scott' (calls kernelboot::bw.scott()).
#' @param bwGlobal logical indicating whether bandwidth smoothing should be derived from all relocations (recommended) or made to vary according to individual sample (point pattern) distributions.
#' @param verbose logical indicating whether messages should be printed
#' @param write2file logical indicating whether results should be written to file
#' @param sproj optional EPSG
#' 
#' @importFrom stats setNames
#'
#' @return list of UD objects (list of x1, x2, fhat) of length equal to length(unique(id))
#' @export
#'
bKDE <- function(xy, id, wts = NULL, ncores = parallelly::availableCores() - 2,
                 userGrid = NULL, bwType = c('pi','silv','scott','user'), 
                 bwGlobal = TRUE, sproj = NULL, write2file = TRUE, verbose = TRUE) {
  
  bwType <- match.arg(bwType, choices = c('pi','silv','scott','user'), several.ok = F)
  bwSelect <- function(xyCoords, ...) {
    # if(bwType == 'user') {
    #   stop('not yet implemented')
    #   # return(bw)
    # }
    if(bwType == 'pi') return(ks::Hpi(xyCoords))
    if(bwType == 'silv') return(kernelboot::bw.silv(xyCoords)) else return(kernelboot::bw.scott(xyCoords))
  }
  
  if(nrow(xy) != length(id)) stop("id must be of length equal to nrow(xy)")
  if(is.null(wts)) {
    wts <- rep(1, nrow(xy))
  } else {
    if(length(wts) != nrow(xy) | !inherits(wts, 'numeric')) stop('wts must be a vector of numeric weights of length equal to nrow(xy)')
  }
  if(is.null(userGrid)) stop("userGrid must be provided")

  if(verbose) {
    if(!is.null(wts)) {
      message("Estimating ", length(unique(id)), " weighted utilization distributions...")
    } else {
      message("Estimating ", length(unique(id)), " unweighted utilization distributions...")
    }
  }
  
  if(ncores == 1) {
    
    kernelUDs <- setNames(terra::rast(
      sapply(1:length(unique(id)), function(m) {
      
      if(verbose) cat(paste(m, "/", length(unique(id)), sep=" "), "\n")
      tempxy <- xy[id %in% unique(id)[m], ]
      
      ## kernel density estimation using with or without spatial weights (scaled to sum to 1)
      wt <- wts[id %in% unique(id)[m]]
      wt <- length(wt) * wt / sum(wt)
      if(bwGlobal) H <- bwSelect(xy) else H <- bwSelect(tempxy)
      kmat <- kde(tempxy, w = wt, H = H,
                  xmin = c(userGrid$range.x[[1]][1], userGrid$range.x[[2]][1]),
                  xmax = c(userGrid$range.x[[1]][2], userGrid$range.x[[2]][2]),
                  gridsize = userGrid$grid.size,
                  density = TRUE)
      
      tf <- paste0(getwd(), '\\tmp\\kernel_', unique(id)[m], '.tif')
      terra::writeRaster(UD2rast(list(x1 = kmat$eval.points[[1]], 
                                      x2 = kmat$eval.points[[2]], 
                                      fhat = kmat$estimate), sproj = sproj), 
                         filename = tf, overwrite = T)
      return(tf)
      
    })), unique(id))

  } else {
    
    cl <- parallelly::makeClusterPSOCK(ncores, default_packages = c('sf','dplyr','terra','ks','wmKDE'))
    parallel::clusterExport(cl, varlist = c('xy','id','wts','bwGlobal','ncores','sproj',
                                'verbose','userGrid','write2file'), envir=environment())
    parallel::clusterEvalQ(cl, terra::terraOptions(memfrac = 0.5 / ncores,
                                                   memmax = 0.5))

    kernelUDs <- setNames(terra::rast(
      
      parallel::parSapplyLB(cl, 1:length(unique(id)), function(m) {
      
      tempxy <- xy[id %in% unique(id)[m], ]
      
      ## Kernel density estimation
      wt <- wts[id %in% unique(id)[m]]
      wt <- length(wt) * wt / sum(wt)
      if(bwGlobal) H <- ks::Hpi(xy) else H <- ks::Hpi(tempxy)
      
      kmat <- ks::kde(tempxy, w = wt, H = H,
                      xmin = c(userGrid$range.x[[1]][1], userGrid$range.x[[2]][1]),
                      xmax = c(userGrid$range.x[[1]][2], userGrid$range.x[[2]][2]),
                      gridsize = userGrid$grid.size,
                      density = T)
      
      tf <- paste0(getwd(), '\\tmp\\kernel_', unique(id)[m], '.tif')

      terra::writeRaster(
        x = wmKDE::UD2rast(
          list(x1 = kmat$eval.points[[1]],
               x2 = kmat$eval.points[[2]],
               fhat = kmat$estimate), 
          sproj = sproj),
        filename = tf, overwrite = T)
        
      rm(kmat)
      gc()
      
      return(tf)
      
    })), unique(id))
    
    parallel::stopCluster(cl)
    # rm(list=ls())
    # gc()
    
  }

  if(verbose) message("Analysis complete !")

  return(kernelUDs)

}
