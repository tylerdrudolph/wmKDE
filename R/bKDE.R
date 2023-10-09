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
#'
#' @return list of UD objects (list of x1, x2, fhat) of length equal to length(unique(id))
#' @export
#'
bKDE <- function(xy, id, wts = NULL, ncores = parallel::detectCores() - 1,
                 userGrid = NULL, bwType = c('pi','silv','scott'), bwGlobal = TRUE, 
                 write2file = FALSE, verbose = TRUE) {
  
  bwType <- match.arg(bwType, choices = c('pi','silv','scott'), several.ok = F)

  if(nrow(xy) != length(id)) stop("id must be of length equal to nrow(xy)")
  if(is.null(wts)) {
    wts <- rep(1, nrow(xy))
  } else {
    if(length(wts) != nrow(xy) | !inherits(wts, 'numeric')) stop('wts must be a vector of numeric weights of length equal to nrow(xy)')
  }
  if(is.null(userGrid)) stop("userGrid must be provided")

  kernelUDs = list()
  levz <- unique(id)

  if(verbose) {
    if(!is.null(wts)) {
      message("Estimating ", length(unique(id)), " weighted utilization distributions...")
    } else {
      message("Estimating ", length(unique(id)), " unweighted utilization distributions...")
    }
  }
  
  if(ncores == 1) {

    for(m in 1:length(levz)) {

      if(verbose) cat(paste(m, "/", length(levz), sep=" "), "\n")
      tempxy <- xy[id %in% levz[m], ]

      ## kernel density estimation using with or without spatial weights (scaled to sum to 1)
      wt <- wts[id %in% levz[m]]
      wt <- length(wt) * wt / sum(wt)
      bwSelect <- function(xyCoords, ...) {
        if(bwType == 'pi') return(ks::Hpi(xyCoords))
        if(bwType == 'silv') return(kernelboot::bw.silv(xyCoords)) else return(kernelboot::bw.scott(xyCoords))
      }
      if(bwGlobal) H <- bwSelect(xy) else H <- bwSelect(tempxy)
      kmat <- kde(tempxy, w=wt, H=H,
                      xmin = c(userGrid$range.x[[1]][1], userGrid$range.x[[2]][1]),
                      xmax = c(userGrid$range.x[[1]][2], userGrid$range.x[[2]][2]),
                      gridsize = userGrid$grid.size,
                  density = TRUE)

      kernelUDs = c(kernelUDs, list(list(x1 = kmat$eval.points[[1]], x2 = kmat$eval.points[[2]], fhat = kmat$estimate)))

    }

  } else {

    cl <- parallelly::makeClusterPSOCK(ncores, default_packages = c('sf','dplyr','terra','ks'))
    parallel::clusterExport(cl, varlist=c('levz','xy','id','wts','bwGlobal',
                                'verbose','userGrid'), envir=environment())

    kernelUDs <- parallel::parLapply(cl, 1:length(levz), function(m) {

      tempxy <- xy[id %in% levz[m], ]

      ## Kernel density estimation
      wt <- wts[id %in% levz[m]]
      wt <- length(wt) * wt / sum(wt)
      if(bwGlobal) H <- ks::Hpi(xy) else H <- ks::Hpi(tempxy)
      kmat <- ks::kde(tempxy, w=wt, H=H,
                  xmin = c(userGrid$range.x[[1]][1], userGrid$range.x[[2]][1]),
                  xmax = c(userGrid$range.x[[1]][2], userGrid$range.x[[2]][2]),
                  gridsize = userGrid$grid.size,
                  density = T)

      return(list(x1=kmat$eval.points[[1]], x2=kmat$eval.points[[2]], fhat=kmat$estimate))

    })

    parallel::stopCluster(cl)

  }

  if(verbose) message("Analysis complete !")

  names(kernelUDs) = levz

  return(kernelUDs)

}

