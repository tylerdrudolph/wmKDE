#' Batch deployment of multiple kernel density estimations based on levels of a specified factor.
#'
#' @param xy 2-column matrix or data.frame defining projected x and y coordinates of individual GPS relocations.
#' @param id vector defining unique sample animals from which individual GPS relocations were derived.
#' @param wts vector defining weights to be attributed to individual relocations. If the sum of relocation weights does not equal the sample size, weights are adjusted accordingly and a warning is issued.
#' @param ncores number of threads (processes) over which to dispatch individual kernel density estimators.
#' @param user.grid dimensions of the matrix over which the kernel will be estimated. Mandatory.
#' @param bw.global logical indicating whether bandwidth smoothing should be derived from all relocations (recommended) or made to vary according to individual sample (point pattern) distributions.
#' @param verbose logical indicating whether messages should be printed
#' @param write2file logical indicating whether results should be written to file
#'
#' @return list of UD objects (list of x1, x2, fhat) of length equal to length(unique(id))
#' @export
#'
#' @examples
#'
bKDE <- function(xy, id, wts = NULL, ncores = parallel::detectCores() - 1,
                 user.grid = NULL, bw.global = TRUE, write2file = FALSE,
                 verbose = TRUE) {

  if(nrow(xy) != length(id)) stop("id must be of length equal to nrow(xy)")
  if(is.null(wts)) {
    wts <- rep(1, nrow(xy))
  } else {
    if(length(wts) != nrow(xy) | !inherits(wts, 'numeric')) stop('wts must be a vector of numeric weights of length equal to nrow(xy)')
  }
  if(is.null(user.grid)) stop("user.grid must be provided")

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
      if(bw.global) H <- Hpi(xy) else H <- Hpi(tempxy)
      kmat <- kde(tempxy, w=wt, H=H,
                      xmin = c(user.grid$range.x[[1]][1], user.grid$range.x[[2]][1]),
                      xmax = c(user.grid$range.x[[1]][2], user.grid$range.x[[2]][2]),
                      gridsize = user.grid$grid.size[1])

      kernelUDs = c(kernelUDs, list(list(x1 = kmat$eval.points[[1]], x2 = kmat$eval.points[[2]], fhat = kmat$estimate)))

    }

  } else {

    cl <- parallelly::makeClusterPSOCK(ncores, default_packages = c('sf','dplyr','terra','ks'))
    parallel::clusterExport(cl, varlist=c('levz','xy','id','wts','bw.global',
                                'verbose','user.grid'), envir=environment())

    kernelUDs <- parallel::parLapply(cl, 1:length(levz), function(m) {

      tempxy <- xy[id %in% levz[m], ]

      ## Kernel density estimation
      wt <- wts[id %in% levz[m]]
      wt <- length(wt) * wt / sum(wt)
      if(bw.global) H <- ks::Hpi(xy) else H <- ks::Hpi(tempxy)
      kmat <- ks::kde(tempxy, w=wt, H=H,
                  xmin = c(user.grid$range.x[[1]][1], user.grid$range.x[[2]][1]),
                  xmax = c(user.grid$range.x[[1]][2], user.grid$range.x[[2]][2]),
                  gridsize = user.grid$grid.size[1])

      return(list(x1=kmat$eval.points[[1]], x2=kmat$eval.points[[2]], fhat=kmat$estimate))

    })

    parallel::stopCluster(cl)

  }

  if(verbose) message("Analysis complete !")

  names(kernelUDs) = levz

  return(kernelUDs)

}

