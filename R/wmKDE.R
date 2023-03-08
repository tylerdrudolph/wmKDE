#' Weighted mean kernel density estimation
#'
#' Estimate the mean utilization distribution (UD) of a wildlife population given GPS collar telemetry relocations derived from individually monitored sample animals.
#' Two levels of weighting are possible:
#' 1) Individual GPS relocations may be weighted during kernel density estimation to account for variation in data quality and/or habitat transformation over time.
#' 2) Individual UDs may be weighted during cell-wise averaging to account for individual variation in sample size and/or quality.
#'
#' @param x 'sf' class object whose coordinates correspond to individual GPS relocations. x must be in a projected coordinate system and must contain a field name assigned to argument 'id' when avg=TRUE.
#' @param id optional character vector indicating column name corresponding to unique collar ID from which each relocated was transmitted. Must be provided when avg = TRUE.
#' @param avg logical indicating whether an averaged UD is desired. When avg=FALSE, only one kernel is generated using all relocations.
#' @param spw optional numeric vector of spatial relocation weights, length of which must be equal to nrow(x). Weights will be rescaled if sum(spw) != nrow(x).
#' @param udw optional 2-column data.frame containing individual UD weights. Must contain a field with name matching argument 'id' (one row per unique record) and a second field of weights entitled 'w'.
#' @param herd.grid optional specification of grid over which to estimate the UD. Default is kernel.grid(locs = st_coordinates(x), exp.range = 3, cell.size = spatres).
#' @param bw.global logical indicating whether bandwidth smoothing should be derived from all relocations (recommended) or made to vary according to individual sample (point pattern) distributions. Default uses plug-in method of bandwidth selection by default (alternative options not yet implemented).
#' @param zscale logical indicating whether individual UD probability densities should be rescaled prior to cellwise averaging (recommended). Only applicable when avg = TRUE.
#' @param spatres vector of length 1 specifying the desired spatial resolution of the output UD in the x & y dimensions. Asymmetrical cells not implemented.
#' @param ktype character vector indicating type of kernel value to return. Options are 'iso' = isopleth contours (default, higher values indicate greater probability); 'prob' = UD probabilities; 'vol' = UD volumes (sum to 1). Multiple arguments are accepted.
#' @param ncores integer indicating the number of parallel processes (threads) over which to execute the estimation of individual UDs. Defaults to detectCores()-1 when avg = TRUE, otherwise 1.
#' @param write2file logical; write results to file?
#' @param ow logical; overwrite existing files?
#' @param fileTag optional string to append to file name when export=TRUE.
#' @param writeDir optional path to desired write folder location. Default is working directory.
#' @param retObj logical; should result be returned into R environment?
#'
#' @return list containing spatRaster object(s) of length equal to length(ktype) and one sf multipolygons object corresponding to isopleth and core area contours
#' @export
#'
#' @examples
#'
wmKDE <- function(x, id = NULL, avg = TRUE, spw = NULL, udw = NULL, herd.grid = NULL,
                  bw.global = TRUE, zscale = TRUE, spatres = 1000, ktype = 'iso',
                  ncores = ifelse(avg, parallel::detectCores() - 1, 1),
                  write2file = FALSE, ow = TRUE, writeDir = getwd(), fileTag = NULL,
                  retObj = TRUE) {

  X <- Y <- w <- layer <- isopleth <- geometry <- plevel <- NULL

  ## Convert from spatial where applicable
  if(inherits(x, 'Spatial')) x <- sf::st_as_sf(x)

  ## validate input arguments
  if(inherits(x, 'sf')) {
    if(sf::st_is_longlat(x)) stop('x is not in a projected coordinate reference system')
    sproj <- sf::st_crs(x)
    xy <- as.data.frame(sf::st_coordinates(x)) %>%
      rename(x = X, y = Y)
  } else {
    stop('x must be a simple features object')
  }

  if(!is.null(id)) {
    id <- match.arg(id, names(x), several.ok = F)
    idvec <- pull(x, id) %>% as.character()
  } else {
    idvec <- rep(1, nrow(x))
  }

  if(!is.null(spw)) {
    spw <- match.arg(spw, names(x), several.ok = F)
    wtvec <- pull(x, spw)
    if(!inherits(wtvec, 'numeric'))  stop('spw not numeric')
  } else {
    wtvec <- rep(1, nrow(x))
  }

  if(avg & !is.null(udw)) {
    if(is.null(id)) stop("cannot apply UD weights when 'id' is not specified")
      if(!id %in% names(udw)) stop(paste0("No '", id, "' field in udw"))
    if(any(!idvec %in% (unique(pull(udw, id) %>% as.character)))) stop(paste0(id, ' values missing from udw'))
  }

  ktype <- match.arg(ktype, c('iso', 'prob', 'vol'), several.ok = T)
  if(!inherits(spatres, 'numeric')) stop('spatres must be numeric')
  if(length(spatres) > 1) {
    message('length(spatres) > 1 but asymmetrical cells are not currently implemented. Only first argument will be used.')
    spatres <- spatres[1]
  }

  ## Log system processing time
  ptime <- system.time({

    ## Generate the spatial grid over which the UD(s) are to be estimated
    if(is.null(herd.grid)) herd.grid <- kernel.grid(xy, exp.range = 3, cell.size = spatres)

    ## Deploy kernel estimation
    if(avg) {

      ## Define the field serving to differentiate UDs
      ncores <- min(ncores, length(unique(idvec)))
      message("Estimating ", length(unique(idvec)), " Utilization Distributions (UD) across ", ncores, ' threads...')

    } else {

      message("Estimating a simple Utilization Distribution (UD)...")

    }

    ## Deploy multiple UD estimations
    udList <- bKDE(xy = xy, id = idvec, wts = wtvec, user.grid = herd.grid,
                   bw.global = bw.global, ncores = ncores, verbose = FALSE)

    if(avg & length(udList) > 1) {

      ## Rescale z values
      if(zscale) {

        message("Rescaling density values...")
        udList <- lapply(udList, function(x) {
          x$fhat <- range01(x$fhat * spatres * spatres)
          return(x)
        })

      }

      ## Derive the mean population UD (no weighting)
      if(!is.null(udw)) fileTag = 'weighted' else fileTag = 'unweighted'
      message(paste0("Deriving the ", fileTag, " mean population UD..."))
      if(spatres != 1000) message("Resampling to ", spatres, "m...")
      if(is.null(udw)) w <- rep(1, length(udList)) else w <- udw$w[match(names(udList), pull(udw, id))]
      wmKern <- wmUD(udList, w = w, sproj = sproj, checksum =! zscale, silent = TRUE)

      if(zscale) wmKern$fhat <- wmKern$fhat / sum(wmKern$fhat) / spatres / spatres

    } else {

      wmKern <- udList[[1]]

    }

    ###########################################
    ## Calculate the core area isopleth
    message("Calculating the core area isopleth...")
    fileTag <- stringr::str_c(Sys.Date(), "_",
                   ifelse(avg, 'mean_', 'simple_'),
                   ifelse(!is.null(spw), 'weighted_kernel_', 'kernel_'),
                   # if(!is.null(fileTag)) stringr::str_c(fileTag, '_'),
                   spatres, "m")

    crit.core.isopleth <- core.area(wmKern)

    ############################################
    ## Extract isopleth polygons, including core/intensive use area
    isopoly <- UD2sf(UD = wmKern, sproj = sproj,
                       probs = sort(c(crit.core.isopleth, c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99, 1))))

    ## Determine the % of points falling within individual isopleth boundaries, including core area
    ftab <- terra::extract(UD2rast(wmKern, sproj), xy)
    names(ftab)[2] <- 'plevel'
    isopoly <- mutate(isopoly, pcntPnts = sapply(isopoly$plevel, function(iso) sum(ftab$plevel >= iso) / nrow(xy)),
                      coreArea = ifelse(isopleth == crit.core.isopleth, TRUE, FALSE), .before = geometry)

    wmKernRast <- wmKern
    wmKernRast$fhat <- 100 - fhat2confin(wmKern$fhat)
    wmKernRast <- UD2rast(wmKernRast, sproj)
    wmKernRast[wmKernRast < 0.05] <- NA
    wmKernRast <- terra::trim(wmKernRast)

    if(write2file) {

      message(paste0("Saving results to ", writeDir, "..."))

      ## Produce READ ME file
      sink(file.path(writeDir, "wmKDE_parameters.txt"))
      cat('wmKDE: WEIGHTED MEAN KERNEL DENSITY ESTIMATION', '\n', '\n')
      cat('Maintained by: Tyler Rudolph (tylerdrudolph@gmail.com)', '\n')
      cat(paste0("System date/time: ", Sys.time()), '\n', '\n')
      cat('MODEL PARAMETERS:', '\n', '\n')
      cat('METHODOLOGY:', '\n')
      cat('Type of analysis:', ifelse(!is.null(spw) | !is.null(udw), 'weighted', ''), ifelse(avg, 'mean', 'Simple'), 'kernel', '\n')
      cat('Spatial weights:', !is.null(spw), '\n')
      cat('UD weights:', !is.null(udw), '\n')
      cat("Number of parallel processes (threads):", ncores, "\n")
      if(!bw.global) cat("Bandwidth smoothing factor varies by individual", "\n") else cat("Global bandwidth smoothing factor", "\n")
      if(zscale) cat("Density (z) values rescaled", "\n") else cat("Density (z) values not rescaled", "\n")
      cat("Spatial resolution =", spatres, "m", "\n", "\n")
      cat('RELOCATION DATA:', '\n')
      cat("Number of unique collared individuals =", sum(!duplicated(idvec)), "\n")
      cat("Total number of relocations =", nrow(xy),  "\n")
      cat("Mean number of relocations =", round(mean(table(idvec))),  "\n")
      cat("Median number of relocations =", stats::median(table(idvec)),  "\n")
      cat("standard deviation =", round(stats::sd(table(idvec)), digits = 2), "\n")
      sink()

    }

    ## Export raster kernel(s)
    if(retObj) {

      outlist <- c()

      if('iso' %in% ktype) {
        terra::writeRaster(wmKernRast, filename = file.path(writeDir, paste0(fileTag, '_iso.tif')), overwrite = ow)
        outlist <- c(outlist, iso = wmKernRast)
      }
      if('prob' %in% ktype) {
        terra::writeRaster(UD2rast(wmKern, sproj), filename = file.path(writeDir, paste0(fileTag, '_prob.tif')), overwrite = ow)
        outlist <- c(outlist, prob = terra::crop(UD2rast(wmKern, sproj), wmKernRast))
      }
      if('vol' %in% ktype) {
        terra::writeRaster(UD2rast(wmKern, sproj) * spatres * spatres, filename = file.path(writeDir, paste0(fileTag, '_vol.tif')), overwrite = ow)
        outlist <- c(outlist, vol = terra::crop(UD2rast(wmKern, sproj) * spatres * spatres, wmKernRast))
      }

      ## Export isopleth contours
      sf::st_write(isopoly, file.path(writeDir, paste0(fileTag, '_isopleth_contours.gpkg')), delete_layer = ow, quiet = TRUE)
      outlist <- list(wmKDE = terra::rast(outlist), isocontours = isopoly)

    }

  })

  if(retObj) return(outlist)

}
