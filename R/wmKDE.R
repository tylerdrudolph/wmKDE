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
#' @param popGrid optional specification of grid over which to estimate the UD. Default is kernelGrid(locs = x, exp.range = 3, cell.size = spatres).
#' @param bwType specify bandwidth selection method. current options include plug-in ('pi' = ks::Hpi(); default argument), 'silv' (calls kernelboot::bw.silv()) and 'scott' (calls kernelboot::bw.scott()).
#' @param bwGlobal logical indicating whether bandwidth smoothing should be derived from all relocations (recommended) or made to vary according to individual sample (point pattern) distributions. Default uses plug-in method of bandwidth selection by default (alternative options not yet implemented).
#' @param zscale logical indicating whether individual UD probability densities should be rescaled prior to cellwise averaging (recommended). Only applicable when avg = TRUE.
#' @param spatres vector of length 1 specifying the desired spatial resolution of the output UD in the x & y dimensions. Asymmetrical cells not implemented.
#' @param ktype character vector indicating type of kernel value to return. Options are 'iso' = isopleth contours (default, higher values indicate greater probability); 'prob' = UD probabilities; 'vol' = UD volumes (sum to 1). Multiple arguments are accepted.
#' @param ncores integer indicating the number of parallel processes (threads) over which to execute the estimation of individual UDs. Defaults to detectCores()-1 when avg = TRUE, otherwise 1.
#' @param trim logical; should outer NA cells be excluded, where applicable, from the kernel raster?
#' @param write2file logical; write results to file?
#' @param ow logical; overwrite existing files?
#' @param fileTag optional string to append to file name when export=TRUE.
#' @param writeDir optional path to desired write folder location. Default is working directory.
#' @param retObj logical; should result be returned into R environment?
#' @param verbose logical; should status messages be printed?
#'
#' @return list containing spatRaster object(s) of length equal to length(ktype) and one sf multipolygons object corresponding to isopleth and core area contours
#' @export
#'
wmKDE <- function(x, id = NULL, avg = TRUE, spw = NULL, udw = NULL, popGrid = NULL,
                  bwType = c('pi', 'silv', 'scott'),
                  bwGlobal = TRUE, zscale = TRUE, spatres = 1000, ktype = 'iso',
                  ncores = ifelse(avg, parallel::detectCores() - 1, 1),
                  trim = TRUE, write2file = FALSE, ow = TRUE, 
                  writeDir = getwd(), fileTag = NULL,
                  retObj = TRUE, verbose = TRUE) {

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

  if(avg & !is.null(id)) {
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
  if(!is.null(spatres) & !inherits(spatres, 'numeric')) stop('spatres must be numeric')
  if(length(spatres) > 1) {
    message('length(spatres) > 1 but asymmetrical cells are not currently implemented. Only first element will be used.')
    spatres <- spatres[1]
  }
  
  bwType <- match.arg(bwType, choices = c('pi','silv','scott','user'), several.ok = F)

  ## Log system processing time
  ptime <- system.time({

    ## Generate the spatial grid over which the UD(s) are to be estimated
    if(is.null(popGrid)) {
      popGrid <- wmKDE::kernelGrid(x, exp.range = 3, cell.size = spatres)
    } else {
      spatres <- res(hgrid$r)[1]
    }

    ## Deploy kernel estimation
    if(avg) {

      ## Define the field serving to differentiate UDs
      ncores <- min(ncores, length(unique(idvec)))
      if(verbose) message("Estimating ", length(unique(idvec)), " Utilization Distributions (UD) across ", ncores, ' threads...')

    } else {

      if(verbose) message("Estimating a simple Utilization Distribution (UD)...")

    }

    ## Deploy multiple UD estimations
    udList <- wmKDE::bKDE(xy = xy, id = idvec, wts = wtvec, userGrid = popGrid, 
                          bwType = bwType, sproj = sproj, bwGlobal = bwGlobal, 
                          ncores = ifelse(avg, ncores, 1), verbose = FALSE)
    
    if(avg & terra::nlyr(udList) > 1) {

      ## Rescale z values
      if(zscale) {
        if(verbose) message("Rescaling density values...")
        udList <- sapp(udList, function(x) {
          x <- x * spatres ^ 2
          mm <- as.vector(terra::minmax(x))
          (x - mm[1]) / (mm[2] - mm[1])
        })
      }

      ## Derive the mean population UD (no weighting)
      if(!is.null(udw)) fileTag = 'weighted' else fileTag = 'unweighted'
      if(verbose) message(paste0("Deriving the ", fileTag, " mean population UD..."))
      if(is.null(udw)) w <- rep(1, length(udList)) else w <- udw$w[match(names(udList), pull(udw, id))]
      if(!zscale) udList <- terra::sapp(udList, wmKDE::finetune)
      wmKern <- terra::weighted.mean(udList, w = w)

      if(zscale) wmKern <- wmKern / unlist(global(wmKern, 'sum', na.rm = T)) / spatres / spatres

    } else {

      wmKern <- udList[[1]]

    }
    
    wmKern <- wmKDE::finetune(wmKern)

    ###########################################
    ## Calculate the core area isopleth
    if(verbose) message("Calculating the core area isopleth...")
    fileTag <- stringr::str_c(Sys.Date(), "_",
                   ifelse(avg, 'mean_', 'simple_'),
                   ifelse(!is.null(spw), 'weighted_kernel_', 'kernel_'),
                   spatres, "m")

    crit.core.isopleth <- wmKDE::core.area(rast2UD(wmKern))

    ############################################
    ## Extract isopleth polygons, including core/intensive use area
    isopoly <- wmKDE::UD2sf(UD = rast2UD(wmKern, sproj = sproj), sproj = sproj,
                            prob = sort(c(crit.core.isopleth, c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99, 1))))
    
    ## Determine the % of points falling within individual isopleth boundaries, including core area
    ftab <- terra::extract(wmKern, xy, fun = 'sum')
    names(ftab)[2] <- 'plevel'
    isopoly <- mutate(isopoly, 
                      pcntPnts = sapply(isopoly$plevel, function(iso) sum(ftab$plevel >= iso) / nrow(xy)),
                      coreArea = ifelse(isopleth == crit.core.isopleth, TRUE, FALSE), .before = geometry) %>%
      sf::st_cast('MULTIPOLYGON')

    wmKernRast <- wmKDE::rast2UD(wmKern, sproj = sproj)
    wmKernRast$fhat <- 100 - wmKDE::fhat2confin(wmKDE::rast2UD(wmKern, sproj = sproj)$fhat)
    wmKernRast <- wmKDE::UD2rast(wmKernRast, sproj)
    
    if(trim) {
      wmKernRast[wmKernRast < 0.05] <- NA
      wmKernRast <- terra::trim(wmKernRast)
    }

    ## Prepare output objects
    outlist <- list(wmKDE = terra::rast(list(iso = wmKernRast,
                                             prob = terra::crop(wmKern, wmKernRast),
                                             vol = terra::crop(wmKern * spatres * spatres, wmKernRast)))[[ktype]],
                    isocontours = isopoly)

    if(write2file) {

      if(verbose) message(paste0("Saving results to ", writeDir, "..."))

      ## Write READ ME file
      sink(file.path(writeDir, "wmKDE_parameters.txt"))
      cat('wmKDE: WEIGHTED MEAN KERNEL DENSITY ESTIMATOR', '\n', '\n')
      cat('Maintained by: Tyler Rudolph (tylerdrudolph@gmail.com)', '\n')
      cat(paste0("System date/time: ", Sys.time()), '\n', '\n')
      cat('MODEL PARAMETERS:', '\n', '\n')
      cat('METHODOLOGY:', '\n')
      cat('Type of analysis:', ifelse(!is.null(spw) | !is.null(udw), 'weighted', ''), ifelse(avg, 'mean', 'Simple'), 'kernel', '\n')
      cat('Spatial weights:', !is.null(spw), '\n')
      cat('UD weights:', !is.null(udw), '\n')
      cat("Number of parallel processes (threads):", ncores, "\n")
      if(!bwGlobal) cat("Bandwidth smoothing factor varies by individual", "\n") else cat("Global bandwidth smoothing factor", "\n")
      if(zscale) cat("Density (z) values rescaled", "\n") else cat("Density (z) values not rescaled", "\n")
      cat("Spatial resolution =", spatres, "m", "\n", "\n")
      cat('RELOCATION DATA:', '\n')
      cat("Number of unique collared individuals =", sum(!duplicated(idvec)), "\n")
      cat("Total number of relocations =", nrow(xy),  "\n")
      cat("Mean number of relocations =", round(mean(table(idvec))),  "\n")
      cat("Median number of relocations =", stats::median(table(idvec)),  "\n")
      cat("standard deviation =", round(stats::sd(table(idvec)), digits = 2), "\n")
      sink()

      ## Export final kernel(s)
      terra::writeRaster(outlist$wmKDE, filename = file.path(writeDir, paste0(fileTag, '.tif')), overwrite = ow)

      ## Export isopleth contours
      sf::st_write(isopoly, file.path(writeDir, paste0(fileTag, '_isopleth_contours.gpkg')), delete_layer = ow, quiet = TRUE)

    }

  })

  if(verbose) {
    message('Procedure complete!')
    message('{ computation time = ', round(unname(ptime[3]) / 60, digits = 2), ' mins }')
  }

  if(retObj) return(outlist) else return(invisible(NULL))

}
