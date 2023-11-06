#' Weighted mean lattice density estimation for acoustic telemetry data
#'
#' @param x 'sf' class object whose coordinates correspond to individual acoustic telemetry locations. x must be in a projected coordinate system and must contain a field name assigned to argument 'id' when avg = TRUE.
#' @param id optional character vector indicating column name corresponding to unique tag ID from which each location was detected by a receptor. Must be provided when avg = TRUE.
#' @param bath 'sf' class POLYGON or MULTIPOLYGON object whose boundaries correspond to the water body under study (required).
#' @param avg logical indicating whether an averaged UD is desired. When avg = FALSE, only one kernel is generated using all detections.
#' @param spatres size of lattice mesh in metres (numeric of length 1)
#' @param udw optional 2-column data.frame containing individual UD weights. Must contain a field with name matching argument 'id' (one row per unique record) and a second field of weights entitled 'w'.
#' @param zscale logical indicating whether individual UD probability densities should be rescaled prior to cellwise averaging (recommended). Only applicable when avg = TRUE.
#' @param ktype character vector indicating type of kernel value to return. Options are 'iso' = isopleth contours (default, higher values indicate greater probability); 'prob' = UD probabilities; and 'vol' = UD volumes (sum to 1). Multiple arguments are accepted.
#' @param ncores integer indicating the number of parallel processes (threads) over which to execute the estimation of individual UDs. Defaults to parallelly::availableCores()-1 when avg = TRUE, otherwise 1.
#' @param write2file logical; write results to file?
#' @param ow logical; overwrite existing files?
#' @param fileTag optional string to append to file name when export=TRUE.
#' @param writeDir optional path to desired write folder location. Default is working directory.
#' @param retObj logical; should result be returned into R environment?
#' @param verbose logical; should status messages be printed?
#' 
#' @return  list containing 1) spatRaster object(s) of length equal to length(ktype), 2) one sf multipolygons object corresponding to isopleth, 3) core area contours and 4) 
#' @export
#' 
wmLDK <- function(x = NULL, id = NULL, bath = NULL, avg = TRUE, spatres = 10, udw = NULL, zscale = T,
                  ktype = c('iso','prob'), ncores = ifelse(avg, parallelly::availableCores() - 1, 1),
                  write2file = F, ow = F, fileTag = NULL, writeDir = getwd(), 
                  retObj = T, verbose = T) {
  
  df <- isopleth <- kernRasterList <- plevel <- NULL
  
  if(avg) {
    if(!is.null(id)) id <- match.arg(id, names(x), several.ok = F) else stop("Argument 'id' was not defined")
  }
  
  ## Create the estimation surface
  if(names(sort(table(st_geometry_type(bath)), decreasing=T)[1]) == 'LINESTRING') bath <- st_cast(bath, 'POLYGON')
  bath <- bath %>% summarize %>% st_make_valid
  bathcoord <- as.data.frame(st_coordinates(bath))
  nodeFillingOutput <- latticeDensity::nodeFilling(poly = sf::st_cast(bath, 'POLYGON') %>% 
                                                     st_make_valid %>% 
                                                     summarize %>% 
                                                     st_coordinates,
                                                   node_spacing = spatres)
  
  ## Create the node neighbourhood structure
  formLatticeOutput <- latticeDensity::formLattice(nodeFillingOutput)
  refbox <- sf::st_bbox(bath)
  fish.grid <- terra::rasterize(x = sf::st_as_sf(bath), 
                                y = terra::rast(xmin = refbox['xmin'], xmax = refbox['xmax'], 
                                                ymin = refbox['ymin'], ymax = refbox['ymax'],
                                                res = rep(spatres, 2), 
                                                crs = sf::st_crs(bath)$wkt),
                                field = 0)
  
  ## Create a temporary raster for each individual point pattern
  # dir.create(file.path(writeDir,'tmp'))
  if(avg) {
    idList <- dplyr::pull(x, id)
  } else {
    idList <- rep(1, nrow(x))
  }
  
  ncores <- min(c(length(unique(idList)), parallelly::availableCores() - 1))
  cl <- parallelly::makeClusterPSOCK(ncores, default_packages = c('terra','sf','dplyr','latticeDensity','Require'))
  parallel::clusterExport(cl, varlist = c('x','id','avg','idList','fish.grid','bath','formLatticeOutput'), envir = environment())
  
  system.time({
    
    latticeKernels <- parallel::parLapply(cl, 1:length(unique(idList)), function(i) {
      
      ## subset individual point pattern
      if(avg) {
        df.id <- dplyr::filter(x, !!as.name(id) == unique(idList[i]))
      } else {
        df.id <- df
      }
      
      ## retrieve # detections per grid cell
      id.grid <- terra::as.data.frame(
        terra::mask(
          terra::rasterize(x = df.id, y = fish.grid, 
                           fun = 'sum', background = 0), 
          terra::vect(bath)), 
        xy = TRUE) %>%
        dplyr::rename(n = sum)
      
      ## Determine the optimal # steps in the random walk
      xval <- latticeDensity::crossvalNparReg(formLatticeOutput,
                              Z = id.grid$n,
                              PointPattern = id.grid[, c('x','y')],
                              M = 0.5, max_steps = 10)
      
      ## Estimate continuous surface using lattice-based non-parametric regression (facultatif)
      tarp <- latticeDensity::createNparReg(formLatticeOutput, 
                                            Z = id.grid$n, 
                                            PointPattern = id.grid[, c("x","y")], 
                                            M = 0.5, k = xval$k)
      
      return(list(tarp = tarp, xval = xval))
      
    })
    
    names(latticeKernels) <- unique(idList)
    parallel::stopCluster(cl)
    
    sumtab <- cbind.data.frame(id = names(latticeKernels), 
                     t(sapply(latticeKernels, function(x) {
                       data.frame(k=x[[2]]$k, Sumsq=min(x[[2]]$SumSq)) 
                     }))) %>%
      tibble::remove_rownames()
    
    latticeKernels <- lapply(latticeKernels, function(x) return(x$tarp))
    
    ## Convert to 'pseudo UD' raster (continuous surface)
    resamp.grid <- terra::buffer(
      terra::rasterize(0, x = sf::st_as_sf(bath), 
                       y = terra::rast(xmin = refbox['xmin'], xmax = refbox['xmax'], 
                                       ymin = refbox['ymin'], ymax = refbox['ymax'],
                                       res = c(1, 1), 
                                       crs = sf::st_crs(bath)$wkt), values = 0), 
      width = 100)
    terra::NAflag(resamp.grid) <- FALSE
    
    ## Convert to spatRaster and resample for viewing quality 
    kernRastList <- terra::rast(lapply(1:length(latticeKernels), function(i) {
      
      x <- latticeKernels[[i]]
      xyres <- unique(c(x$EW_locs[2]-x$EW_locs[1], x$NS_locs[2]-x$NS_locs[1]))
      z <- x$NparRegMap
      
      iKern <- terra::mask(
        terra::resample(
          terra::cover(
            terra::rasterize(x = dplyr::bind_cols(z = z[as.numeric(row.names(x$nodes))], x$nodes) %>%
                               sf::st_as_sf(., coords = c('x','y'), crs = sf::st_crs(bath)$wkt), 
                             y = fish.grid, 
                             field = 'z'), 
            terra::classify(terra::buffer(fish.grid, 10), rcl = cbind(c(0,1), c(NA,0)))), 
          resamp.grid, method = 'cubicspline'), 
        terra::vect(bath))
      
      names(iKern) <- names(latticeKernels)[i]
      
      ## normalize values to sum to 1
      ival <- c(terra::values(iKern))
      ival[!is.na(ival)] <- ival[!is.na(ival)] / sum(ival, na.rm=T)
      terra::values(iKern) <- ival  
      
      return(iKern)
      
    }))
    
    ## Calculate (weighted) mean Utilization DistributionD
    if(is.null(udw)) {
      mKern <- terra::mean(kernRastList, na.rm = TRUE)
    } else {
      mKern <- terra::weighted.mean(x = kernRasterList, 
                                    w = udw$w[match(names(kernRasterList), udw[,id])], 
                                    na.rm = TRUE)
    }
    
    ## Estimate core area isopleth and generate isolines
    isocore <- core.area(rast2UD(mKern))
    isolines <- UD2sf(UD = rast2UD(mKern, sproj = sf::st_crs(df)$wkt), 
                      sproj = sf::st_crs(df), 
                      probs = c(isocore, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1)), 
                      spType = 'line') %>%
      dplyr::mutate(coreArea = ifelse(isopleth == isocore, TRUE, FALSE), .after = plevel) %>%
      dplyr::arrange(isopleth)
    
    mKern <- c(iso = mKern, prob = mKern, vol = mKern)
    terra::values(mKern$iso) <- 100 - fhat2confin(terra::values(mKern$iso))
    terra::values(mKern$prob) <- terra::values(mKern$vol) / terra::res(mKern$vol)[1] / terra::res(mKern$vol)[2]
    
    ## Visualize (optional)
    if(FALSE) {
      plot(mKern, breaks = 10)
      plot(isolines %>% st_geometry, col = ifelse(isolines$coreArea, 'red', 'darkgrey'), add = T)
    }
    
    return(list(mkern = mKern, isolines = isolines, xval = sumtab))
    
  })
  
}