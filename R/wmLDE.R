#' Weighted mean lattice density estimation for acoustic telemetry data
#'
#' @param x 'sf' class object whose coordinates correspond to individual acoustic telemetry locations. x must be in a projected coordinate system and must contain a field name assigned to argument 'id' when avg = TRUE.
#' @param id optional character vector indicating column name corresponding to unique tag ID from which each location was detected by a receptor. Must be provided when avg = TRUE.
#' @param bath 'sf' class POLYGON or MULTIPOLYGON object whose boundaries correspond to the water body under study (required).
#' @param spatres size of lattice mesh in metres (numeric of length 1)
#' @param udw optional 2-column data.frame containing individual UD weights. Must contain a field with name matching argument 'id' (one row per unique record) and a second field of weights entitled 'w'.
#' @param zscale logical indicating whether individual UD probability densities should be rescaled prior to cellwise averaging (recommended). Only applicable when avg = TRUE.
#' @param ktype character vector indicating type of kernel value to return. Options are 'iso' = isopleth contours (default, higher values indicate greater probability); 'prob' = UD probabilities; 'vol' = UD volumes (sum to 1). Multiple arguments are accepted.
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
wmLDK <- function(x, id = NULL, bath, spatres = 10, udw = NULL, zscale = T,
                  ktype = c('iso','prob','vol'), ncores = ifelse(avg, parallelly::availableCores() - 1, 1),
                  write2file = F, ow = F, fileTag = NULL, writeDir = getwd(), 
                  retObj = T, verbose = T) {
  
  Require::Require(c('dplyr','stringr','sf','latticeDensity','terra','wmKDE'))
  
  if(!is.null(id)) id <- match.arg(id, names(x), several.ok = F)
  
  ## Create the estimation surface
  nodeFillingOutput <- nodeFilling(poly = slot(slot(slot(as_Spatial(bath), "polygons")[[1]],"Polygons")[[1]],"coords"), 
                                   node_spacing = spatres)
  
  ## Create the node neighbourhood structure
  formLatticeOutput <- formLattice(nodeFillingOutput)
  refbox <- st_bbox(bath)
  fish.grid <- rasterize(x = st_as_sf(bath), 
                         y = rast(xmin = refbox['xmin'], xmax = refbox['xmax'], 
                                  ymin = refbox['ymin'], ymax = refbox['ymax'],
                                  res = rep(spatres, 2), 
                                  crs = st_crs(bath)$wkt),
                         field = 0)
  
  ## Create a temporary raster for each individual point pattern
  # dir.create(file.path(writeDir,'tmp'))
  idList <- unique(pull(x, id))
  # lapply(1:length(idList), function(i) {
  #   if(verbose) cat('i = ', i, '\n')
  #   saveRDS(as.data.frame(mask(rasterize(x = dplyr::filter(x, !!as.name(id) == idList[i]),
  #                                        y = fish.grid, fun = 'sum', background = 0), vect(bath)), xy = TRUE) %>%
  #             rename(n = sum), file = paste0(writeDir, '/tmp/tempRast_', i, '.rds'))
  # })
  
  cl <- parallelly::makeClusterPSOCK(ncores, default_packages = c('terra','sf','dplyr','latticeDensity','Require'))
  parallel::clusterExport(cl, varlist = c('x','id','idList','fish.grid','bath','formLatticeOutput'), envir = environment())
  
  system.time({
    
    latticeKernels <- parallel::parLapply(cl, 1:length(idList), function(i) {
      
      ## subset individual point pattern
      df.id <- dplyr::filter(df, !!as.name(id) == idList[i])
      
      ## retrieve points per grid cell
      # id.grid <- readRDS(paste0(writeDir, '/tmp/tempRast_', i, '.rds'))
      id.grid <- as.data.frame(mask(rasterize(x = dplyr::filter(x, !!as.name(id) == idList[i]),
                                              y = fish.grid, fun = 'sum', background = 0), vect(bath)), xy = TRUE) %>%
        rename(n = sum)
      
      ## Determine the optimal # steps in the random walk
      xval <- crossvalNparReg(formLatticeOutput,
                              Z = id.grid$n,
                              PointPattern = id.grid[, c('x','y')],
                              M = 0.5, max_steps = 10)
      
      ## Estimate continuous surface using lattice-based non-parametric regression (facultatif)
      tarp <- createNparReg(formLatticeOutput, Z = id.grid$n, PointPattern = id.grid[, c("x","y")], M = 0.5, k = xval$k)
      
      return(list(tarp = tarp, xval = xval))
      
    })
    
    names(latticeKernels) <- idList
    parallel::stopCluster(cl)
    
    sumtab <- cbind.data.frame(id = names(latticeKernels), 
                     t(sapply(latticeKernels, function(x) data.frame(k=x[[2]]$k, Sumsq=min(x[[2]]$SumSq))))) %>%
      tibble::remove_rownames()
    
    latticeKernels <- lapply(latticeKernels, function(x) return(x$tarp))
    
    ## Convert to 'pseudo UD' raster (continuous surface)
    resamp.grid <- buffer(rasterize(0, x = st_as_sf(bath), 
                                    y = rast(xmin = refbox['xmin'], xmax = refbox['xmax'], 
                                             ymin = refbox['ymin'], ymax = refbox['ymax'],
                                             res = c(1, 1), 
                                             crs = st_crs(bath)$wkt), values = 0), width=100)
    NAflag(resamp.grid) <- FALSE
    
    ##
    kernRastList <- rast(lapply(1:length(latticeKernels), function(i) {
      
      x <- latticeKernels[[i]]
      xyres <- unique(c(x$EW_locs[2]-x$EW_locs[1], x$NS_locs[2]-x$NS_locs[1]))
      z <- x$NparRegMap
      
      iKern <- mask(
        resample(
          cover(
            rasterize(x = bind_cols(z = z[as.numeric(row.names(x$nodes))], x$nodes) %>%
                        st_as_sf(., coords = c('x','y'), crs = st_crs(bath)$wkt), 
                      y = fish.grid, 
                      field = 'z'), 
            classify(buffer(fish.grid, 10), rcl = cbind(c(0,1), c(NA,0)))), 
          resamp.grid, method = 'cubicspline'), 
        vect(bath))
      
      names(iKern) <- names(latticeKernels)[i]
      
      ## normalize values to sum to 1
      ival <- c(values(iKern))
      ival[!is.na(ival)] <- ival[!is.na(ival)] / sum(ival, na.rm=T)
      values(iKern) <- ival  
      
      return(iKern)
      
    }))
    
    ## ! A FAIRE !
    # wKernz <- lapply(unique(rtab$week), function(p) {
    #   x <- lapply(c('SAFO','SAOQ'), function(spp) {
    #     mKern <- weighted.mean(x = kernRasterList[[rtab$week == p & rtab$Espece == spp]], w = rtab$poids[rtab$week == p & rtab$Espece == spp], na.rm = T)
    #     isocore <- core.area(wmKDE::rast2UD(mKern))
    #     isolines <- UD2sf(UD = rast2UD(mKern, sproj = st_crs(df)$wkt), sproj = st_crs(df), 
    #                       probs = c(isocore, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1)), 
    #                       spType = 'line') %>%
    #       mutate(coreArea = ifelse(isopleth == isocore, TRUE, FALSE), .after = plevel) %>%
    #       arrange(isopleth)
    #     values(mKern) <- 100 - fhat2confin(values(mKern))
    #     # plot(mKern, breaks = 10, main = paste(p, '-', spp))
    #     # plot(isolines %>% st_geometry, col = ifelse(isolines$coreArea, 'red', 'darkgrey'), add = T)
    #     return(list(mkern = mKern, isolines = isolines))
    #   })
    #   names(x) <- c('SAFO','SAOQ')
    #   return(x)
    # })
    # names(wKernz) <- unique(rtab$week)
    
    
  })
  
}