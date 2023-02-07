#' Mean weighted kernel density estimation
#'
#' Estimate the mean utilization distribution (UD) of a wildlife population given GPS collar telemetry relocations derived from individually monitored sample animals.
#' Two levels of weighting are possible:
#' 1) Individual GPS relocations may be weighted during kernel density estimation to account for variation in data quality and/or habitat transformation over time.
#' 2) Individual UDs may be weighted during cell-wise averaging to account for individual variation in sample size and/or quality.
#'
#' @param xy a Spatial* or 'sf' class object, the first 2 columns of which correspond to projected GPS relocation coordinates (x, y).
#' @param avg logical indicating whether an averaged UD is desired. When avg=FALSE, only one kernel is generated using all relocations.
#' @param id vector of length equal to nrow(xy) specifying the unique collared animal from which each relocated was transmitted. Only necessary when avg=TRUE.
#' @param spwt logical indicating whether weights should be applied to individual relocations.
#' @param wts vector equal in length to nrow(xy) that defines the spatial weights applied to individual relocations. If the sum of relocation weights does not equal the sample size, weights are adjusted accordingly and a warning is issued. Only applicable when spwt=TRUE.
#' @param bw.global logical indicating whether bandwidth smoothing should be derived from all relocations (recommended) or made to vary according to individual sample (point pattern) distributions.
#' @param zscale logical indicating whether individual UD density values should be rescaled prior to cellwise averaging (recommended)0. Only applicable when avg=TRUE.
#' @param spatres vector specifying the desired spatial resolution of the output UD in the x & y dimensions. Arguments of length 1 are recycled.
#' @param ncores integer indicating the number of parallel processes (threads) over which to execute the estimation of individual UDs. Defaults to detectCores()-1 when avg=TRUE, otherwise 1.

#' @param ow logical; overwrite existing files?
#' @param titre optional title to print to plot
#' @param export logical; write results to file?
#' @param obj.ret logical; should results be returned?
#' @param herd.grid optional specification of grid over which to estimate the UD. Generated with kernel.grid().
#' @param fappend optional string to append to file name when export=TRUE.
#'
#' @return
#' @export
#'
#' @examples
mwKDE <- function(xy, avg=TRUE, id, spwt=TRUE, wts=NULL, bw.global=TRUE,
                  zscale=TRUE, spatres=1000, ncores=ifelse(avg, parallel::detectCores() - 1, 1),
                  ow=TRUE, titre=NULL, export=TRUE, obj.ret=FALSE, herd.grid=NULL,
                  pwm=NULL, fappend=NULL) {

  ## Convert from spatial where applicable
  if(inherits(xy, 'Spatial')) {
    sproj <- sp::proj4string(xy)
    xy <- xy@data[,c("x","y")]
  }
  if(inherits(xy, 'sf')) {
    #sproj <- st_crs(xy)
    sproj = sp::CRS("+init=epsg:32198")
    xy <- st_coordinates(xy)
  }

  if(inherits(xy, 'matrix')) xy <- as.data.frame(xy)
  if(is.factor(id)) id <- droplevels(id)

  ## Validate inputs
  if(ncol(xy)!=2) {
    stop('xy must contain two columns')
  } else {
    if(!all(apply(xy, 2, is.numeric))) stop('xy columns must be numeric')
    names(xy)[1:2] <- c('x','y')
  }

  ## Log system processing time
  ptime <- system.time({

    raster::rasterOptions(progress='text')
    options(warn=-1)

    ## Generate the spatial grid over which the UD(s) are to be estimated
    if(is.null(herd.grid)) herd.grid <- kernel.grid(xy, exp.range=3, cell.size=spatres)

    ## Define spatial weights if applicable
    if(spwt & length(wts) != nrow(xy)) stop("Must provide a vector 'wts' corresponding to relocation weights")

    ## Deploy kernel estimation
    if(avg) {

      ## Define the field serving to differentiate UDs
      message("Estimating ", length(unique(id)), " Utilization Distributions (UD) across ", ncores, ' threads...')

    } else {

      id <- rep(1, nrow(xy))
      message("Estimating a simple Utilization Distribution (UD)...")

    }

    ## Deploy multiple UD estimations
    UD_list <- bKDE(xy=xy, id=id, wts=wts, user.grid=herd.grid, bw.global=bw.global, ncores=ncores, verbose=FALSE)

    if(avg & length(UD_list) > 1) {

      ## Rescale z values
      if(zscale) {

        message("Rescaling density values...")
        UD_list <- lapply(UD_list, function(x) {
          x$fhat <- range01(x$fhat * spatres * spatres)
          return(x)
        })

        ## Derive the mean population UD (no weighting)
        message("Deriving the unweighted mean population UD...")
        if(spatres != 1000) message("Échantillonnage au ", spatres, "m...")
        if(!is.null(pwm)) pwm <- pwm else pwm = rep(1, length(UD_list))
        mwKern <- mwUD(UD_list, w=pwm, dcrs=sproj, checksum=!zscale, silent=TRUE)
        mwKern$fhat <- mwKern$fhat / sum(mwKern$fhat) / spatres / spatres

      } else {

        ## NOT YET IMPLEMENTED - MUST DEFINE 'w'.
        message("Deriving the weighted mean population UD...")
        if(spatres != 1000) message("Échantillonnage au ", spatres, "m...")
        mwKern <- mwUD(UD_list, w=rep(1, length(UD_list)), dcrs=sproj, dres=spatres, checksum=!zscale, silent=TRUE)

      }

    } else {

      mwKern <- UD_list[[1]]

    }

    ###########################################
    ## Calculate the core area isopleth
    message("Calculating the core area isopleth...")
    fname <- stringr::str_c(Sys.Date(), "_",
                   ifelse(avg, 'mean_', 'simple_'),
                   ifelse(spwt, 'weighted_kernel_', 'kernel_'),
                   if(!is.null(fappend)) str_c(fappend, '_'),
                   spatres, "m")

    crit.core.isopleth <- data.frame(fname, critval = core.area(mwKern))

    ############################################
    ## Extraire polygones correspondant aux ZUI
    core.poly = UD2sp(UD=mwKern, PID=fname, probs=crit.core.isopleth$critval, proj=sproj)

    ## Déterminer le % des localisations se trouvant au sein de la ZUI
    cua.percent = percent.pts.in.poly(xy=xy, poly=core.poly)
    tdata = cbind.data.frame(fname, confin=crit.core.isopleth$critval, cua.percent)
    row.names(tdata) <- fname
    core.area.spdf <- sp::SpatialPolygonsDataFrame(core.poly, tdata)
    methods::slot(core.area.spdf, "polygons") <- lapply(methods::slot(core.area.spdf, "polygons"), maptools::checkPolygonsHoles)

    ##################################
    ## Calculer les contours de probabilités
    message('Extracting the isopleths...')
    probz <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99, 1)
    library(maptools) ## write this properly into UD2sp !!
    ci_contours <- do.call(rbind,
                           lapply(probz, FUN=function(p) {
                             if(p < 1) {
                               x <- UD2sp(mwKern, PID=p, probs=p, proj=sproj)
                             } else {
                               x <- UD2sp(mwKern, PID=p, levels=1e-13, proj=sproj)
                             }
                             xdata <- data.frame(fname, prob=p)
                             row.names(xdata) <- p
                             return(sp::SpatialPolygonsDataFrame(x, xdata))
                           }))

    ############################
    ## Sauvegarder le résultat
    if(export) {

      ## define empty temp folder
      if(file.exists(tempdir())) unlink(str_c(tempdir(), '/', dir(tempdir())), recursive=T)
      dir.create(tempdir())

      message("Saving results...")
      ## Produire lisez-moi.txt file
      sink(stringr::str_c(tempdir(), "/mwKDE_parametres.txt"))
      cat('mwKDE: MEAN WEIGHTED KERNEL DENSITY ESTIMATION', '\n', '\n')
      cat('Developed by Tyler Rudolph, M.Sc. (tylerdrudolph@gmail.com)', '\n')
      cat(stringr::str_c("System date/time: ", Sys.time()), '\n', '\n')
      cat('MODEL PARAMETERS:', '\n', '\n')
      cat('METHODOLOGY:', '\n')
      cat("Type of analysis:", ifelse(avg, 'Mean', 'Simple'), ifelse(spwt, 'weighted kernel', 'kernel'), '\n')
      cat("Number of parallel processes (threads):", ncores, "\n")
      if(!bw.global) cat("Bandwidth smoothing factor varies by individual", "\n") else cat("Global bandwidth smoothing factor", "\n")
      if(zscale) cat("Density (z) values rescaled", "\n") else cat("Density (z) values not rescaled", "\n")
      cat("Spatial resolution =", spatres, "m", "\n", "\n")
      cat('RELOCATION DATA:', '\n')
      cat("Number of unique collared individuals =", sum(!duplicated(id)), "\n")
      cat("Total number of relocations =", nrow(xy),  "\n")
      cat("Mean number of relocations =", mean(table(id)),  "\n")
      cat("Median number of relocations =", median(table(id)),  "\n")
      cat("standard deviation =", sd(table(id)), "\n")
      sink()

      ## Sauvegarder le kernel
      save(mwKern, file=stringr::str_c(tempdir(), '/', fname, ".RData"))

    }

    #######################################################
    ## Valider l'ajustement des kernels individuels à la ZUI globale
    validtab <- data.frame(table(id))
    names(validtab) <- c('IDAnimal','n.total')
    validtab$pct.in.zui <- sapply(as.character(validtab$IDAnimal), function(i) {
     percent.pts.in.poly(xy=xy[id==i,], poly=core.poly)
    })

     if(spwt) {

       validtab <- cbind(validtab, wpct.in.zui = sapply(as.character(validtab$IDAnimal), function(i) {
         percent.pts.in.poly(xy=xy[id==i,], w=wts[id==i], poly=core.poly)
       }))

       validtab$net.diff <- validtab$wpct.in.zui - validtab$pct.in.zui

      validtab$IDAnimal <- as.character(validtab$IDAnimal)
      #validtab <- rbind(c('TOUS', sum(validtab$n.total), percent.pts.in.poly(xy=locs, poly=core.poly),
      #                  percent.pts.in.poly(xy=locs, w=wts, poly=core.poly), NA), validtab)

      vstats <- cbind.data.frame(c("Mean", "Std.Dev"),
                                 rbind(mean=round(colMeans(validtab[,c(3:5)]), digits=2),
                                       sd=round(apply(validtab[,c(3:5)], 2, stats::sd), digits=2)))

      names(vstats) <- c("Parameter", '%_in_coreArea', '%_in_coreArea_wtd', 'net.diff')
      names(validtab)[3:4] <- c('%_in_coreArea', '%_in_coreArea_wtd')

    } else {

      vstats <- cbind.data.frame(c("Mean", "Std.Dev"),
                                 rbind(mean=round(mean(validtab$pct.in.zui), digits=2),
                                       sd=round(stats::sd(validtab$pct.in.zui), digits=2)))
      names(vstats) <- c("Parameter", '%_in_coreArea')
      names(validtab)[3] <- '%_in_coreArea'

    }

    if(export) saveRDS(list(vtab=validtab[order(validtab$`%_in_coreArea`, decreasing=T),], stats=vstats), file=str_c(tempdir(), '/validtab.rds'))

    #################################
    ## Exporter les contours de probabilités en format vectoriel
    #methods::slot(ci_contours, "polygons") <- lapply(methods::slot(ci_contours, "polygons"), maptools::checkPolygonsHoles)
    if(export) st_write(sf::st_buffer(sf::st_as_sf(ci_contours), 0.0), dsn=tempdir(), layer=fname, driver="ESRI Shapefile", delete_layer=ow)
    #if(export) rgdal::writeOGR(ci_contours, dsn=tempdir(), layer=fname, driver="ESRI Shapefile", overwrite_layer=ow)

    #######################
    ## Exporter la ZUI en format fichier de forme
    if(export) rgdal::writeOGR(core.area.spdf, dsn=tempdir(), layer=stringr::str_c("ZUI_", fname), driver="ESRI Shapefile", overwrite=ow)

    ###################################
    ## Exporter le kernel pondéré en format raster (.tifs), valeurs de 0 à 100
    if(obj.ret) outKern <- mwKern
    mwKern$fhat <- fhat2confin(mwKern$fhat)
    if(export) raster::writeRaster(UD2rast(mwKern), filename=stringr::str_c(tempdir(), "/", fname, ".tif"), format="GTiff", overwrite=ow)

    ###################################
    ## Exporter graphique du résultat
    if(export) {
      grDevices::png(file=stringr::str_c(tempdir(), "/", fname, ".png"), width = 8, height = 6, units='in', res=300, bg = "white")
      raster::plot(raster::crop(raster.invert(UD2rast(mwKern)),
                hone.extent(spatpol=ci_contours[ci_contours$prob==1,], r=spatres*2)),
           main=ifelse(is.null(titre), fname, titre))
      plot(ci_contours[ci_contours$prob %in% c(0.1,0.5,0.75,0.95,1.00),], border='lightgrey', lwd=0.5, add=T)
      plot(core.area.spdf, border="red", lwd=0.5, add=T)
      grDevices::dev.off()

      zip(zipfile = str_c('output/', fname), files = str_c(tempdir(), '/', dir(tempdir())))
    }

    message("Analysis complete.")

  })

  message(stringr::str_c("Processing time: ", round(unname(ptime[3]), digits=2), ' seconds'))

  if(export) saveRDS(ptime, str_c(tempdir(), '/ptime.rds'))

  if(obj.ret) return(list(mwKern=outKern, ci_contours=ci_contours, coreArea=core.area.spdf, herd.grid=herd.grid))

}

