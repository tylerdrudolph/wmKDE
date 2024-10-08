# # Require::Require(c('dplyr','stringr','sf','latticeDensity','terra','wmKDE'))
# # setwd('/mnt/DATA/Huron-Wendat/SAFO-SAOQ')
# # 
# # x <- readRDS('inputs/detections_rarefaction_lac_Martel.rds') %>%
# #   filter(Mois %in% c('Sep','Oct','Nov')) %>%
# #   # fixer identifiant unique par semaine
# #   mutate(week = lubridate::week(DATETIME))
# # 
# # id <- 'DETECTEDID'
# # 
# # bath <- st_read('inputs/raw/Bathy_Martel', 'bathy_martel') %>%
# #   st_cast(., 'POLYGON') %>%
# #   st_buffer(., 0.0) %>%
# #   st_union(., by_feature = F)
# # 
# # spatres = 10
# # 
# # writeDir = 'outputs'
# # 
# # latticeKernels <- readRDS(paste0('outputs/latticeKernels_byTagWeek_automne_', 10, 'm.rds'))

Require::Require(c('dplyr','sf','terra','stringr','wmKDE'))
cattab <- readRDS('archive/cattab.rds')
msrXpop <- readRDS('archive/newdata.msr.msrxpop.coefs.rds')$coefs
# msrXpop <- do.call(rbind, lapply(1:length(msrXpop), function(i) {
#   bind_cols(POP = names(msrXpop)[i], msrXpop[[i]])
# }))

cartab <- readRDS('archive/MAJ_2024/cartab_2024-01-26.rds') %>%
  rename(POPULATION = POP_ASSIGN,
         dburn = burn,
         dcut = cut,
         droad = road) %>%
  filter(!POPULATION %in% c('Acquisition_BJ','Migrateur-TRAF','Indetermine')) %>%
  mutate(uval = paste(ifelse(dburn > 0, 1, 0), 
                      ifelse(dcut > 0, 1, 0), 
                      ifelse(droad > 0, 1, 0)),
         cat = cattab$cat[match(uval, cattab$uval)])

do.call(rbind, lapply(split(cartab, cartab$POPULATION), function(x) {
  msr <- msrXpop[[unique(x$POPULATION)]]
  x$wi <- msr$wi[match(x$cat, msr$cat)]
  x
})) %>%
  ## weights based on observed change * appropriate Manly coefficient (pop-specific)
  mutate(p1 = ifelse(dburn + dcut + droad > 317, 317, dburn + dcut + droad),
         p1 = wmKDE::range01(ifelse(p1 == 0, 317, p1 * wi)),
         ## binary weights where any environmental change since relocation merits 0, otherwise 1
         p2 = ifelse(dcut + dburn + droad == 0, 1, 0)) %>%
  st_as_sf(., coords = c('X_LAMB','Y_LAMB'), crs = 32198) %>%
  st_transform(6623) %>%
  st_write(., dsn = 'input/cartab.gpkg', delete_layer = T)

cell.size = 250
uPop <- (st_read('input/cartab.gpkg', query = 'select distinct(POPULATION) from cartab', quiet = TRUE) %>% pull)[-c(1:2)]

for(pop in uPop) {

  cat('Population =', pop, '\n')
  
  ptab <- st_read('input/cartab.gpkg', 
                  query = paste0("select * from cartab where POPULATION is '", pop, "'"), 
                  quiet = TRUE)
  
  ## Ensure sufficient locs per IDANIMAL
  ntab <- table(ptab$IDANIMAL)
  ptab <- filter(ptab, IDANIMAL %in% names(ntab)[ntab >= 50])
  rm(ntab)
  
  ## 2 c) Define common grid parameters
  hgrid <- kernelGrid(locs = ptab,
                      cell.size = cell.size)
  
  ## 2 d) Estimate simple kernel using all points in the population
  dtrain <- filter(ptab, AN < 2022)
  skUD <- setNames(wmKDE(dtrain, id = 'IDANIMAL', avg = F, trim = F, popGrid = hgrid, ktype = 'iso')$wmKDE, 'skUD')
  writeRaster(skUD, filename = paste0('output/skUD_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  rm(skUD)
  gc()
  
  ##################################################################################
  ## 3 a) No spatial weighting (k0)
  sapp(rast(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                 id = dtrain$IDANIMAL,
                 wts = NULL,
                 sproj = st_crs(dtrain),
                 bwGlobal = T,
                 userGrid = hgrid,
                 ncores = 6, 
                 verbose = T)),
       fun = fhat2confin, 
       filename = 'tmp/tUD0_1.tif')
  gc()
  
  ## 3 b) Mean population kernel
  wmKDE(dtrain, id = 'IDANIMAL', avg = T, trim = F, ncores = 6, popGrid = hgrid, 
        ktype = 'iso', write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD0_2', retObj = FALSE)
  gc()
  
  ## 3 c) Weighted mean population kernel
  wmKDE(dtrain, id = 'IDANIMAL', avg = T, trim = F, ncores = 6, popGrid = hgrid, ktype = 'iso',
        udw = group_by(dtrain, IDANIMAL) %>%
          summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
          st_drop_geometry %>%
          select(IDANIMAL, w), 
        write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD0_3', retObj = FALSE)
  gc()
  
  ## 3 d) Write to file
  f <- list.files('tmp', pattern = 'tUD0_', full.names = T)
  f <- f[str_sub(f, start = -4L)=='.tif']
  writeRaster(c(setNames(rast('tmp/tUD0_1.tif'), unique(dtrain$IDANIMAL)),
                setNames(rast(f[str_detect(f, 'tUD0_2')]), 'mUD'),
                setNames(rast(f[str_detect(f, 'tUD0_3')]), 'wmUD')),
              filename = paste0('output/tUD0_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  gc()

  ##########################################################################################
  # 4 a) Poids = delta values (max = 317) * pop-specific Manly coefficients (k1)
  sapp(rast(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                 id = dtrain$IDANIMAL,
                 wts = dtrain$p1,
                 sproj = st_crs(dtrain),
                 bwGlobal = T,
                 userGrid = hgrid,
                 ncores = 6, verbose = T)),
       fun = fhat2confin, filename = 'tmp/tUD1_1.tif')
  gc()
  
  ## 4 b) Mean population kernel
  wmKDE(dtrain, id = 'IDANIMAL', spw = 'p1', avg = T, 
        trim = F, ncores = 6, popGrid = hgrid, ktype = 'iso', 
        write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD1_2', 
        retObj = FALSE)
  gc()
  
  ## 4 c) Weighted mean population kernel
  wmKDE(dtrain, id = 'IDANIMAL', spw = 'p1', avg = T, trim = F, ncores = 6, popGrid = hgrid, ktype = 'iso',
        udw = group_by(dtrain, IDANIMAL) %>%
          summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
          st_drop_geometry %>%
          select(IDANIMAL, w), 
        write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD1_3', retObj = FALSE)
  gc()
  
  ## 4 d) Write to file
  f <- list.files('tmp', pattern = 'tUD1_', full.names = T)
  f <- f[str_sub(f, start = -4L)=='.tif']
  writeRaster(c(setNames(rast('tmp/tUD1_1.tif'), unique(dtrain$IDANIMAL)),
                setNames(rast(f[str_detect(f, 'tUD1_2')]), 'mUD'),
                setNames(rast(f[str_detect(f, 'tUD1_3')]), 'wmUD')),
              filename = paste0('output/tUD1_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  gc()

  ############################################################3
  # 3 c i) Poids binaires = 'tolérance zéro' (k2)
  tUD2_list <- c(setNames(sapp(rast(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                                         id = dtrain$IDANIMAL,
                                         wts = dtrain$p2,
                                         sproj = st_crs(dtrain),
                                         bwGlobal = T,
                                         userGrid = hgrid,
                                         ncores = 6, verbose = T)),
                               fun = fhat2confin, filename = tempfile(fileext = '.tif')), unique(dtrain$IDANIMAL)),
                 ## 3 b ii) Add mean and weighted mean population kernels (k2)
                 setNames(wmKDE(dtrain, id = 'IDANIMAL', spw = 'p2', avg = T, trim = F, 
                                ncores = 6, popGrid = hgrid, ktype = 'iso', write2file = TRUE, writeDir = 'tmp')$wmKDE, 'mUD'),
                 setNames(wmKDE(dtrain, id = 'IDANIMAL', spw = 'p2', avg = T, trim = F, 
                                ncores = 6, popGrid = hgrid, ktype = 'iso',
                                udw = group_by(dtrain, IDANIMAL) %>%
                                  summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
                                  st_drop_geometry %>%
                                  select(IDANIMAL, w), 
                                write2file = TRUE, writeDir = 'tmp')$wmKDE, 'wmUD'))
  
  writeRaster(tUD2_list, filename = paste0('output/tUD2_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  rm(tUD2_list)
  gc()

  ##############################################################
  # 4 a) Generate validation kernels and write to file
  vtrain <- filter(ptab, AN >= 2022)

  ## Ensure sufficient locs per IDANIMAL
  ntab <- table(vtrain$IDANIMAL)
  vtrain <- filter(vtrain, IDANIMAL %in% names(ntab)[ntab >= 50])
  rm(ntab)

  xy <- as.data.frame(st_coordinates(vtrain)) %>% rename(x=X, y=Y)
  vUD_list <- c(setNames(sapp(rast(bKDE(xy = xy,
                                   id = vtrain$IDANIMAL, 
                                   wts = NULL, 
                                   sproj = st_crs(vtrain),
                                   bwGlobal = T,
                                   userGrid = hgrid, 
                                   ncores = 4, 
                                   verbose = T)), 
                              fun = fhat2confin,
                              filename = tempfile(fileext = '.tif')), unique(vtrain$IDANIMAL)),
                ## 4 b) Calculate and append mean (unweighted) population kernel
                setNames(wmKDE(vtrain, id = 'IDANIMAL', avg = T, trim = F, popGrid = hgrid, ktype = 'iso')$wmKDE, 'mUD'))
  
  ## 5) Write to file
  writeRaster(vUD_list, filename = paste0('output/vUD_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  rm(vUD_list)
  rm(ptab)
  gc()
  
}




#####################

