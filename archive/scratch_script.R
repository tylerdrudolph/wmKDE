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

# readRDS('/mnt/DATA/MFFP/GitHub/Outil-caribou-wmKDE/data/cartab.rds') %>%
readRDS('archive/MAJ_2024/cartab_2024-01-26.rds') %>%
  # rename(IDANIMAL = IDAnimal,
         # AN = An) %>%
  # mutate(POPULATION = as.character(Pop)) %>%
  rename(POPULATION = POP_ASSIGN) %>%
  # filter(!str_detect(POPULATION, 'Acquisition_BJ')) %>%
  filter(!POPULATION %in% c('Acquisition_BJ','Migrateur-TRAF','Indetermine')) %>%
  st_as_sf(., coords = c('X_LAMB','Y_LAMB'), crs = 32198) %>%
  st_transform(6623) %>%
  st_write(., dsn = 'input/cartab.gpkg', delete_layer = T)

cell.size = 250
uPop <- (st_read('input/cartab.gpkg', query = 'select distinct(POPULATION) from cartab', quiet = TRUE) %>% pull)[-c(1:2)]

for(pop in uPop) {

  cat('Population =', pop, '\n')

  ptab <- st_read('input/cartab.gpkg', query = paste0("select * from cartab where POPULATION is '", pop, "'"), quiet = TRUE)

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

  ## 3 a) No spatial weighting (k0)
  # setNames(
    wmKDE(dtrain, id = 'IDANIMAL', avg = T, trim = F, popGrid = hgrid, ktype = 'iso',
                 write2file = TRUE,
                 writeDir = 'tmp')
    # $wmKDE, 'mUD')
    
    ## c) Weighted mean population kernel
    # setNames(
      wmKDE(dtrain, id = 'IDANIMAL', avg = T, trim = F, popGrid = hgrid, ktype = 'iso',
                   udw = group_by(dtrain, IDANIMAL) %>%
                     summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
                     st_drop_geometry %>%
                     select(IDANIMAL, w),
            write2file = T, 
            writeDir = 'tmp')
      # $wmKDE, 'wmUD')

  
  tUD0_list <- c(setNames(sapp(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                                    id = dtrain$IDANIMAL,
                                    wts = NULL,
                                    sproj = st_crs(dtrain),
                                    bwGlobal = T,
                                    userGrid = hgrid,
                                    ncores = parallelly::availableCores() - 2, verbose = T),
                               fhat2confin), unique(dtrain$IDANIMAL)),
                          ## b) Mean population kernel
                          setNames(wmKDE(dtrain, id = 'IDANIMAL', avg = T, trim = F, popGrid = hgrid, ktype = 'iso')$wmKDE, 'mUD'),
                          ## c) Weighted mean population kernel
                          setNames(wmKDE(dtrain, id = 'IDANIMAL', avg = T, trim = F, popGrid = hgrid, ktype = 'iso',
                                         udw = group_by(dtrain, IDANIMAL) %>%
                                           summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
                                           st_drop_geometry %>%
                                           select(IDANIMAL, w))$wmKDE, 'wmUD'))
  
  writeRaster(tUD0_list, filename = paste0('output/tUD0_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  rm(tUD0_list)
  gc()

  ##########################################################################################
  # 3 b i) Poids = delta values (max = 317) * pop-specific Manly coefficients (k1)
  tUD1_list <- c(setNames(sapp(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                                    id = dtrain$IDANIMAL,
                                    wts = dtrain$p1,
                                    sproj = st_crs(dtrain),
                                    bwGlobal = T,
                                    userGrid = hgrid,
                                    ncores = 7, verbose = T),
                               fhat2confin), unique(dtrain$IDANIMAL)),
                 ## 3 b ii) Add mean and weighted mean population kernels (k1)
                 setNames(wmKDE(dtrain, id = 'IDANIMAL', spw = 'p1', avg = T, trim = F, popGrid = hgrid, ktype = 'iso')$wmKDE, 'mUD'),
                 setNames(wmKDE(dtrain, id = 'IDANIMAL', spw = 'p1', avg = T, trim = F, popGrid = hgrid, ktype = 'iso',
                                udw = group_by(dtrain, IDANIMAL) %>%
                                  summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
                                  st_drop_geometry %>%
                                  select(IDANIMAL, w))$wmKDE, 'wmUD'))
  
  writeRaster(tUD1_list, filename = paste0('output/tUD1_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  rm(tUD1_list)
  gc()

  ############################################################3
  # 3 c i) Poids binaires = 'tolérance zéro' (k2)
  tUD2_list <- c(setNames(sapp(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                                    id = dtrain$IDANIMAL,
                                    wts = dtrain$p2,
                                    sproj = st_crs(dtrain),
                                    bwGlobal = T,
                                    userGrid = hgrid,
                                    ncores = 7, verbose = T),
                               fhat2confin), unique(dtrain$IDANIMAL)),
                 ## 3 b ii) Add mean and weighted mean population kernels (k2)
                 setNames(wmKDE(dtrain, id = 'IDANIMAL', spw = 'p2', avg = T, trim = F, popGrid = hgrid, ktype = 'iso')$wmKDE, 'mUD'),
                 setNames(wmKDE(dtrain, id = 'IDANIMAL', spw = 'p2', avg = T, trim = F, popGrid = hgrid, ktype = 'iso',
                                udw = group_by(dtrain, IDANIMAL) %>%
                                  summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
                                  st_drop_geometry %>%
                                  select(IDANIMAL, w))$wmKDE, 'wmUD'))
  
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
  vUD_list <- c(setNames(sapp(bKDE(xy = xy,
                                   id = vtrain$IDANIMAL, 
                                   wts = NULL, 
                                   sproj = st_crs(vtrain),
                                   bwGlobal = T,
                                   userGrid = hgrid, 
                                   ncores = 7, 
                                   verbose = T), 
                              fhat2confin), unique(vtrain$IDANIMAL)),
                ## 4 b) Calculate and append mean (unweighted) population kernel
                setNames(wmKDE(vtrain, id = 'IDANIMAL', avg = T, trim = F, popGrid = hgrid, ktype = 'iso')$wmKDE, 'mUD'))
  
  ## 5) Write to file
  writeRaster(vUD_list, filename = paste0('output/vUD_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  rm(vUD_list)
  rm(ptab)
  gc()
  
}




#####################

