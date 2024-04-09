Require::Require(c('dplyr','sf','terra','stringr','wmKDE'))
devtools::load_all()
cattab <- readRDS('archive/cattab.rds')
msrXpop <- readRDS('archive/newdata.msr.msrxpop.coefs.rds')$coefs
pertval <- data.table::setDT(readRDS('archive/MAJ_2024/Perturbations/extract_allyear.rds'))
idref <- read.csv('archive/MAJ_2024/Validation_methode/ID_cartab_pourValidation.csv', sep = ';')

## 1) ----
cartab <- readRDS('archive/MAJ_2024/cartab_2024-01-26.rds') %>%
  rename(POPULATION = POP_ASSIGN) %>% 
  filter(!POPULATION %in% c('Acquisition_BJ','Migrateur-TRAF','Indetermine'))
cartab$dburn <- (pertval[Donnees=='Feux'][match(cartab$ID_TELEMETRIE,ID_TELEMETRIE),'Annee_2021'] %>% pull()) - cartab$burn
cartab$dcut <- (pertval[Donnees=='Coupe'][match(cartab$ID_TELEMETRIE,ID_TELEMETRIE),'Annee_2021'] %>% pull()) - cartab$cut
cartab$droad <- (pertval[Donnees=='Route'][match(cartab$ID_TELEMETRIE,ID_TELEMETRIE),'Annee_2021'] %>% pull()) - cartab$road

cartab <- mutate(cartab,
                 uval = paste(ifelse(dburn > 0, 1, 0), 
                              ifelse(dcut > 0, 1, 0), 
                              ifelse(droad > 0, 1, 0)),
                 cat = cattab$cat[match(uval, cattab$uval)]) %>% 
  filter(IDANIMAL %in% unique(idref$IDAnimal[idref$n > 50]))

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

## 2) ----
cell.size = 250
uPop <- (st_read('input/cartab.gpkg', query = 'select distinct(POPULATION) from cartab', quiet = TRUE) %>% pull)[-c(1:2)]

for(pop in uPop) {

  cat('Population =', pop, '\n')
  
  ## 3 a) ----
  ptab <- st_read('input/cartab.gpkg', 
                  query = paste0("select * from cartab where POPULATION is '", pop, "' and IDANIMAL in(",
                                 str_c("'", str_c(idref[idref$Pop == pop,'IDAnimal'], collapse = "','"), "')")), 
                  quiet = TRUE)
  
  ## Ensure sufficient locs per IDANIMAL ----
  ntab <- table(ptab$IDANIMAL)
  ptab <- filter(ptab, IDANIMAL %in% names(ntab)[ntab >= 50])
  rm(ntab)
  
  ## 3 c) Define common grid parameters ----
  hgrid <- kernelGrid(locs = ptab,
                      cell.size = cell.size)

  ## 3 d) Estimate simple kernel using all points in the population ----
  dtrain <- filter(ptab, AN < 2021)
  skUD <- setNames(wmKDE(dtrain, id = 'IDANIMAL', avg = F, trim = F, popGrid = hgrid, ktype = 'iso')$wmKDE, 'skUD')
  writeRaster(skUD, filename = paste0('output/skUD_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  rm(skUD)
  gc()
  
  ##################################################################################
  ## 4 a) No spatial weighting (k0) ----
  sapp(rast(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                 id = dtrain$IDANIMAL,
                 wts = NULL,
                 sproj = st_crs(dtrain),
                 bwGlobal = T,
                 userGrid = hgrid,
                 ncores = 6, 
                 verbose = T)),
       fun = fhat2confin, 
       filename = 'tmp/tUD0_1.tif', overwrite = T)
  gc()
  
  ## 4 b) Mean population kernel ----
  wmKDE(dtrain, id = 'IDANIMAL', avg = T, trim = F, ncores = 6, popGrid = hgrid, 
        ktype = 'iso', write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD0_2', retObj = FALSE)
  gc()
  
  ## 4 c) Weighted mean population kernel ----
  wmKDE(dtrain, id = 'IDANIMAL', avg = T, trim = F, ncores = 6, popGrid = hgrid, ktype = 'iso',
        udw = group_by(dtrain, IDANIMAL) %>%
          summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
          st_drop_geometry %>%
          select(IDANIMAL, w), 
        write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD0_3', retObj = FALSE)
  gc()
  
  ## 4 d) Write to file ----
  f <- list.files('tmp', pattern = 'tUD0_', full.names = T)
  f <- f[str_sub(f, start = -4L)=='.tif']
  writeRaster(c(setNames(rast('tmp/tUD0_1.tif'), unique(dtrain$IDANIMAL)),
                setNames(rast(f[str_detect(f, 'tUD0_2')]), 'mUD'),
                setNames(rast(f[str_detect(f, 'tUD0_3')]), 'wmUD')),
              filename = paste0('output/tUD0_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  gc()
  file.remove(list.files('tmp', pattern = 'tUD0_', full.names = T))

  ##########################################################################################
  # 5 a) Poids = delta values (max = 317) * pop-specific Manly coefficients (k1) ----
  sapp(rast(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                 id = dtrain$IDANIMAL,
                 wts = dtrain$p1,
                 sproj = st_crs(dtrain),
                 bwGlobal = T,
                 userGrid = hgrid,
                 ncores = 6, verbose = T)),
       fun = fhat2confin, filename = 'tmp/tUD1_1.tif', overwrite = T)
  gc()
  
  ## 4 b) Mean population kernel ----
  wmKDE(dtrain, id = 'IDANIMAL', spw = 'p1', avg = T, 
        trim = F, ncores = 6, popGrid = hgrid, ktype = 'iso', 
        write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD1_2', 
        retObj = FALSE)
  gc()
  
  ## 5 c) Weighted mean population kernel ----
  wmKDE(dtrain, id = 'IDANIMAL', spw = 'p1', avg = T, trim = F, ncores = 6, popGrid = hgrid, ktype = 'iso',
        udw = group_by(dtrain, IDANIMAL) %>%
          summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
          st_drop_geometry %>%
          select(IDANIMAL, w), 
        write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD1_3', retObj = FALSE)
  gc()
  
  ## 5 d) Write to file ----
  f <- list.files('tmp', pattern = 'tUD1_', full.names = T)
  f <- f[str_sub(f, start = -4L)=='.tif']
  writeRaster(c(setNames(rast('tmp/tUD1_1.tif'), unique(dtrain$IDANIMAL)),
                setNames(rast(f[str_detect(f, 'tUD1_2')]), 'mUD'),
                setNames(rast(f[str_detect(f, 'tUD1_3')]), 'wmUD')),
              filename = paste0('output/tUD1_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  gc()
  file.remove(list.files('tmp', pattern = 'tUD1_', full.names = T))
  
  ############################################################3
  ## 6 a) Poids binaires = 'tolérance zéro' (k2) ----
  sapp(rast(bKDE(xy = as.data.frame(st_coordinates(dtrain)) %>% rename(x=X, y=Y),
                 id = dtrain$IDANIMAL,
                 wts = dtrain$p2,
                 sproj = st_crs(dtrain),
                 bwGlobal = T,
                 userGrid = hgrid,
                 ncores = 6, verbose = T)),
       fun = fhat2confin, filename = 'tmp/tUD2_1.tif', overwrite = T)
  gc()
  
  ## 6 b) Mean population kernel ----
  wmKDE(dtrain, id = 'IDANIMAL', spw = 'p2', avg = T, trim = F, 
        ncores = 6, popGrid = hgrid, ktype = 'iso', write2file = TRUE, 
        writeDir = 'tmp', fileTag = 'tUD2_2', retObj = FALSE)
  gc()

  ## 6 c) Weighted mean population kernel ----
  wmKDE(dtrain, id = 'IDANIMAL', spw = 'p2', avg = T, trim = F, 
        ncores = 6, popGrid = hgrid, ktype = 'iso',
        udw = group_by(dtrain, IDANIMAL) %>%
          summarize(w = 1 - (sum(dburn > 0 | dcut > 0 | droad > 0) / length(dburn))) %>%
          st_drop_geometry %>%
          select(IDANIMAL, w), 
        write2file = TRUE, writeDir = 'tmp', fileTag = 'tUD2_3', retObj = FALSE)
  gc()

  ## 6 d) Write to file ----
  f <- list.files('tmp', pattern = 'tUD2_', full.names = T)
  f <- f[str_sub(f, start = -4L)=='.tif']
  writeRaster(c(setNames(rast('tmp/tUD2_1.tif'), unique(dtrain$IDANIMAL)),
                setNames(rast(f[str_detect(f, 'tUD2_2')]), 'mUD'),
                setNames(rast(f[str_detect(f, 'tUD2_3')]), 'wmUD')),
              filename = paste0('output/tUD2_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  gc()
  file.remove(list.files('tmp', pattern = 'tUD2_', full.names = T))

  ##############################################################
  ## 7 a) Generate validation kernels and write to file ----
  vtrain <- filter(ptab, AN >= 2021)

  ## 7 b) Ensure sufficient locs per IDANIMAL ----
  ntab <- table(vtrain$IDANIMAL)
  vtrain <- filter(vtrain, IDANIMAL %in% names(ntab)[ntab >= 50])
  rm(ntab)

  xy <- as.data.frame(st_coordinates(vtrain)) %>% rename(x=X, y=Y)
  
  ## 7 c)  ----
  sapp(rast(bKDE(xy = xy,
                 id = vtrain$IDANIMAL, 
                 wts = NULL, 
                 sproj = st_crs(vtrain),
                 bwGlobal = T,
                 userGrid = hgrid, 
                 ncores = 6, 
                 verbose = T)), 
       fun = fhat2confin,
       filename = 'tmp/vUD_1.tif', overwrite = T)
  
  ## 7 d) Calculate and append mean (unweighted) population kernel ----
  wmKDE(vtrain, id = 'IDANIMAL', avg = T, trim = F, popGrid = hgrid, ktype = 'iso',
        write2file = TRUE, writeDir = 'tmp', fileTag = 'vUD_2', retObj = FALSE)

  ## 7 e) Write to file ----
  f <- list.files('tmp', pattern = 'vUD_', full.names = T)
  f <- f[str_sub(f, start = -4L)=='.tif']
  writeRaster(c(setNames(rast('tmp/vUD_1.tif'), unique(vtrain$IDANIMAL)),
                setNames(rast(f[str_detect(f, 'vUD_2')]), 'mUD')),
              filename = paste0('output/vUD_', pop, '_', cell.size, 'm.tif'), overwrite = T)
  gc()
  file.remove(list.files('tmp', pattern = 'vUD_', full.names = T))
}




## Test outcomes (e.g.)
x1 <- F1score(pred = rast('output/vUD_Caniapiscau_250m.tif', lyr = 'mUD'),
              obs = rast('output/tUD0_Caniapiscau_250m.tif', lyr = 'mUD'))

x2 <- F1score(pred = rast('output/vUD_Caniapiscau_250m.tif', lyr = 'mUD'),
              obs = rast('output/tUD1_Caniapiscau_250m.tif', lyr = 'mUD'))

x3 <- F1score(pred = rast('output/vUD_Caniapiscau_250m.tif', lyr = 'mUD'),
              obs = rast('output/tUD2_Caniapiscau_250m.tif', lyr = 'mUD'))


y1 <- F1score(pred = rast('output/vUD_Caniapiscau_250m.tif', lyr = 'mUD'),
              obs = vtrain, 
              id = 'IDANIMAL')

                                          

