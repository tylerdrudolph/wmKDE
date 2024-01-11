# Require::Require(c('wmKDE','dplyr','sf','terra'))

## 1) same n + sd

## 2) same n, different sd

## 3) different n, same sd

## 4) different n, different sd


# train <- bind_rows(data.frame(idName = 'id1', col=8, x = rnorm(n=500), y = rnorm(n=500)),
#                    data.frame(idName = 'id2', col=4, x = rnorm(n=100, mean=-3.5), y = rnorm(n=100, mean=-3.5))) %>%
#   st_as_sf(., coords=c('x','y'), crs=st_crs(32198))

simSet <- function(nID = 5, equal.n = TRUE, equal.sd = TRUE, multimodal = F) {
  
  if(equal.n) nobs <- sample(seq(50, 1000, 50), 1)
  if(equal.sd) sdObs <- sample(seq(1, 2, 0.05), 1)
  
  do.call(rbind, lapply(1:nID, function(i) {
    
    if(!equal.n) nobs <- sample(seq(50, 1000, 50), 1)
    if(!equal.sd) sdObs <- sample(seq(1, 2, 0.05), 1)
    
    x <- data.frame(idName = paste0('id', i), 
               x = rnorm(n = nobs, 
                         mean = runif(n = 1, min = -1, max = 1),
                         sd = sdObs),
               y = rnorm(n = nobs, 
                         mean = runif(n = 1, min = -1, max = 1),
                         sd = sdObs)) %>%
      st_as_sf(., coords = c('x','y'), crs = st_crs(32198))
    
    x$set <- 'train'
    x$set[sample.int(n = nrow(x), size = nrow(x) * 0.2, replace = F)] <- 'valid'
    
    return(x)
    
  }))
  
}

vMix <- lapply(1:50, function(i) {

  cat(i, '\n')

  ## pooled kernel
  x <- simSet(equal.n = F, equal.sd = F)
  spatres <- 0.5
  bkern <- wmKDE(filter(x, set=='train'), avg = F, spatres = spatres, verbose = F, trim = F)$wmKDE

  ## mean kernel (zscaling improves fit)
  mkern0 <- resample(wmKDE(filter(x, set=='train'), id = 'idName', ncores = 1, zscale = F, avg = T, spatres = spatres, verbose = F, trim = F)$wmKDE, bkern)
  mkern1 <- resample(wmKDE(filter(x, set=='train'), id = 'idName', ncores = 1, zscale = T, avg = T, spatres = spatres, verbose = F, trim = F)$wmKDE, bkern)

  ## validation kernel
  vkern <- resample(wmKDE(filter(x, set=='valid'), id = 'idName', ncores = 1, avg = T, spatres = spatres, verbose = F, trim = F)$wmKDE, bkern)

  list(F1score.rawPts = setNames(c(F1score(pred = bkern, obs = filter(x, set == 'valid'))$FScore,
                                   F1score(pred = mkern0, obs = filter(x, set == 'valid'))$FScore,
                                   F1score(pred = mkern1, obs = filter(x, set == 'valid'))$FScore),
                                 c('bkern','mkern0','mkern1')),

       F1score.wtdPts = setNames(c(F1score(pred = bkern, obs = filter(x, set == 'valid'), id = 'idName') %>%
                                     filter(id == 'wtdMean') %>%
                                     pull(FScore),
                                   F1score(pred = mkern0, obs = filter(x, set == 'valid'), id = 'idName') %>%
                                     filter(id == 'wtdMean') %>%
                                     pull(FScore),
                                   F1score(pred = mkern1, obs = filter(x, set == 'valid'), id = 'idName') %>%
                                     filter(id == 'wtdMean') %>%
                                     pull(FScore)),
                                 c('bkern','mkern0','mkern1')),

       F1score.mkerns = setNames(c(F1score(pred = bkern, obs = vkern)$FScore,
                                   F1score(pred = mkern0, obs = vkern)$FScore,
                                   F1score(pred = mkern1, obs = vkern)$FScore),
                                 c('bkern','mkern0','mkern1')))

})

# sort(sapply(vNull, function(x) names(which.max(x$F1score.wtdPts))) %>% table(), decreasing = T)
# sort(sapply(vNull, function(x) names(which.max(x$F1score.rawPts))) %>% table(), decreasing = T)
# sort(sapply(vNull, function(x) names(which.max(x$F1score.mkerns))) %>% table(), decreasing = T)
# 
# sort(sapply(vN, function(x) names(which.max(x$F1score.wtdPts))) %>% table(), decreasing = T)
# sort(sapply(vN, function(x) names(which.max(x$F1score.rawPts))) %>% table(), decreasing = T)
# sort(sapply(vN, function(x) names(which.max(x$F1score.mkerns))) %>% table(), decreasing = T)
# 
# sort(sapply(vSD, function(x) names(which.max(x$F1score.wtdPts))) %>% table(), decreasing = T)
# sort(sapply(vSD, function(x) names(which.max(x$F1score.rawPts))) %>% table(), decreasing = T)
# sort(sapply(vSD, function(x) names(which.max(x$F1score.mkerns))) %>% table(), decreasing = T)
# 
# sort(sapply(vMix, function(x) names(which.max(x$F1score.wtdPts))) %>% table(), decreasing = T)
# sort(sapply(vMix, function(x) names(which.max(x$F1score.rawPts))) %>% table(), decreasing = T)
# sort(sapply(vMix, function(x) names(which.max(x$F1score.mkerns))) %>% table(), decreasing = T)
# 
saveRDS(setNames(list(vNull, vN, vSD, vMix), c('vNull','vN','vSD','vMix')), 
        file = 'output/simResults.rds')

simOut <- do.call(rbind, lapply(c('vNull','vN','vSD','vMix'), function(i) {
  simout <- do.call(rbind, lapply(readRDS('output/simResults.rds')[[i]], function(x) {
    cbind.data.frame(simType = i, 
                     kType = c('bkern','mkern0','mkern1'), t(do.call(rbind, x)))
  })) %>% 
    tibble::remove_rownames()
})) %>%
  group_by(simType, kType) %>%
  summarize(rawPts = mean(F1score.rawPts),
            wtdPts = mean(F1score.wtdPts),
            mkerns = mean(F1score.mkerns), .groups = 'drop') %>%
  
  




# 
# 
# ## 1) validate using equal # observations sampled from each unique ID)
# valid = validSet()
# F1score(pred = bkern, obs = slice_sample(valid, n = 100, by = 'idName', replace = F),
#         id = 'idName', threshVal = 70, binWidth = 100)
# F1score(pred = mkern, obs = slice_sample(valid, n = 100, by = 'idName', replace = F),
#         id = 'idName', threshVal = 70, binWidth = 100)
# 
# ## 2) Validate using different # observations that are inversely weighted
# valid = validSet()
# F1score(pred = bkern, obs = valid, id = 'idName', threshVal = 70, binWidth = 100)
# F1score(pred = mkern, obs = valid, id = 'idName', threshVal = 70, binWidth = 100)
# 
# ## 3) validate using resampled observations of equivalent number across unique IDs
# valid = validSet()
# b <- sapply(1:50, function(i) F1score(pred = bkern, obs = slice_sample(valid, n = 25, by = 'idName', replace = T),
#                                       id = 'idName', threshVal = 70, binWidth = 100) %>%
#               filter(id == 'wtdMean') %>%
#               pull(FScore))
# m <- sapply(1:50, function(i) F1score(pred = mkern, obs = slice_sample(valid, n = 25, by = 'idName', replace = T),
#                                       id = 'idName', threshVal = 70, binWidth = 100) %>%
#               filter(id == 'wtdMean') %>%
#               pull(FScore))
# plot(density(m), col='red')
# lines(density(b))
