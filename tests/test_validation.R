Require::Require(c('wmKDE','dplyr','sf','terra'))

train <- bind_rows(data.frame(idName = 'id1', col=8, x = rnorm(n=500), y = rnorm(n=500)),
                   data.frame(idName = 'id2', col=4, x = rnorm(n=100, mean=-3.5), y = rnorm(n=100, mean=-3.5))) %>%
  st_as_sf(., coords=c('x','y'), crs=st_crs(32198))

validSet <- function() {
  bind_rows(data.frame(idName = 'id1', col=8, x = rnorm(n=500), y = rnorm(n=500)),
            data.frame(idName = 'id2', col=4, x = rnorm(n=100, mean=-3.5), y = rnorm(n=100, mean=-3.5))) %>%
    st_as_sf(., coords=c('x','y'), crs=st_crs(32198))
}

## pooled kernel
bkern <- wmKDE(train, avg = F, spatres = 0.1, verbose = F, trim = T)$wmKDE

## mean kernel (different results)
mkern <- resample(wmKDE(train, id = 'id', ncores = 1, avg = T, spatres = 0.1, verbose = F, trim = T)$wmKDE, bkern2)



## validate
valid = validSet()
F1score(pred = bkern, obs = valid, id = 'idName', threshVal = 70, binWidth = 100)
F1score(pred = mkern, obs = valid, id = 'idName', threshVal = 70, binWidth = 100)

valid = validSet()
F1score(pred = bkern, obs = group_by(valid, idName) %>% slice_head(n = 100),
id = 'idName', threshVal = 70, binWidth = 100)
F1score(pred = mkern, obs = group_by(valid, idName) %>% slice_head(n = 100),
id = 'idName', threshVal = 70, binWidth = 100)

