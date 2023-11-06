# Require::Require(c('wmKDE','dplyr','sf','terra'))
# 
locs2 <- bind_rows(data.frame(id = 1, col=8, x = rnorm(n=500), y = rnorm(n=500)),
                   data.frame(id = 2, col=4, x = rnorm(n=100, mean=-3.5), y = rnorm(n=100, mean=-3.5))) %>%
  st_as_sf(., coords=c('x','y'), crs=st_crs(32198))

## pooled kernel
bkern2 <- wmKDE(locs2, avg=F, spatres=0.1, verbose=F, trim = F)

## mean kernel (different results)
mkern2 <- wmKDE(locs2, id='id', ncores=1, avg=T, spatres=0.1, verbose=F, trim = F)

# ## validate
# relevantCats = c('(0.7–0.8]','(0.8–0.9]','(0.9–1]')
# # relevantCats <- '(0.5–1]'
# nbins = 10
# 
# cuts = c(0, 0.6, 1)
# v <- crosstab(c(classify(setNames(bkern2$wmKDE, 'pred'), 
#                          rcl = cuts), #seq(0, 1, 1/nbins)),
#                 classify(setNames(mkern2$wmKDE, 'true'), 
#                          rcl = cuts))) #seq(0, 1, 1/nbins))))
# recall(v, relevant = row.names(v)[2])#relevantCats)
# precision(v, relevant = row.names(v)[2])#relevantCats)
# F_meas(v, relevant = row.names(v)[2])#relevantCats, beta = 1)
# 
# caret::F_meas(v)
# caret::confusionMatrix(v)
# 
# layerCor(c(bkern2$wmKDE, mkern2$wmKDE),
#          fun = 'pearson')
# 
