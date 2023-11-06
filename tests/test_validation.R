# Require::Require(c('wmKDE','dplyr','sf','terra'))
# 
# locs2 <- bind_rows(data.frame(id = 1, col=8, x = rnorm(n=500), y = rnorm(n=500)),
#                    data.frame(id = 2, col=4, x = rnorm(n=100, mean=-3.5), y = rnorm(n=100, mean=-3.5))) %>%
#   st_as_sf(., coords=c('x','y'), crs=st_crs(32198))
# 
# ## pooled kernel
# bkern2 <- wmKDE(locs2, avg=F, spatres=0.1, verbose=F, trim = F)
# 
# ## mean kernel (different results)
# mkern2 <- wmKDE(locs2, id='id', ncores=1, avg=T, spatres=0.1, verbose=F, trim = F)
# 
# ## validate
# F1score(bkern2$wmKDE, mkern2$wmKDE, threshVal = 90)