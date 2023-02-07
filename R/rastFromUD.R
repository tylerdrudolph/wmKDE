rastFromUD <- function(UD, wkt=st_crs(32198)$wkt) {

  return(flip(rast(x = t(UD$fhat),
                   extent = ext(min(UD$x1), max(UD$x1), min(UD$x2), max(UD$x2)),
                   crs=wkt)))

  #
  # rast(nrows = nrow(UD$fhat)+1,
  # ncols = ncol(UD$fhat) + 1,
  # xmin = min(UD$x1),
  # xmax = max(UD$x1),
  # ymin = min(UD$x2),
  # ymax = max(UD$x2),
  # crs = st_crs(32198)$wkt,
  # resolution = c(UD$x1[2] - UD$x1[1],
  #                UD$x2[2] - UD$x2[1]))#))

}



# xygrid <- expand.grid(x=UD$x1, y=UD$x2)
# xygrid$z <- UD$fhat[xygrid$x, xygrid$y]
#
# xyindex$z <- UD$fhat[xyindex$x, xyindex$y]
# nrow(xyindex)
