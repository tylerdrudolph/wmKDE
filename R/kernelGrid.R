#' Define a spatial grid corresponding to the distribution area of the population of interest.
#'
#' Function called by mwKDE to define the common grid over which individual UDs are estimated.
#'
#' @param locs data.frame with columns 'x' and 'y' corresponding to projected spatial coordinates of relocations
#' @param exp.range single integer defining the range expansion parameter ensuring sufficient area is available to fully estimate the population UD (3 by default).
#' @param cell.size desired cell resolution (dx=dy)
#'
#' @return Generates range.x and grid.size arguments required by KernSmooth::bkde2D.
#' @export
#'
#' @examples
#'
kernelGrid <- function(locs, exp.range = 3, cell.size) {
  
  if(is.null(cell.size)) stop("must specify cell size!")
  
  r <- terra::rast(xmin = terra::xmin(ext(locs)) - exp.range * stats::sd(sf::st_coordinates(locs)[,1]), 
                   xmax = terra::xmax(ext(locs)) + exp.range * stats::sd(sf::st_coordinates(locs)[,1]),
                   ymin = terra::ymin(ext(locs)) - exp.range * stats::sd(sf::st_coordinates(locs)[,2]), 
                   ymax = terra::ymax(ext(locs)) + exp.range * stats::sd(sf::st_coordinates(locs)[,2]),
                   res = cell.size,
                   crs = terra::crs(locs))
              
  outmat <- cbind.data.frame(x = terra::xFromCol(r, col = 1:max(dim(r)[-3])), y = terra::yFromRow(r, row = 1:max(dim(r)[-3])))
  grid.size <- c(sum(!is.na(outmat$x)), sum(!is.na(outmat$y)))
  range.x <- list(c(min(outmat$x, na.rm=T), max(outmat$x, na.rm=T)), c(min(outmat$y, na.rm=T), max(outmat$y, na.rm=T)))

  # if(sum(is.na(locs$x) | is.na(locs$y)) > 0) message(sum(is.na(locs$x) | is.na(locs$y)), " cases of missing coordinates excluded")
  # locs <- locs[!is.na(locs$x) & !is.na(locs$y),]
  # 
  # min.grid.x = min(locs$x) - exp.range * stats::sd(locs$x)
  # max.grid.x = max(locs$x) + exp.range * stats::sd(locs$x)
  # xlength <- ceiling((max.grid.x - min.grid.x) / cell.size)
  # 
  # min.grid.y = min(locs$y) - exp.range * stats::sd(locs$y)
  # max.grid.y = max(locs$y) + exp.range * stats::sd(locs$y)
  # ylength <- ceiling((max.grid.y - min.grid.y) / cell.size)
  # 
  # X = seq(from=min.grid.x, by=cell.size, length.out=max(xlength, ylength))
  # max.grid.x <- max(X)
  # Y = seq(from=min.grid.y, by=cell.size, length.out=max(xlength, ylength))
  # max.grid.y <- max(Y)
  # 
  # i <- which(X >= min.grid.x - 0.5 * cell.size & X < max.grid.x + 0.5 * cell.size &
  #              Y >= min.grid.y - 0.5 * cell.size & Y < max.grid.y + 0.5 * cell.size)
  # 
  # xi <- findInterval(X[i], X - 0.5 * cell.size)
  # yi <- findInterval(Y[i], Y - 0.5 * cell.size)
  # 
  # outmat <- cbind.data.frame(x=X[xi], y=Y[yi])
  # grid.size <- rep(nrow(outmat), 2)
  # range.x <- list(c(min(outmat$x, na.rm=T), max(outmat$x, na.rm=T)), c(min(outmat$y, na.rm=T), max(outmat$y, na.rm=T)))

  return(list(grid.size = grid.size, range.x = range.x, mat = outmat, r = r))

}
