#' Redistribute UD probability densities proportionally so as to ensure that the sum of cell volumes equals one.
#'
#' @param r spatRaster containing kernel density probabilities (utilization distribution)
#' @param out integer parameter specifying whether output cell values are aspatial (out=1) or volume-based (out=2).
#'
#' @return spatRaster with probabilities adjusted such that corresponding volumes sum to exactly 1
#' @export
#'
finetune <- function(r, tin = 'prob', tout = 'prob') {
  
  if(tin == 'vol') {
    r <- r / prod(res(r))
  }
  cell.area <- prod(res(r))
  v <- r * cell.area
  tot <- unlist(terra::global(v, 'sum'))
  
  if(tot < 1) {
    if(1 - tot > .005) warning("Sum of cell values not equal to 1 (", tot, "); pls verify")
    rem <- 1 - tot
    weights <- v / tot
    # unlist(terra::global(weights * rem,'sum')) == rem
    v <- v + (weights * rem)
  }
  
  if(tot > 1) {
    if(tot - 1 > .005) warning("Sum of cell values not equal to 1 (", round(tot, 5), "); pls verify")
    rem <- tot - 1
    weights <- v / tot
    # unlist(terra::global(weights * rem,'sum')) == rem
    v <- v - (weights * rem)
  }
  
  if(tout == 'prob') {
    return(v / cell.area)
  } else {
    return(v)
  }
  
}
