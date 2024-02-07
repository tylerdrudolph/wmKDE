#' Title
#'
#' @param r 
#' @param out 
#'
#' @return
#' @export
#'
#' @examples
finetune <- function(r, out = 1) {
  
  cell.area <- prod(res(r))
  v <- r * cell.area
  tot <- unlist(global(v, 'sum'))
  
  if(tot < 1) {
    if(1 - tot > .005) warning("Sum of cell values not equal to 1 (", tot, "); pls verify")
    rem <- 1 - tot
    weights <- v / tot
    # unlist(global(weights * rem,'sum')) == rem
    v <- v + (weights * rem)
  }
  
  if(tot > 1) {
    if(tot - 1 > .005) warning("Sum of cell values not equal to 1 (", round(tot, 5), "); pls verify")
    rem <- tot - 1
    weights <- v / tot
    # unlist(global(weights * rem,'sum')) == rem
    v <- v - (weights * rem)
  }
  
  if(out == 1) {
    return(v / cell.area)
  } else {
    return(v)
  }
  
}
