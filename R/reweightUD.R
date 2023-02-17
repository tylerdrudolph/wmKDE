#' Redistribute UD probability densities proportionally so as to ensure that the sum of cell volumes equals one.
#'
#' @param UD a list containing 3 elements: 'x1' (x coordinates), 'x2' (y coordinates), 'fhat' (UD probability densities).
#' @param out integer parameter specifying whether output cell values are aspatial (out=1) or volume-based (out=2).
#'
#' @return original input UD with probabilities adjusted such that corresponding volumes sum to exactly 1
#' @export
#'
#' @examples
#'
reweightUD <- function(UD, out=1) {
  cell.area = (sort(unique(UD$x1))[2] - sort(unique(UD$x1))[1]) * (sort(unique(UD$x2))[2] - sort(unique(UD$x2))[1])# (UD$x1[2]-UD$x1[1]) * (UD$x2[2]-UD$x2[1])
  mat = UD$fhat
  v = mat * cell.area
  if(sum(v) < 1) {
    tot = sum(v)
    if(1 - tot > .005) warning("Sum of cell values not equal to 1 (", tot, "); pls verify")
    rem = 1 - sum(v)
    weights = v / sum(v)
    sum(weights * rem) == rem
    v = v + (weights * rem)
  }
  if(sum(v) > 1) {
    tot = sum(v)
    if(tot -1 > .005) warning("Sum of cell values not equal to 1 (", round(tot, 5), "); pls verify")
    rem = sum(v) - 1
    weights = v / sum(v)
    sum(weights*rem) == rem
    v = v - (weights * rem)
  }
  if(out==1) {
    UD$fhat = v / cell.area
  } else {
      UD$fhat = v
  }
  return(UD)
}


