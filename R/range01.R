#' Rescale values of a vector to between 0 and 1
#'
#' @param x numeric vector
#'
#' @return
#' @export
#'
#' @examples
range01 = function(x) {
  (x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
}
