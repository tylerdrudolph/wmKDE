#' Convert UD density estimates into relative probabilities
#'
#' @param z numeric vector representing UD density values (fhat)
#'
#' @return
#' @export
#'
#' @examples
fhat2confin <- function(z) {
  dimxy = dim(z)
  if(sum(z)!=1) z = z/sum(z)
  index  = 1:length(z)
  vord   = z[order(z, decreasing = TRUE)]
  indord = index[order(z, decreasing = TRUE)]
  vsu    = cumsum(vord)
  vreord = vsu[order(indord)] * 100
  return(matrix(vreord,ncol=dimxy[2]))
}
