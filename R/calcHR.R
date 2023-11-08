#' calcHR
#'
#' Function to calculate home range size and isopleth values corresponding to UD probability contours. Adapted from code shared by John Fieberg.
#'
#' @param est UD object, list containing 3 elements: x (x1) and y (x2) coordinates and corresponding density estimates (fhat)
#' @param p vector of probabilities at which to derive home range size
#' @param silent logical indicating whether or not to print messages (default = TRUE)
#'
#' @return probability density bin values associated with user-defined isopleths
#' @export
#'
calcHR <- function(est, p = c(0.5, 0.95), silent = TRUE) {

  ##  p = quantiles for estimation and display purposes

  x <- est$x1
  y <- est$x2
  z <- est$fhat

  # Calculate HR size
  #  Get size of grid
  dx <- unique(x)[2] - unique(x)[1]
  dy <- unique(y)[2] - unique(y)[1]

  # Estimate total volume and make sure = 1 or close to
  totp <- sum(z, na.rm=T) * dx * dy

  if(round(totp, 4) < 0.95) {
    warning("ERROR, total probability not equal to 1")
    message(paste("totp = ", totp))
  } else { 
    # to ensure integrates to 1
    z <- z / totp 
  }

  #  Get minimum number of squares to encompass pHR% probability
  nps <- length(p)
  zsort <- sort(as.vector(z))
  pt <- cumsum(zsort) / sum(zsort)

  #  p's for plotting and number of cells for calculating HR
  pplot <- matrix(0, nps, 1)
  ncells <- matrix(0, nps, 1)
  for(i in 1:nps) {
    if(silent==FALSE) { cat(p[i],"\n") }
    pplot[i] <- min(zsort[pt >= (1-p[i])])
    ncells[i] <- length(pt[pt >= (1-p[i])])
  }

  # Home range area
  HR <- dx * dy * ncells
  out <- list(HR, pplot, totp)
  names(out) <- c("HR", "ps", "totp")

  return(out)

}

