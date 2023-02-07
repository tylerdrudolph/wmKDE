#' Derive the core area isopleth from a utilization distribution using the method described in Vander Wal & Rogers (2012).
#'
#' @param UD list containing 3 elements: 'x1' (x coordinates), 'x2' (y coordinates) and 'fhat' (probability densities)
#' @param silent logical indicating whether or not to print messages
#' @param plot logical indicating whether or not to print result
#'
#' @return
#' @export
#'
#' @examples
core.area <- function(UD, plot=FALSE, silent=TRUE) {

  ## Generate sequence of isopleth volumes
  IVseq <- seq(0.01, 0.99, 0.01)

  ## Estimate contour levels and size within each isopleth
  hrsum <- calcHR(UD, p = IVseq, silent=silent)
  kernel.areas <- c(hrsum$HR * 1e-6)

  ## Create a dataframe with the percent area (PA) and isopleth volume (IV) scaled from 0-1
  df <- as.data.frame(cbind(IVseq, kernel.areas / max(kernel.areas)))

  ## Name the columns
  colnames(df) <- c("IV","PA")

  ## CAUTION: Starting parameters may differ for your data.
  nls.fit <- stats::nls(PA~(b0) * (exp(b1^IV)), data=df, start=list(b0=0.01, b1=4.2), na.action="na.omit", model=TRUE)

  ## b0 coefficient
  b0 <- summary(nls.fit)$coefficients[1,1]

  ## b1 coefficient
  b1 <- summary(nls.fit)$coefficients[2,1]

  ## Isopleth volume where the curve's slope = 1
  cutval <- (-log(b0*b1)/b1)

  ## Plot
  if(plot) {
    graphics::plot(y=df[,2], x=df[,1], type="l", xlab="Confidence Interval", ylab="Probability")
    graphics::abline(v=cutval, lty="dotted", col="red")
  }

  return(cutval)
}
