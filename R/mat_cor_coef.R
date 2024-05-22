mat_cor_coef <- function(pred, obs, threshVal = 70) {
  
  mltools::mcc(terra::as.matrix(wmKDE::fhat2confin(pred) >= threshVal)[,1],
               terra::as.matrix(wmKDE::fhat2confin(obs) >= threshVal)[,1])

}
