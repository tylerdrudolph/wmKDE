#' F1 Score for binned continuous data (kernel utilization distributions)
#'
#' @param pred spatRaster object corresponding to predicted probabilities of occurrence (wmKDE()$iso) from fitted wmKDE (training set) model (0 <= x >= 100).
#' @param obs spatRaster object corresponding to observed probabilities of occurrence (wmKDE()$iso) from fitted wmKDE (validation set) model (0 <= x >= 100).
#' @param threshVal minimum probability value considered as a 'presence' (TRUE).
#' @param nbins number of bins to create for the confusion matrix (default = 10).
#' @param beta relative weight of recall relative to precision. Default is equal weight to both.
#' 
#' @return list containing 1) recall, 2) precision, and 3) F1 score for the evaluated data (predicted vs observed) and input parameters.
#' @export
#' 

F1score <- function(pred, obs, threshVal = 50, nbins = 10, beta = 1) {
  
  cuts <- seq(0, 100, nbins)
  
  v <- terra::crosstab(c(terra::classify(stats::setNames(pred, 'pred'),
                           rcl = cuts), 
                  terra::classify(stats::setNames(obs, 'true'),
                           rcl = cuts)))
  
  list(recall = caret::recall(v, relevant = row.names(v)[(cuts > threshVal)[-1]]),
       precision = caret::precision(v, relevant = row.names(v)[(cuts > threshVal)[-1]]),
       FScore = caret::F_meas(v, relevant = row.names(v)[(cuts > threshVal)[-1]], beta = 1))
  
  # confusionMatrix(v, mode = 'prec_recall',
  #                 positive = row.names(v)[(cuts > threshVal)[-1]])
  
}
