#' Compute a suite of common metrics serving to assess predictive capacity of an empirical/training model given an observed/validation sample
#'
#' @param r a spatRaster consisting of exactly 2 layers, 1) continuous probabilities estimated from training data; 2) continuous probabilities obtained from out-of-sample (validation) data
#' @param probs probability levels consisting of a 'successful outcome' (e.g. presence vs. absence) from which to compute validation metrics
#'
#' @return data.frame with tabulated results
#' @export
#'
validate <- function(r, probs = c(50,70)) {
  
  stopifnot(terra::nlyr(r) == 2)
  
  contab <- lapply(probs, function(p) {
    values(stats::setNames(sapp(r, fhat2confin) >= p, c('pred','true')), dataframe = TRUE)
  })
  
  data.frame(prob = probs,
             rmse = sapply(contab, function(x) mltools::rmse(pred = x[,'pred'], actuals = x[,'true'])),
             auc_roc = sapply(contab, function(x) mltools::auc_roc(preds = x[,'pred'], actuals = x[,'true'])),
             precision = sapply(contab, function(x) caret::precision(data = factor(x[,'pred']), reference = factor(x[,'true']), relevant = 'TRUE')),
             recall = sapply(contab, function(x) caret::recall(data = factor(x[,'pred']), reference = factor(x[,'true']), relevant = 'TRUE')),
             F05 = sapply(contab, function(x) caret::F_meas(data = factor(x[,'pred']), reference = factor(x[,'true']), relevant = 'TRUE', beta = 0.5)),
             F1 = sapply(contab, function(x) caret::F_meas(data = factor(x[,'pred']), reference = factor(x[,'true']), relevant = 'TRUE', beta = 1)),
             mcc = sapply(contab, function(x) mltools::mcc(preds = x[,'pred'], actuals = x[,'true'])))
  
}