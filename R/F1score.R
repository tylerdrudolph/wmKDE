F1score <- function(ud1, ud2, threshVal = 50, nbins = 10, beta = 1) {
  
  cuts <- seq(0, 100, nbins)
  
  v <- crosstab(c(classify(setNames(ud1, 'pred'),
                           rcl = cuts), 
                  classify(setNames(ud2, 'true'),
                           rcl = cuts)))
  
  list(recall = recall(v, relevant = row.names(v)[(cuts > threshVal)[-1]]),
       precision = precision(v, relevant = row.names(v)[(cuts > threshVal)[-1]]),
       FScore = F_meas(v, relevant = row.names(v)[(cuts > threshVal)[-1]], beta = 1))
  
  confusionMatrix(v, mode = 'prec_recall',
                  positive = row.names(v)[(cuts > threshVal)[-1]])
  
}
