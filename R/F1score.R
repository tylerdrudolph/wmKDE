#' F1 Score for binned continuous data (kernel utilization distributions)
#'
#' @param pred spatRaster object corresponding to predicted probabilities of occurrence (wmKDE()$iso) from fitted wmKDE (training set) model (0 <= x >= 100).
#' @param obs spatRaster object corresponding to observed probabilities of occurrence (wmKDE()$iso) from fitted wmKDE (validation set) model (0 <= x >= 100). Alternatively, can be an sf object ('POINT' or 'MULTIPOINT') corresponding to 'true' observations.
#' @param id for weighted measures; character vector of length one corresponding to the field in obs (if SpatialPoints) distinguishing factor levels, in which case a confusion table will be constructed for each level of the factor and weighted according to the number of points in each level.
#' @param threshVal minimum probability value considered as a 'presence' (TRUE).
#' @param binWidth range of values per bin on a scale of 100. Default value of 100 will result in two bins: >= threshVal & > threshVal.
#' @param beta relative weight of recall relative to precision. Default is equal weight to both.
#' 
#' @importFrom stats setNames
#' @importFrom caret recall
#' @importFrom caret precision
#' @importFrom caret F_meas
#' 
#' @return either a list containing 1) recall, 2) precision, and 3) F1 score for the evaluated data (predicted vs observed) and input parameters, or a data.frame if obs is binary and !is.null(id).
#' @export
#' 
F1score <- function(pred, obs, id = NULL, threshVal = 70, binWidth = 100, beta = 1) {

  FScore <- NULL
  
  browser()
  
  cuts <- sort(c(seq(0, 100, binWidth), threshVal))
  cuts <- cuts[!duplicated(cuts)]

  Fcalc <- function(v, beta = beta) {
    list(recall = caret::recall(v, relevant = row.names(v)[(cuts > threshVal)[-1]]),
         precision = caret::precision(v, relevant = row.names(v)[(cuts > threshVal)[-1]]),
         FScore = caret::F_meas(v, relevant = row.names(v)[(cuts > threshVal)[-1]], beta = 1))
  }
  
  if(inherits(obs, 'SpatRaster')) {
    
    return(Fcalc(terra::crosstab(c(terra::classify(stats::setNames(pred, 'pred'),
                                                   rcl = cuts), 
                                   terra::classify(stats::setNames(obs, 'true'),
                                                   rcl = cuts)))))
  } else {
    
    if(is.null(id)) {
      
      return(Fcalc(terra::crosstab(c(terra::classify(stats::setNames(pred, 'pred'),
                                                     rcl = cuts), 
                                     terra::classify(stats::setNames(
                                       terra::rasterize(x = obs, y = pred, 
                                                        fun = function(x) 100, background = 0.1), 'true'), rcl = cuts)))))
      
    } else {
      
      Uid <- pull(obs, id) %>% unique
      nobs <- table(pull(obs, !!as.name(id)))
      
      v <- setNames(lapply(Uid, function(i) {
        terra::crosstab(c(terra::classify(stats::setNames(pred, 'pred'),
                                          rcl = cuts), 
                          terra::classify(stats::setNames(
                            terra::rasterize(x = dplyr::filter(obs, !!as.name(id) == i), 
                                             y = pred, 
                                             fun = function(x) 100, background = 0.1), i), rcl = cuts)))
      }), Uid)
      
      ftab <- do.call(rbind, lapply(v, function(k) t(Fcalc(k)))) %>%
        as.data.frame %>%
        dplyr::mutate(dplyr::across(recall:FScore, ~ unlist(.x)),
                      dplyr::across(recall:FScore, ~ ifelse(is.na(.x), 0, .x)),
                      id = names(v), .before = recall) %>% 
        dplyr::bind_cols(nobs = as.vector(nobs)) %>%
        dplyr::mutate(w = (nobs / sum(nobs))[order(nobs, decreasing = T)])
      
      return(dplyr::bind_rows(ftab,
                              data.frame(id = 'wtdMean',
                                         # rbind(
                                         # t(apply(ftab[,c(2:4)], 2, function(x) mean(x, na.rm = T))),
                                         t(apply(ftab[,c(2:4)], 2, function(x) stats::weighted.mean(x, w = ftab$w))), 
                                         # n = NA,
                                         w = sum(ftab$w))))
      
    }
    
  }
  
  # caret::confusionMatrix(v, mode = 'everything',
  #                        positive = row.names(v)[(cuts > threshVal)[-1]])
  
}

