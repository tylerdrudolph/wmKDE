pgoi <- function(pred, obs) {
  
  # gs <- as.vector(global(pred, 'sum'))[[1]]
  # pred <- pred / gs
  # 
  # gs <- as.vector(global(obs, 'sum'))[[1]]
  # obs <- obs / gs
  
  x <- values(pred)[,1]
  y <- values(obs)[,1]
  (2 - sum(pmax(x / sum(x), y / sum(y)))) / 1
  
  # (2 - global(app(c(pred, obs), fun='max'), 'sum')) / 1
  
}
