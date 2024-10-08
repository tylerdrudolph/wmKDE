############################################################################################
## Calculate overlap indices between different UDs in a list

olapCalc <- function (ud_list, method = c("HR", "PHR", "VI", "BA", "UDOI", "HD"),
                      percent = 95, conditional = TRUE, ...) {
  if(length(ud_list) == 1)
    stop("multiple UDs are needed for this function")
  if(round(sum(ud_list[[1]]$fhat * spatres^2)) != 1)
    stop("x should not be a volume under UD")
  ## get volumes in object 'vol'
  vol <- lapply(ud_list, function(x) {
    x$fhat <- fhat2confin(x$fhat)
    return(x)
  })
  ## sort coordinates of non-volume UD
  ud_list <- lapply(ud_list, function(ud) {
    coo <- cbind(ud$x1, ud$x2)
    as.vector(ud$fhat[order(coo[,1], coo[,2]), ])
  })
  ## sort coordinates of volume UD
  vol <- lapply(vol, function(ud) {
    coo <- cbind(ud$x1, ud$x2)
    as.vector(ud$fhat[order(coo[,1], coo[,2]), ])
  })
  
  ## output matrix to stock
  res <- matrix(0, ncol = length(ud_list), nrow = length(ud_list))
  
  for (i in 1:length(ud_list)) {
    for (j in 1:i) {
      if (method == "HR") {
        vi <- vol[[i]]
        vj <- vol[[j]]
        vi[vi <= percent] <- 1
        vi[vi > percent] <- 0
        vj[vj <= percent] <- 1
        vj[vj > percent] <- 0
        vk <- vi * vj
        res[i, j] <- sum(vk)/sum(vi)
        res[j, i] <- sum(vk)/sum(vj)
      }
      if (method == "PHR") {
        vi <- ud_list[[i]]
        vj <- ud_list[[j]]
        ai <- vol[[i]]
        aj <- vol[[j]]
        ai[ai <= percent] <- 1
        ai[ai > percent] <- 0
        aj[aj <= percent] <- 1
        aj[aj > percent] <- 0
        if (conditional) {
          vi <- vi * ai
          vj <- vj * aj
          res[j, i] <- sum(vi * aj) * (spatres^2)
          res[i, j] <- sum(vj * ai) * (spatres^2)
        }
        else {
          res[j, i] <- sum(vi * aj) * (spatres^2)
          res[i, j] <- sum(vj * ai) * (spatres^2)
        }
      }
      if (method == "VI") {
        vi <- ud_list[[i]]
        vj <- ud_list[[j]]
        ai <- vol[[i]]
        aj <- vol[[j]]
        ai[ai <= percent] <- 1
        ai[ai > percent] <- 0
        aj[aj <= percent] <- 1
        aj[aj > percent] <- 0
        if (conditional) {
          vi <- vi * ai
          vj <- vj * aj
          res[i, j] <- res[j, i] <- sum(pmin(vi, vj)) *
            (spatres^2)
        }
        else {
          res[i, j] <- res[j, i] <- sum(pmin(vi, vj)) *
            (spatres^2)
        }
      }
      if (method == "BA") {
        vi <- ud_list[[i]]
        vj <- ud_list[[j]]
        ai <- vol[[i]]
        aj <- vol[[j]]
        ai[ai <= percent] <- 1
        ai[ai > percent] <- 0
        aj[aj <= percent] <- 1
        aj[aj > percent] <- 0
        if (conditional) {
          vi <- vi * ai
          vj <- vj * aj
          res[j, i] <- res[i, j] <- sum(sqrt(vi) * sqrt(vj)) *
            (spatres^2)
        }
        else {
          res[j, i] <- res[i, j] <- sum(sqrt(vi) * sqrt(vj)) *
            (spatres^2)
        }
      }
      if (method == "UDOI") {
        vi <- ud_list[[i]]
        vj <- ud_list[[j]]
        ai <- vol[[i]]
        aj <- vol[[j]]
        ai[ai <= percent] <- 1
        ai[ai > percent] <- 0
        aj[aj <= percent] <- 1
        aj[aj > percent] <- 0
        if (conditional) {
          vi <- vi * ai
          vj <- vj * aj
          ak <- sum(ai * aj) * (spatres^2)
          res[j, i] <- res[i, j] <- ak * sum(vi * vj) *
            (spatres^2)
        }
        else {
          ak <- sum(ai * aj) * (spatres^2)
          res[j, i] <- res[i, j] <- ak * sum(vi * vj) *
            (spatres^2)
        }
      }
      if (method == "HD") {
        vi <- ud_list[[i]]
        vj <- ud_list[[j]]
        ai <- vol[[i]]
        aj <- vol[[j]]
        ai[ai <= percent] <- 1
        ai[ai > percent] <- 0
        aj[aj <= percent] <- 1
        aj[aj > percent] <- 0
        if (conditional) {
          vi <- vi * ai
          vj <- vj * aj
          res[j, i] <- res[i, j] <- sqrt(sum((sqrt(vi) -
                                                sqrt(vj))^2 * (spatres^2)))
        }
        else {
          res[j, i] <- res[i, j] <- sqrt(sum((sqrt(vi) -
                                                sqrt(vj))^2 * (spatres^2)))
        }
      }
    }
  }
  rownames(res) <- names(ud_list)
  colnames(res) <- names(ud_list)
  return(res)
}

