rankdiff <- function(rank, rank_p, t) {

  Xirank <- cbind(rank, rank_p)
  
  Xirank <- data.frame(Xirank)
  Xirank$diff <- Xirank[, 2] - Xirank[, 1]
  Xirank$label <- NA
  
  for (i in 1:nrow(Xirank)) {
    diff <- Xirank[i, 2] - Xirank[i, 1]
    if (diff < (-1 * t)) {
      Xirank[i, 4] <- -1
    }
    if (diff > t) {
      Xirank[i, 4] <- 1
    }
    if (abs(diff) <= t) {
      Xirank[i, 4] <- 0
    }
  }
  return(Xirank)
}
