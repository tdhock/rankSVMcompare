prediction_rankdiff <- function(Xirank,tau){

  #Gives outcome predictions based on given tau thresholds, taking a matrix of rank differences  

  for (i in 1:nrow(Xirank)) {
    diff <- Xirank[i, 2] - Xirank[i, 1]
    if (diff < -1 * tau) {
      Xirank[i, 4] <- -1
    }
    if (diff > tau) {
      Xirank[i, 4] <- 1
    }
    if (abs(diff) <= tau) {
      Xirank[i, 4] <- 0
    }
  }
  return(Xirank)
}