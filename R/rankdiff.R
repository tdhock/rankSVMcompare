rankdiff <- function(qp, matrix_Xi, matrix_Xip, t) {
  Xirank <- qp$rank(X = matrix_Xi)
  Xiprank <- qp$rank(X = matrix_Xip)
  Xirank <- cbind(Xirank, Xiprank)
  
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
