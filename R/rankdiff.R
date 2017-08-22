rankdiff <- function(qp, matrix_Xi, matrix_Xip) {

  #Computes the rankdiff given rsvm objects  

  Xirank <- qp$rank.scaled(X = matrix_Xi)
  Xiprank <- qp$rank.scaled(X = matrix_Xip)
  
  Xirank <- cbind(Xirank, Xiprank)
  Xirank <- data.frame(Xirank)
  
  Xirank$difference <- abs(Xirank$X2 - Xirank$X1)
  return(Xirank)
}
