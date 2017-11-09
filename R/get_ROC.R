get_ROC <- function(rankdiff_vector, labls) {
  ROC <- compute_threshold(rankdiff_vector, labls)
  return(ROC)
}
