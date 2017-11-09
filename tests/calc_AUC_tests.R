## TESTS calcAUC
test_rankdiff_vector <- c(0.1, 0.2, -0.3, 0.4)
test_labels <- c(1, 0, -1, 0)
AUC_test <- calc_AUC(test_rankdiff_vector, test_labels)

if (AUC_test != 0.25) {
  stop()
}

test_rankdiff_vector <- c(0,1)
test_labels <- c(0,1)
AUC_test <- calc_AUC(test_rankdiff_vector, test_labels)
if (AUC_test != 1) {
  stop()
}

# TESTS Baseline
test_labels <- c(1,0,-1,0,1,1)
AUC_test <- calc_Baseline(test_labels)
if(AUC_test != 0.375){
  stop()
}






