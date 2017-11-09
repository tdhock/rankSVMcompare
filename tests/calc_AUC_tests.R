## TESTS
test_rankdiff_vector <- c(0.1, 0.2, -0.3, 0.4)
test_labels <- c(1, 0, -1, 0)
if (0.25 != calc_AUC(test_rankdiff_vector, test_labels)) {
  stop()
}

test_rankdiff_vector <- c(0,1)
test_labels <- c(0,1)
if (1 != calc_AUC(test_rankdiff_vector, test_labels)) {
  stop()
}
