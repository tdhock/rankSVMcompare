library(testthat)
context("convert pairs")
library(rankSVMcompare)

p <- list(Xi=cbind(var=c(3,0,1)),
          Xip=cbind(var=c(0,-2,0)),
          yi=as.integer(c(-1,1,0)))
result <- pairs2svmData(p)
## Inequality pairs such that yi=1 or -1 are mapped to 1, and
## equality pairs such that yi=0 are duplicated and mapped to -1.
test_that("inequality => 1, equality => -1", {
  expect_equal(result$yi, c(1,1,-1,-1))
})
test_that("The duplicate equality features are flipped", {
  expect_equal(result$Xi[3], result$Xip[4])
  expect_equal(result$Xi[4], result$Xip[3])
})

colnames(p$Xi) <- NULL

test_that("error when no column names", {
  expect_error({
    pairs2svmData(p)
  }, "input feature matrix columns must have names")
})
