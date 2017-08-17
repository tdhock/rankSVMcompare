library(testthat)
context("soft margin")
library(rankSVMcompare)

data(separable)

no.var <- cbind(zero=rep(0, nrow(separable$Xi)),
                one=rep(1, nrow(separable$Xi)))
bad.pairs <- list(Xi=no.var, Xip=no.var, yi=separable$yi)

test_that("bad.pairs data set errors", {
  expect_error({
    softCompareQP(bad.pairs)
  }, "variance=0 for all input features")
})

good.bad <- with(separable, {
  list(Xi=cbind(Xi, no.var),
       Xip=cbind(Xip, no.var),
       yi=yi)
})

test_that("good.bad data set can be trained", {
  fit <- softCompareQP(good.bad)
  expect_is(fit, "list")
})

## Test data that does not have all the required columns.
not.enough <- with(separable, {
  list(Xi=cbind(Xi[,1,drop=FALSE], no.var),
       Xip=cbind(Xip[,1,drop=FALSE], no.var),
       yi=yi)
})
## Test data that does not have named columns.
not.named <- good.bad
for(i in seq_along(not.named)){
  colnames(not.named[[i]]) <- NULL
}

test_that("not.named data errors", {
  expect_error({
    softCompareQP(not.named)
  }, "input feature matrix columns must have names")
})


fit <- softCompareQP(separable)
test_that("separable data set can be trained", {
  expect_is(fit, "list")
})

test_that("prediction works for good.bad data", {
  pred.vec <- with(good.bad, fit$predict(Xi, Xip))
  expect_true(is.numeric(pred.vec))
  expect_equal(length(pred.vec), nrow(good.bad$Xi))
  expect_true(all(pred.vec %in% c(-1, 0, 1)))
})

test_that("prediction works for separable data", {
  pred.vec <- with(separable, fit$predict(Xi, Xip))
  expect_true(is.numeric(pred.vec))
  expect_equal(length(pred.vec), nrow(separable$Xi))
  expect_true(all(pred.vec %in% c(-1, 0, 1)))
})

test_that("prediction errors for bad.pairs", {
  expect_error({
    with(bad.pairs, fit$predict(Xi, Xip))
  }, "missing features distance, angle")
})

test_that("prediction errors for not.enough", {
  expect_error({
    with(not.enough, fit$predict(Xi, Xip))
  }, "missing features angle")
})

test_that("prediction errors for not.named", {
  expect_error({
    with(not.named, fit$predict(Xi, Xip))
  }, "missing features distance, angle")
})

