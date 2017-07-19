library(testthat)

context('curtail function calls')

##public functions

test_that("dsnb_stacked function works", {
  s <- 5
  t <- 7
  resp <- dsnb_stacked(x = min(s,t):(s+t-1), p = .2, s = s, t = t)
  expect_equal(resp[, "x"], 5:11)
  expect_equal(round(resp[, "s"], 6), c(0.00032, 0.00128, 0.003072, 0.005734, 0.009175, 0.013212, 0.017616))
  expect_equal(round(resp[, "t"], 6), c(0, 0, 0.209715, 0.293601, 0.234881, 0.140929, 0.070464))
  expect_is(resp, "matrix")
})
test_that("cdsnb_stacked function works", {
  s <- 5
  t <- 7
  resp <- cdsnb_stacked(x=min(s,t):(s+t-1), shape=c(1, 1), s=s, t=t)
  expect_equal(resp[, "k"], 5:11)
  expect_equal(round(resp[, "s"], 6), c(0.166667, 0.119048, 0.089286, 0.069444, 0.055556, 0.045455, 0.037879))
  expect_equal(round(resp[, "t"], 6), c(0, 0, 0.125, 0.097222, 0.077778, 0.063636, 0.053030))
  expect_is(resp, "data.frame")
})

test_that("esnb function works", {
  resp <- esnb(p = .4, s = 7, t=11)
  expect_equal(round(resp, 5), 14.50153)
  expect_is(resp, "numeric")
})

test_that("vsnb function works", {
  resp <- vsnb(p = .4, s = 7, t=11)
  expect_equal(round(resp, 6), 4.544609)
  expect_is(resp, "numeric")
})

test_that("ecsnb function works", {
  resp <- ecsnb(shape=c(1, 1), s=7, t=11)
  expect_equal(round(resp, 5), 11.54329)
  expect_is(resp, "numeric")
})

test_that("vcsnb function works", {
  resp <- vcsnb(shape=c(1, 1), s=7, t=11)
  expect_equal(round(resp, 6), 9.209055)
  expect_is(resp, "numeric")
})

test_that("criticalValues function works for two-stage design", {
  resp <- criticalValues(n=c(6, 30), p=c(.8, .2), pearly = .1, alpha = .1)
  expect_equal(resp[1], 4)
  expect_equal(resp[2], 11)
  expect_is(resp, 'numeric')
})

test_that("criticalValues function works for one-stage design", {
  resp <- criticalValues(n=36, p=.2, alpha = .1)
  expect_equal(resp,11)
  expect_is(resp, 'numeric')
})

test_that("probEarlyStop function works", {
  resp <- probEarlyStop(p=0.8, n = 6, r = 4)
  expect_equal(round(resp, 5), 0.09888)
  expect_is(resp, 'numeric')
})

test_that("expectedStage1SampleSize function works for two-stage design", {
  resp <- expectedStage1SampleSize(p = c(0.8, .2), n = c(6, 30), r = c(4, 11))
  expect_equal(resp$expectation, 4.76)
  expect_equal(round(resp$standardDeviation, 7), 0.7797435)
 expect_is(resp, 'list')
})

test_that("expectedStage1SampleSize function works for one-stage design", {
  resp <- expectedStage1SampleSize(p = 0.8, n = 6, r = 4)
  expect_equal(resp$expectation, 4.76)
  expect_equal(round(resp$standardDeviation, 7), 0.7797435)
  expect_is(resp, 'list')
})


test_that("expectedTotalSampleSize function works", {
  resp <- expectedTotalSampleSize(p = c(0.8, 0.2), n= c(6, 30), r = c(4, 11))
  expect_equal(round(resp, 5), 29.32231)
  expect_is(resp, 'numeric')
})

### optimal and minimax functions

test_that("bestDesigns function works", {
  resp <- bestDesigns(p = c(0.8, 0.2), ntot = 36, pearly = .1, alpha = .1)
  expect_equal(resp["Optimal", "p1"], 0.8)
  expect_equal(resp["Optimal", "n1"], 6)
  expect_equal(resp["Optimal", "r1"], 4)
  expect_equal(resp["Optimal", "p2"], 0.2)
  expect_equal(resp["Optimal", "n2"], 30)
  expect_equal(resp["Optimal", "r2"], 11)
  expect_equal(resp["Optimal", "Alpha"], 0.1)
  expect_equal(round(resp["Optimal", "PET"], 5), 0.09888)
  expect_equal(round(resp["Optimal", "ECSS"], 5), 29.32231)
  expect_equal(round(resp["Optimal", "P(MaxSS)"], 8), 0.06584612)

  expect_equal(resp["Minimax", "p1"], 0.8)
  expect_equal(resp["Minimax", "n1"], 6)
  expect_equal(resp["Minimax", "r1"], 4)
  expect_equal(resp["Minimax", "p2"], 0.2)
  expect_equal(resp["Minimax", "n2"], 30)
  expect_equal(resp["Minimax", "r2"], 11)
  expect_equal(resp["Minimax", "Alpha"], 0.1)
  expect_equal(round(resp["Minimax", "PET"], 5), 0.09888)
  expect_equal(round(resp["Minimax", "ECSS"], 5), 29.32231)
  expect_equal(round(resp["Minimax", "P(MaxSS)"], 8), 0.06584612)
  
  expect_is(resp, 'data.frame')
})

test_that("minimaxDesign function works", {
  resp <- minimaxDesign(p = c(.8, .2), n = c(3,33), r = criticalValues(n=c(3,33), p=c(.8, .2), pearly = .1, alpha =.1))
  expect_equal(round(resp, 5), 0.07063)
  expect_is(resp, 'numeric')
})

test_that("probRejectTraditional function works", {
  resp <- probRejectTraditional(p = c( .8, .2), n = c(12, 24), r = c(8, 11))
  expect_equal(round(resp, 5), 0.08578)
  expect_is(resp, 'numeric')
})

test_that("probReject function works", {
  resp <- probReject(p = c( .8, .2), n = c(12, 24), r = c(8, 11))
  expect_equal(round(resp, 5), 0.09044)
  expect_is(resp, 'numeric')
})

test_that("allOptimalDesigns function works", {
  resp <- allOptimalDesigns(p = c( .8, .2), n = 10)
  expect_equal(resp$n1, c(2, 6, 5, 4, 3, 7, 8, 9))
  expect_equal(resp$n2, c(8, 4, 5, 6, 7, 3, 2, 1))
  expect_equal(round(resp$`Expected Sample Size`, 6), c(7.214133, 7.234434, 7.272286, 7.329535, 7.383349, 7.426047, 7.587027, 7.986812))
  expect_is(resp, 'data.frame')
})

test_that("allMinimaxDesigns function works", {
  resp <- allMinimaxDesigns(p = c( .8, .2), n = 10)
  expect_equal(resp$n1, c(6, 2, 5, 9, 8, 4, 7, 3))
  expect_equal(resp$n2, c(4, 8, 5, 1, 2, 6, 3, 7))
  expect_equal(round(resp$`Probability of Maximum Sample Size`, 8), c(0.06404506, 0.06491341, 0.06498202, 0.06502810, 0.06545818, 0.06559949, 0.06574490, 0.06593741))
  expect_is(resp, 'data.frame')
})


#### private functions

test_that("ph2snb function works", {
  resp <- ph2snb(p = .8, s = 3, t = 5)
  expect_equal(round(resp[1], 5), 0.00000)
  expect_equal(round(resp[2], 5), 0.00000)
  expect_equal(round(resp[3], 5), 0.51200)
  expect_equal(round(resp[4], 5), 0.30720)
  expect_equal(round(resp[5], 5), 0.12320)
  expect_equal(round(resp[6], 5), 0.04224)
  expect_equal(round(resp[7], 5), 0.01536)
  expect_is(resp, 'numeric')
})

test_that("ph2valid function works", {
  resp <- ph2valid(p = c( .2, .3), n = c(10, 10), r = c(5, 5))
  expect_equal(resp, FALSE)
  expect_is(resp, 'logical')
})

test_that("ph2tnb function works", {
  resp <- ph2tnb(p = c( .8, .2), n = c(6, 30), r = c(4, 11))
  expect_equal(resp[1], 0)
  expect_equal(resp[2], 0)
  expect_equal(resp[3], 0)
  expect_equal(round(resp[4], 7), 0.4545455)
  expect_equal(round(resp[5], 7), 0.3636364)
  expect_equal(round(resp[6], 7), 0.1818182)
  expect_is(resp, 'numeric')
})

test_that("is.wholenumber function works", {
  resp <- is.wholenumber(2.365)
  expect_equal(resp, FALSE)
  expect_is(resp, 'logical')
})
