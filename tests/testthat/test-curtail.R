context('curtail function calls')

##public functions

test_that("criticalValues function works", {
  resp <- criticalValues(n=c(6, 30), p=c(.8, .2), pearly = .1, alpha = .1)
  expect_equal(resp[1], 4)
  expect_equal(resp[2], 11)
  expect_is(resp, 'vector')
})


test_that("probEarlyStop function works", {
  resp <- probEarlyStop(p=0.8, n = 6, r = 4)
  expect_equal(round(resp, 5), 0.09888)
  expect_is(resp, 'numeric')
})


test_that("expectedStage1SampleSize function works", {
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
  resp <- minimaxDesign(p = c(.8, .2), n = c(3,33), r = ph2crit(n=c(3,33), p=c(.8, .2), pearly = .1, alpha =.1))
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
  expect_equal(round(resp[1,], 6), c(2, 8, 7.214133))
  expect_equal(round(resp[2,], 6), c(6, 4, 7.234434))
  expect_equal(round(resp[3,], 6), c(5, 5, 7.272286))
  expect_equal(round(resp[4,], 6), c(4, 6, 7.329535))
  expect_equal(round(resp[5,], 6), c(3, 7, 7.383349))
  expect_equal(round(resp[6,], 6), c(7, 3, 7.426047))
  expect_equal(round(resp[7,], 6), c(8, 2, 7.587027))
  expect_equal(round(resp[8,], 6), c(9, 1, 7.986812))
  expect_equal(round(resp[9,], 6), c(10, 0, 8.150780))
  expect_is(resp, 'data.frame')
})

test_that("allMinimaxDesigns function works", {
  resp <- allMinimaxDesigns(p = c( .8, .2), n = 10)
  expect_equal(round(resp[1,], 8), c(6, 4, 0.06404506))
  expect_equal(round(resp[2,], 8), c(2, 8, 0.06491341))
  expect_equal(round(resp[3,], 8), c(5, 5, 0.06498202))
  expect_equal(round(resp[4,], 8), c(9, 1, 0.06502810))
  expect_equal(round(resp[5,], 8), c(8, 2, 0.06545818))
  expect_equal(round(resp[6,], 8), c(4, 6, 0.06559949))
  expect_equal(round(resp[7,], 8), c(7, 3, 0.06574490))
  expect_equal(round(resp[8,], 8), c(3, 7, 0.06593741))
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
  expect_is(resp, 'vector')
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
  expect_is(resp, 'vector')
})

test_that("is.wholenumber function works", {
  resp <- is.wholenumber(2.365)
  expect_equal(resp, FALSE)
  expect_is(resp, 'logical')
})