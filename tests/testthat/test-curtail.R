context('ph2 function calls')

test_that("ph2crit function works", {
  resp <- ph2crit(n=c(6, 30), p=c(.8, .2), pearly = .1, alpha = .1)
  expect_equal(resp[1], 4)
  expect_equal(resp[2], 11)
  expect_is(resp, 'vector')
})


test_that("ph2early function works", {
  resp <- ph2early(p=0.8, n = 6, r = 4)
  expect_equal(round(resp, 5), 0.09888)
  expect_is(resp, 'numeric')
})


test_that("ph2Eearly function works", {
  resp <- ph2Eearly(p = 0.8, n = 6, r = 4)
  expect_equal(resp[1], 4.76)
  expect_equal(round(resp[2], 7), 0.7797435)
 expect_is(resp, 'vector')
})


test_that("ph2Ess function works", {
  resp <- ph2Ess(p = c(0.8, 0.2), n= c(6, 30), r = c(4, 11))
  expect_equal(round(resp, 5), 29.32231)
  expect_is(resp, 'numeric')
})

### optimal and minimax functions

test_that("ph2designs function works", {
  resp <- ph2designs(p = c(0.8, 0.2), ntot = 36, pearly = .1, alpha = .1)
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

test_that("ph2mmax function works", {
  resp <- ph2mmax(p = c(.8, .2), n = c(3,33), r = ph2crit(n=c(3,33), p=c(.8, .2), pearly = .1, alpha =.1))
  expect_equal(round(resp, 5), 0.07063)
  expect_is(resp, 'numeric')
})

test_that("ph2reject function works", {
  resp <- ph2reject(p = c( .8, .2), n = c(12, 24), r = c(8, 11))
  expect_equal(round(resp, 5), 0.08578)
  expect_is(resp, 'numeric')
})

test_that("ph2rejcs function works", {
  resp <- ph2rejcs(p = c( .8, .2), n = c(12, 24), r = c(8, 11))
  expect_equal(round(resp, 5), 0.09044)
  expect_is(resp, 'numeric')
})

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

test_that("ph2rejcs function works", {
  resp <- ph2valid(p = c( .2, .3), n = c(10, 10), r = c(5, 5))
  expect_equal(resp, FALSE)
  expect_is(resp, 'logical')
})

