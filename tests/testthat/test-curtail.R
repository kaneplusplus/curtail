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
  resp <- ph2Eearly (p = 0.8, n = 6, r = 4)
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