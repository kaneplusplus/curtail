context('ph2 function calls')

test_that("ph2crit function works", {
  resp <- ph2crit(n=c(6, 30), p=c(.8, .2), pearly = .1, alpha = .1)
  expect_equal(resp[1], 4)
  expect_equal(resp[2], 11)
  expect_is(resp, 'vector')
})


test_that("ph2early function works", {
  resp <- ph2early(p=0.8, n = 6, r = 4)
  expect_equal(resp, 0.09888)
  expect_is(resp, 'numeric')
})