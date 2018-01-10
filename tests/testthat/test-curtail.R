library(testthat)

context('curtail function calls')

#' Single stage designs
test_that("single_stage_curtail_trial works when s and t are parameters", {
p_null <- 0.2
p_alt <- 0.4
s <- 7
t <- 11
trial <- single_stage_curtail_trial(p_null=p_null, p_alt=p_alt, s=s, t=t)

expect_equal(as.numeric(trial[1:4]), c(p_null, p_alt, s, t))
expect_equal(trial$power, single_stage_power(pAlt = p_alt, s = s, t = t))
expect_equal(trial$significance, single_stage_significance(pNull = p_null, s=s, t=t))
expect_equal(trial$mean_ss, single_stage_expected_sample_size(p = p_null, s=s, t=t)$expectation)

expect_is(trial, "data.frame")
})

#trials <- single_stage_curtail_trial(p_null=0.2, p_alt=0.4, n=17)
test_that("single_stage_curtail_trial works when n is a parameter", {
  p_null <- 0.2
  p_alt <- 0.4
  n <- 17
  trials <- single_stage_curtail_trial(p_null=p_null, p_alt=p_alt, n= n)
  s <- seq_len(n-1)
  t <- n - s + 1
  
  df <- data.frame(p_null=rep(p_null, length(s)),
                   p_alt=rep(p_alt, length(s)), s=s, t=t)
  
  
  expect_equal(data.frame(trials[1:4]), df)
  expect_equal(trials$power, apply(df, MARGIN = 1, function(x) single_stage_power(x["p_alt"], x["s"], x["t"])))
  expect_equal(trials$significance, apply(df, MARGIN=1, function(x) single_stage_significance(x["p_null"], x["s"], x["t"])))
  expect_equal(trials$mean_ss, apply(df, MARGIN=1, function(x) single_stage_expected_sample_size(x["p_null"], x["s"], x["t"])$expectation))
  
  expect_is(trials, "data.frame")
})


#' Two-stage curtailed designs
#' Case 1 - User inputs p, n, and r
#' trial <- two_stage_curtail_trial(p1_null = 0.8, p2_null = 0.2, 
#' p1_alt = 0.8, p2_alt = 0.4, n1 = 6, n2 = 30, r1 = 4, r2 = 11)
test_that("two_stage_curtail_trial works for case 1", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  r1 <- 4
  r2 <- 11
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
     p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, r1=r1, r2=r2)
  
  expect_equal(as.numeric(trial[1:8]), c(p1_null, p2_null, p1_alt, p2_alt, n1, n2, r1, r2))
  expect_equal(trial$power, two_stage_power(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$significance, two_stage_significance(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$stage1_mean_ss, expected_stage1_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2))$expectation)
  expect_equal(trial$mean_ss_null, expected_total_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$mean_ss_alt, expected_total_sample_size(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$PET, prob_early_stop(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$minimax_prob, minimax_design_old(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  
  expect_is(trial, "data.frame")
})

#' Test significance function
test_that("significance function works for two-stage curtail trial", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  r1 <- 4
  r2 <- 11
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, r1=r1, r2=r2)
  sig <- significance(trial)
  expect_equal(sig, trial$significance)
  expect_is(sig, "numeric")
})

#' Test significance function
test_that("power function works for two-stage curtail trial", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  r1 <- 4
  r2 <- 11
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, r1=r1, r2=r2)
  pow <- power(trial)
  expect_equal(pow, trial$power)
  expect_is(pow, "numeric")
})

#' Test stage1_sample_size function
test_that("stage1_sample_size function works for two-stage curtail trial", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  r1 <- 4
  r2 <- 11
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, r1=r1, r2=r2)
  stage1 <- stage1_sample_size(trial)
  expect_equal(stage1, trial$stage1_mean_ss)
  expect_is(stage1, "numeric")
})


#' Test expected_sample_size function
test_that("expected_sample_size function works for two-stage curtail trial", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  r1 <- 4
  r2 <- 11
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, r1=r1, r2=r2)
  ess <- expected_sample_size(trial)
  expect_equal(ess, trial$mean_ss_null)
  expect_is(ess, "numeric")
})


#' Test expected_sample_size_alt function
test_that("expected_sample_size_alt function works for two-stage curtail trial", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  r1 <- 4
  r2 <- 11
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, r1=r1, r2=r2)
  essA <- expected_sample_size_alt(trial)
  expect_equal(essA, trial$mean_ss_alt)
  expect_is(ess, "numeric")
})

#' Test PET function
test_that("PET function works for two-stage curtail trial", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  r1 <- 4
  r2 <- 11
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, r1=r1, r2=r2)
  pet <- PET(trial)
  expect_equal(pet, trial$PET)
  expect_is(pet, "numeric")
})

#' Test minimax_probability function
test_that("minimax_probability function works for two-stage curtail trial", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  r1 <- 4
  r2 <- 11
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, r1=r1, r2=r2)
  minimax <- minimax_probability(trial)
  expect_equal(minimax, trial$minimax_prob)
  expect_is(minimax, "numeric")
})



#' Create a two_stage_curtail_trial
#' Case 2:  User inputs p, n
#' Using default values of alpha and prob_early
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30))
test_that("two_stage_curtail_trial works for case 2", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30

  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2)
  r1 <- two_stage_critical_values(n=c(n1, n2), p=c(p1_null, p2_null))[1]
  r2 <- two_stage_critical_values(n=c(n1, n2), p=c(p1_null, p2_null))[2]
  
  expect_equal(as.numeric(trial[1:8]), c(p1_null, p2_null, p1_alt, p2_alt, n1, n2, r1, r2))
  expect_equal(trial$power, two_stage_power(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$significance, two_stage_significance(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$stage1_mean_ss, expected_stage1_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2))$expectation)
  expect_equal(trial$mean_ss_null, expected_total_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$mean_ss_alt, expected_total_sample_size(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$PET, prob_early_stop(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$minimax_prob, minimax_design_old(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  
  expect_is(trial, "data.frame")
})

#' Create a two_stage_curtail_trial
#' Case 3:  User inputs p, n, prob_early, alpha
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30), prob_early=0.1, alpha=0.1)
test_that("two_stage_curtail_trial works for case 3", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  prob_early <- 0.1
  alpha <- 0.1

  
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, 
                                   prob_early=prob_early, alpha=alpha)
  r1 <- two_stage_critical_values(n=c(n1, n2), p=c(p1_null, p2_null), pearly=prob_early, alpha=alpha)[1]
  r2 <- two_stage_critical_values(n=c(n1, n2), p=c(p1_null, p2_null), pearly=prob_early, alpha=alpha)[2]
  
  expect_equal(as.numeric(trial[1:8]), c(p1_null, p2_null, p1_alt, p2_alt, n1, n2, r1, r2))
  expect_equal(trial$power, two_stage_power(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$significance, two_stage_significance(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$stage1_mean_ss, expected_stage1_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2))$expectation)
  expect_equal(trial$mean_ss_null, expected_total_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$mean_ss_alt, expected_total_sample_size(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$PET, prob_early_stop(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$minimax_prob, minimax_design_old(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  
  expect_is(trial, "data.frame")
})

#' Create a two_stage_curtail_trial
#' Case 3a:  User inputs p, n, prob_early
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30), prob_early=0.1)
test_that("two_stage_curtail_trial works for case 3a", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  prob_early <- 0.1

  
  
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, 
                                   prob_early=prob_early)
  r1 <- two_stage_critical_values(n=c(n1, n2), p=c(p1_null, p2_null), pearly=prob_early)[1]
  r2 <- two_stage_critical_values(n=c(n1, n2), p=c(p1_null, p2_null), pearly=prob_early)[2]
  
  expect_equal(as.numeric(trial[1:8]), c(p1_null, p2_null, p1_alt, p2_alt, n1, n2, r1, r2))
  expect_equal(trial$power, two_stage_power(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$significance, two_stage_significance(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$stage1_mean_ss, expected_stage1_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2))$expectation)
  expect_equal(trial$mean_ss_null, expected_total_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$mean_ss_alt, expected_total_sample_size(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$PET, prob_early_stop(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$minimax_prob, minimax_design_old(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  
  expect_is(trial, "data.frame")
})

#' Create a two_stage_curtail_trial
#' Case 3b:  User inputs p, n, alpha
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30), alpha=0.1)
test_that("two_stage_curtail_trial works for case 3b", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n1 <- 6
  n2 <- 30
  alpha <- 0.05
  
  
  
  trial <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n1=n1, n2=n2, 
                                   alpha=alpha)
  r1 <- two_stage_critical_values(n=c(n1, n2), p=c(p1_null, p2_null), alpha=alpha)[1]
  r2 <- two_stage_critical_values(n=c(n1, n2), p=c(p1_null, p2_null), alpha=alpha)[2]
  
  expect_equal(as.numeric(trial[1:8]), c(p1_null, p2_null, p1_alt, p2_alt, n1, n2, r1, r2))
  expect_equal(trial$power, two_stage_power(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$significance, two_stage_significance(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$stage1_mean_ss, expected_stage1_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2))$expectation)
  expect_equal(trial$mean_ss_null, expected_total_sample_size(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$mean_ss_alt, expected_total_sample_size(p = c(p1_alt, p2_alt), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$PET, prob_early_stop(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  expect_equal(trial$minimax_prob, minimax_design_old(p = c(p1_null, p2_null), n = c(n1, n2), r = c(r1, r2)))
  
  expect_is(trial, "data.frame")
})

#' Create a two_stage_curtail_trial
#' Case 4:  User inputs p, n_total, prob_early, alpha
#' trials <- two_stage_curtail_trial(p1_null = 0.8, p2_null=0.2, 
#' p1_alt = 0.8, p2_alt = 0.4, n_total=36, prob_early=0.1, alpha=0.1)
test_that("two_stage_curtail_trial works for case 4", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n_total <- 36
  prob_early=0.2
  alpha <- 0.05
  

  trials <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n_total = n_total, 
                                   prob_early=prob_early, alpha=alpha)

  
  n1 <- seq_len(n_total-1)
  n2 <- n_total - n1
  n <- cbind(n1, n2)
  r <- matrix(apply(n, 1, function(x){
    two_stage_critical_values(x, c(p1_null, p2_null),
                              pearly = prob_early,alpha = alpha)
  }), nrow=length(n1), ncol=2, byrow=TRUE)
  
  r1 <- r[,1]
  r2 <- r[,2]
  df <- data.frame(p1_null=rep(p1_null, length(n1)),
                    p2_null=rep(p2_null, length(n1)),
                    p1_alt=rep(p1_alt, length(n1)),
                    p2_alt=rep(p2_alt, length(n1)),
                    n1=n1, n2=n2, r1=r[,1], r2=r[,2])
  if(min(df$r1)==0){
    df <- subset(df, r1>0)
    rownames(df) <- NULL
  }
  
  expect_equal(data.frame(trials[,1:8]), df)
  
  expect_equal(trials$power, apply(df, MARGIN = 1, function(x) two_stage_power(p = c(x["p1_alt"], x["p2_alt"]), 
                                                    n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$significance, apply(df, MARGIN = 1, function(x) two_stage_significance(p = c(x["p1_null"], x["p2_null"]), 
                                                    n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$stage1_mean_ss, apply(df, MARGIN = 1, function(x) expected_stage1_sample_size(p = c(x["p1_null"], x["p2_null"]), 
                                                                                        n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))$expectation))
  expect_equal(trials$mean_ss_null, apply(df, MARGIN = 1, function(x) expected_total_sample_size(p = c(x["p1_null"], x["p2_null"]), 
                                                                                                 n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$mean_ss_alt, apply(df, MARGIN = 1, function(x) expected_total_sample_size(p = c(x["p1_alt"], x["p2_alt"]), 
                                                                                                n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$PET, apply(df, MARGIN = 1, function(x) prob_early_stop(p = c(x["p1_null"], x["p2_null"]), 
                                                                             n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$minimax_prob, apply(df, MARGIN = 1, function(x) minimax_design_old(p = c(x["p1_null"], x["p2_null"]), 
                                                                         n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  
  expect_is(trials, "data.frame")
})

#' Create a two_stage_curtail_trial
#' Case 4a:  User inputs p, n_total, prob_early
#' trials <- two_stage_curtail_trial(p1_null = 0.8, p2_null=0.2, 
#' p1_alt = 0.8, p2_alt = 0.4, n_total=36, prob_early=0.1, alpha=0.1)
test_that("two_stage_curtail_trial works for case 4a", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n_total <- 36
  prob_early=0.2

  trials <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                   p1_alt = p1_alt, p2_alt = p2_alt, n_total=n_total, 
                                   prob_early=prob_early)
  n1 <- seq_len(n_total-1)
  n2 <- n_total - n1
  n <- cbind(n1, n2)
  r <- matrix(apply(n, 1, function(x){
    two_stage_critical_values(x, c(p1_null, p2_null),
                              pearly = prob_early)
  }), nrow=length(n1), ncol=2, byrow=TRUE)
  
  r1 <- r[,1]
  r2 <- r[,2]
  df <- data.frame(p1_null=rep(p1_null, length(n1)),
                   p2_null=rep(p2_null, length(n1)),
                   p1_alt=rep(p1_alt, length(n1)),
                   p2_alt=rep(p2_alt, length(n1)),
                   n1=n1, n2=n2, r1=r[,1], r2=r[,2])
  if(min(df$r1)==0){
    df <- subset(df, r1>0)
    rownames(df) <- NULL
  }
  
  expect_equal(data.frame(trials[,1:8]), df)
  
  expect_equal(trials$power, apply(df, MARGIN = 1, function(x) two_stage_power(p = c(x["p1_alt"], x["p2_alt"]), 
                                                                               n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$significance, apply(df, MARGIN = 1, function(x) two_stage_significance(p = c(x["p1_null"], x["p2_null"]), 
                                                                                             n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$stage1_mean_ss, apply(df, MARGIN = 1, function(x) expected_stage1_sample_size(p = c(x["p1_null"], x["p2_null"]), 
                                                                                                    n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))$expectation))
  expect_equal(trials$mean_ss_null, apply(df, MARGIN = 1, function(x) expected_total_sample_size(p = c(x["p1_null"], x["p2_null"]), 
                                                                                                 n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$mean_ss_alt, apply(df, MARGIN = 1, function(x) expected_total_sample_size(p = c(x["p1_alt"], x["p2_alt"]), 
                                                                                                n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$PET, apply(df, MARGIN = 1, function(x) prob_early_stop(p = c(x["p1_null"], x["p2_null"]), 
                                                                             n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$minimax_prob, apply(df, MARGIN = 1, function(x) minimax_design_old(p = c(x["p1_null"], x["p2_null"]), 
                                                                                     n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  
  expect_is(trials, "data.frame")
})

#' Create a two_stage_curtail_trial
#' Case 4b:  User inputs p, n_total, alpha
#' trials <- two_stage_curtail_trial(p1_null = 0.8, p2_null=0.2, 
#' p1_alt = 0.8, p2_alt = 0.4, n_total=36, alpha=0.1)
test_that("two_stage_curtail_trial works for case 4b", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n_total <- 36
  alpha <- 0.05
  
  trials <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                    p1_alt = p1_alt, p2_alt = p2_alt, n_total=n_total, 
                                    alpha=alpha)
  n1 <- seq_len(n_total-1)
  n2 <- n_total - n1
  n <- cbind(n1, n2)
  r <- matrix(apply(n, 1, function(x){
    two_stage_critical_values(x, c(p1_null, p2_null),
                              alpha = alpha)
  }), nrow=length(n1), ncol=2, byrow=TRUE)
  
  r1 <- r[,1]
  r2 <- r[,2]
  df <- data.frame(p1_null=rep(p1_null, length(n1)),
                   p2_null=rep(p2_null, length(n1)),
                   p1_alt=rep(p1_alt, length(n1)),
                   p2_alt=rep(p2_alt, length(n1)),
                   n1=n1, n2=n2, r1=r[,1], r2=r[,2])
  if(min(df$r1)==0){
    df <- subset(df, r1>0)
    rownames(df) <- NULL
  }
  
  expect_equal(data.frame(trials[,1:8]), df)
  
  expect_equal(trials$power, apply(df, MARGIN = 1, function(x) two_stage_power(p = c(x["p1_alt"], x["p2_alt"]), 
                                                                               n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$significance, apply(df, MARGIN = 1, function(x) two_stage_significance(p = c(x["p1_null"], x["p2_null"]), 
                                                                                             n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$stage1_mean_ss, apply(df, MARGIN = 1, function(x) expected_stage1_sample_size(p = c(x["p1_null"], x["p2_null"]), 
                                                                                                    n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))$expectation))
  expect_equal(trials$mean_ss_null, apply(df, MARGIN = 1, function(x) expected_total_sample_size(p = c(x["p1_null"], x["p2_null"]), 
                                                                                                 n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$mean_ss_alt, apply(df, MARGIN = 1, function(x) expected_total_sample_size(p = c(x["p1_alt"], x["p2_alt"]), 
                                                                                                n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$PET, apply(df, MARGIN = 1, function(x) prob_early_stop(p = c(x["p1_null"], x["p2_null"]), 
                                                                             n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$minimax_prob, apply(df, MARGIN = 1, function(x) minimax_design_old(p = c(x["p1_null"], x["p2_null"]), 
                                                                                     n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  
  expect_is(trials, "data.frame")
})

#' Create a two_stage_curtail_trial
#' Case 5:  User inputs p, n_total
#' Using default values of prob_early and alpha, with n_total
#' trials <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n_total=36)
test_that("two_stage_curtail_trial works for case 5", {
  p1_null <- 0.8
  p2_null <- 0.2
  p1_alt <- 0.8
  p2_alt <- 0.4
  n_total <- 36

  trials <- two_stage_curtail_trial(p1_null = p1_null, p2_null = p2_null, 
                                    p1_alt = p1_alt, p2_alt = p2_alt, n_total=n_total)
                                  
  n1 <- seq_len(n_total-1)
  n2 <- n_total - n1
  n <- cbind(n1, n2)
  r <- matrix(apply(n, 1, function(x){
    two_stage_critical_values(x, c(p1_null, p2_null))
  }), nrow=length(n1), ncol=2, byrow=TRUE)
  
  r1 <- r[,1]
  r2 <- r[,2]
  df <- data.frame(p1_null=rep(p1_null, length(n1)),
                   p2_null=rep(p2_null, length(n1)),
                   p1_alt=rep(p1_alt, length(n1)),
                   p2_alt=rep(p2_alt, length(n1)),
                   n1=n1, n2=n2, r1=r[,1], r2=r[,2])
  if(min(df$r1)==0){
    df <- subset(df, r1>0)
    rownames(df) <- NULL
  }
  
  expect_equal(data.frame(trials[,1:8]), df)
  
  expect_equal(trials$power, apply(df, MARGIN = 1, function(x) two_stage_power(p = c(x["p1_alt"], x["p2_alt"]), 
                                                                               n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$significance, apply(df, MARGIN = 1, function(x) two_stage_significance(p = c(x["p1_null"], x["p2_null"]), 
                                                                                             n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$stage1_mean_ss, apply(df, MARGIN = 1, function(x) expected_stage1_sample_size(p = c(x["p1_null"], x["p2_null"]), 
                                                                                                    n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))$expectation))
  expect_equal(trials$mean_ss_null, apply(df, MARGIN = 1, function(x) expected_total_sample_size(p = c(x["p1_null"], x["p2_null"]), 
                                                                                                 n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$mean_ss_alt, apply(df, MARGIN = 1, function(x) expected_total_sample_size(p = c(x["p1_alt"], x["p2_alt"]), 
                                                                                                n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$PET, apply(df, MARGIN = 1, function(x) prob_early_stop(p = c(x["p1_null"], x["p2_null"]), 
                                                                             n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  expect_equal(trials$minimax_prob, apply(df, MARGIN = 1, function(x) minimax_design_old(p = c(x["p1_null"], x["p2_null"]), 
                                                                                     n = c(x["n1"], x["n2"]), r = c(x["r1"], x["r2"]))))
  
  expect_is(trials, "data.frame")
})



##private functions?
# 
# test_that("single_stage_significance function works", {
#   resp <- single_stage_significance(.2, 5, 7) 
#   expect_equal(round(resp, 8), 0.05040957)
#   expect_is(resp, "numeric")
# })
# 
# test_that("single_stage_power function works", {
#   resp <- single_stage_power(.6, 5, 7) 
#   expect_equal(round(resp, 7), 0.9006474)
#   expect_is(resp, "numeric")
# })
# 
# test_that("single_stage_expected_sample_size function works", {
#   resp <- single_stage_expected_sample_size(.2, 7, 11)
#   expect_equal(round(resp$expectation, 5), 13.61483)
#   expect_equal(round(resp$standardDeviation, 6), 1.627825)
#   expect_is(resp, 'list')
# })
# 
# test_that("dsnb_stacked function works", {
#   s <- 5
#   t <- 7
#   resp <- dsnb_stacked(x = min(s,t):(s+t-1), p = .2, s = s, t = t)
#   expect_equal(resp[, "x"], 5:11)
#   expect_equal(round(resp[, "s"], 6), c(0.00032, 0.00128, 0.003072, 0.005734, 0.009175, 0.013212, 0.017616))
#   expect_equal(round(resp[, "t"], 6), c(0, 0, 0.209715, 0.293601, 0.234881, 0.140929, 0.070464))
#   expect_is(resp, "matrix")
# })
# test_that("cdsnb_stacked function works", {
#   s <- 5
#   t <- 7
#   resp <- cdsnb_stacked(x=min(s,t):(s+t-1), shape=c(1, 1), s=s, t=t)
#   expect_equal(resp[, "k"], 5:11)
#   expect_equal(round(resp[, "s"], 6), c(0.166667, 0.119048, 0.089286, 0.069444, 0.055556, 0.045455, 0.037879))
#   expect_equal(round(resp[, "t"], 6), c(0, 0, 0.125, 0.097222, 0.077778, 0.063636, 0.053030))
#   expect_is(resp, "data.frame")
# })
# 
# test_that("esnb function works", {
#   resp <- esnb(p = .4, s = 7, t=11)
#   expect_equal(round(resp, 5), 14.50153)
#   expect_is(resp, "numeric")
# })
# 
# test_that("vsnb function works", {
#   resp <- vsnb(p = .4, s = 7, t=11)
#   expect_equal(round(resp, 6), 4.544609)
#   expect_is(resp, "numeric")
# })
# 
# test_that("ecsnb function works", {
#   resp <- ecsnb(shape=c(1, 1), s=7, t=11)
#   expect_equal(round(resp, 5), 11.54329)
#   expect_is(resp, "numeric")
# })
# 
# test_that("vcsnb function works", {
#   resp <- vcsnb(shape=c(1, 1), s=7, t=11)
#   expect_equal(round(resp, 6), 9.209055)
#   expect_is(resp, "numeric")
# })
# 
# test_that("critical_values function works for two-stage design", {
#   resp <- critical_values(n=c(6, 30), p=c(.8, .2), pearly = .1, alpha = .1)
#   expect_equal(resp[1], 4)
#   expect_equal(resp[2], 11)
#   expect_is(resp, 'numeric')
# })
# 
# test_that("critical_values function works for one-stage design", {
#   resp <- critical_values(n=36, p=.2, alpha = .1)
#   expect_equal(resp,11)
#   expect_is(resp, 'numeric')
# })
# 
# test_that("prob_early_stop function works", {
#   resp <- prob_early_stop(p=0.8, n = 6, r = 4)
#   expect_equal(round(resp, 5), 0.09888)
#   expect_is(resp, 'numeric')
# })
# 
# test_that("expected_stage1_sample_size function works", {
#   resp <- expected_stage1_sample_size(p = c(0.8, .2), n = c(6, 30), r = c(4, 11))
#   expect_equal(resp$expectation, 4.76)
#   expect_equal(round(resp$standardDeviation, 7), 0.7797435)
#  expect_is(resp, 'list')
# })
# 
# test_that("expected_total_sample_size function works", {
#   resp <- expected_total_sample_size(p = c(0.8, 0.2), n= c(6, 30), r = c(4, 11))
#   expect_equal(round(resp, 5), 29.32231)
#   expect_is(resp, 'numeric')
# })
# 
# ### optimal and minimax functions
# 
# test_that("best_designs function works", {
#   resp <- best_designs(p = c(0.8, 0.2), ntot = 36, pearly = .1, alpha = .1)$designs
#   expect_equal(resp["Optimal", "p1"], 0.8)
#   expect_equal(resp["Optimal", "n1"], 6)
#   expect_equal(resp["Optimal", "r1"], 4)
#   expect_equal(resp["Optimal", "p2"], 0.2)
#   expect_equal(resp["Optimal", "n2"], 30)
#   expect_equal(resp["Optimal", "r2"], 11)
#   expect_equal(resp["Optimal", "Alpha"], 0.1)
#   expect_equal(round(resp["Optimal", "PET"], 5), 0.09888)
#   expect_equal(round(resp["Optimal", "ECSS"], 5), 29.32231)
#   expect_equal(round(resp["Optimal", "P(MaxSS)"], 8), 0.06584612)
# 
#   expect_equal(resp["Minimax", "p1"], 0.8)
#   expect_equal(resp["Minimax", "n1"], 6)
#   expect_equal(resp["Minimax", "r1"], 4)
#   expect_equal(resp["Minimax", "p2"], 0.2)
#   expect_equal(resp["Minimax", "n2"], 30)
#   expect_equal(resp["Minimax", "r2"], 11)
#   expect_equal(resp["Minimax", "Alpha"], 0.1)
#   expect_equal(round(resp["Minimax", "PET"], 5), 0.09888)
#   expect_equal(round(resp["Minimax", "ECSS"], 5), 29.32231)
#   expect_equal(round(resp["Minimax", "P(MaxSS)"], 8), 0.06584612)
#   
#   expect_is(resp, 'data.frame')
# })
# 
# test_that("minimax_design_old function works", {
#   resp <- minimax_design_old(p = c(.8, .2), n = c(3,33), r = critical_values(n=c(3,33), p=c(.8, .2), pearly = .1, alpha =.1))
#   expect_equal(round(resp, 5), 0.07063)
#   expect_is(resp, 'numeric')
# })
# 
# 
# test_that("two_stage_significance function works", {
#   resp <- two_stage_significance(p = c( .8, .2), n = c(12, 24), r = c(8, 11))
#   expect_equal(round(resp, 5), 0.08578)
#   expect_is(resp, 'numeric')
# })
# 
# test_that("two_stage_power function works", {
#   resp <- two_stage_power(p = c( .8, .4), n = c(12, 24), r = c(8, 11))
#   expect_equal(round(resp, 6), 0.850723)
#   expect_is(resp, 'numeric')
# })
# 
# test_that("all_optimal_designs function works", {
#   resp <- all_optimal_designs(p = c( .8, .2), n = 10)
#   expect_equal(resp$n1, c(2, 6, 5, 4, 3, 7, 8, 9))
#   expect_equal(resp$n2, c(8, 4, 5, 6, 7, 3, 2, 1))
#   expect_equal(round(resp$`Expected Sample Size`, 6), c(7.214133, 7.234434, 7.272286, 7.329535, 7.383349, 7.426047, 7.587027, 7.986812))
#   expect_is(resp, 'data.frame')
# })
# 
# test_that("all_minimax_designs function works", {
#   resp <- all_minimax_designs(p = c( .8, .2), n = 10)
#   expect_equal(resp$n1, c(6, 2, 5, 9, 8, 4, 7, 3))
#   expect_equal(resp$n2, c(4, 8, 5, 1, 2, 6, 3, 7))
#   expect_equal(round(resp$`Probability of Maximum Sample Size`, 8), c(0.06404506, 0.06491341, 0.06498202, 0.06502810, 0.06545818, 0.06559949, 0.06574490, 0.06593741))
#   expect_is(resp, 'data.frame')
# })
# 
# test_that("ph2snb function works", {
#   resp <- ph2snb(p = .8, s = 3, t = 5)
#   expect_equal(round(resp[1], 5), 0.00000)
#   expect_equal(round(resp[2], 5), 0.00000)
#   expect_equal(round(resp[3], 5), 0.51200)
#   expect_equal(round(resp[4], 5), 0.30720)
#   expect_equal(round(resp[5], 5), 0.12320)
#   expect_equal(round(resp[6], 5), 0.04224)
#   expect_equal(round(resp[7], 5), 0.01536)
#   expect_is(resp, 'numeric')
# })
# 
# test_that("ph2valid function works", {
#   resp <- ph2valid(p = c( .2, .3), n = c(10, 10), r = c(5, 5))
#   expect_equal(resp, FALSE)
#   expect_is(resp, 'logical')
# })
# 
# test_that("ph2tnb function works", {
#   resp <- ph2tnb(p = c( .8, .2), n = c(6, 30), r = c(4, 11))
#   expect_equal(resp[1], 0)
#   expect_equal(resp[2], 0)
#   expect_equal(resp[3], 0)
#   expect_equal(round(resp[4], 7), 0.4545455)
#   expect_equal(round(resp[5], 7), 0.3636364)
#   expect_equal(round(resp[6], 7), 0.1818182)
#   expect_is(resp, 'numeric')
# })
# 
# test_that("is.wholenumber function works", {
#   resp <- is.wholenumber(2.365)
#   expect_equal(resp, FALSE)
#   expect_is(resp, 'logical')
# })
# 
# test_that("prob_reject_traditional function works", {
#   resp <- prob_reject_traditional(p = c( .8, .2), n = c(12, 24), r = c(8, 11))
#   expect_equal(round(resp, 5), 0.08578)
#   expect_is(resp, 'numeric')
# })
# 
# test_that("prob_reject function works", {
#   resp <- prob_reject(p = c( .8, .2), n = c(12, 24), r = c(8, 11))
#   expect_equal(round(resp, 5), 0.08578)
#   expect_is(resp, 'numeric')
# })
