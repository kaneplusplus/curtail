#' Find critical values (r1, r2) 
#' 
#' Finds r1, the minimum number of Stage 1 successes to continue to 
#' Stage 2, and r2, the minimum number of Stage 2 successes to 
#' reject the null hypothesis.
#' 
#' Find r1 so that the probability of early stopping
#' is less than or equal to pearly. 
#' Find r2 from Binomial model with no early stopping
#' and significance level alpha.
#'
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'          and Stage 2 (n2).
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2).
#' @param pearly desired probability of early stopping (default = .1).
#' @param alpha desired significance level (default = .1).
#' @examples
#' ph2crit( n = c( 5, 31), p = c(.8,.2), pearly = .1, alpha = .1)
#' @export
ph2crit = function(n, p, pearly = .1, alpha = .1) {}

#' Stopped negative binomial distribution mass function
#'
#' Returns the probability density function for minimum number of 
#' events for either s successes or t failures.
#'
#' @param p probability of success in each trial. 
#' @param s number of successes.
#' @param t number of failures.
#' @examples
#' ph2dspb(p = .8, s = 3, t = 5)
#' @export
ph2dspb = function(p, s, t) {}

#' Probability of stopping early 
#'
#' Probability of stopping the trial after Stage 1
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) or a scalar containing p1
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'          and Stage 2 (n2) or a scalar containing n1
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2) or a scalar containing r2
#' @examples
#' ph2early(p = .8, n = 25, r = 18)
#' ph2early(p = c(.8, .2), n = c(25, 9), r = 18)
#' @export
ph2early = function( p, n, r) {}

#' Expected sample size and SD for decision to continue to Stage 2
#'
#' Mean and SD of the minimum number of Stage 1 patients necessary to be able
#' to decide whether to either continue to Stage 2 or else terminate early
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) or a scalar containing p1
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) or a scalar containing n1
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2) or a scalar containing r1
#' @examples
#' ph2Eearly( p = .8, n = 5, r = 3)
#' ph2early(p = c(.8, .2), n = c(5, 31), r = 3)
#' @export
ph2Eearly = function(p, n1, r1) {}

#' Expected curtailed sample size
#'
#' The expected number of patients who are enrolled and followed to their
#' endpoint before critical endpoints in Stage 1 and Stage 2 are achieved
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
#' @examples
#' ph2Ess(p = c(.8, .2), n=c(6, 30), r=c(4, 11))
#' ph2Ess(p = c(.8, .2), n = c(18, 18), r = c(12, 11))
#' @export
ph2Ess = function(p, n, r) {}

#' Evaluate the probability that the maximum sample size is needed
#'
#' Evaluate the probability that the maximum sample size n1+n2 is
#' required to complete the trial
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
#' @examples
#' ph2mmax( p = c(.8, .2), n = c(3,33), r = ph2crit(n=c(3,33), p=c(.8, .2), pearly = .1, alpha =.1))
#' @export
ph2mmax = function(p, n, r) {}

#' Probability of rejecting the null hypothesis under traditional sampling
#'
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
#' @examples
#' ph2reject(p = c( .8, .2), n = c(12, 24), r = c(8, 11))
#' ph2reject(p = c( .8, .2), n = c(12, 24), r = c(0, 11))
ph2reject = function(p, n, r) {}

#' Probability of rejecting the null hypothesis under curtailed sampling
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
#' @examples
#' # Put exmple code here.
#' @export
ph2rejcs = function(p, n, r) {}

#' Truncated negative binomial distribution mass function
#'
#' Returns the probability density function for the number of events
#' to reach r= r1 successes in Stage 1 with p = p1 probability of success.
#' Can also be the number of events to reach r failures with p = (1-p1) probability
#' of failure.
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) or a scalar containing p1
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) or a scalar containing n1
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2) or a scalar containing r1
#' @examples
#' ph2tnb(p = .8, n = 6, r = 4) 
#` ph2tnb(p = .2, n = 6, r = 3) 
#' @export
ph2tnb = function(p, n, r) {}

#' Check validity of parameter values
#'
#' DESCRIPTION HERE
#'
#' @param p p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2)
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
#' @examples
#' ph2valid(p = c(.2, .1), n = c(10, 10), r = c(5, 5))
#' ph2valid(p = c( .2, .3), n = c(10, 10), r = c(5, 5))
#' @export
ph2valid = function(p,n,r) {}

