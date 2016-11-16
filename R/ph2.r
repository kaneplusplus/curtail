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

#' Expected sample size for decision to continue to Stage 2
#'
#' DESCRIPTION HERE
#'
#' @param p
#' @param n1
#' @param r1
#' @examples
#' # Put exmple code here.
#' @export
ph2Eearly = function(p, n1, r1) {}

#' Expected curtailed sample size
#'
#' DESCRIPTION HERE
#'
#' @param p p
#' @param n n
#' @param r r
#' @examples
#' # Put exmple code here.
#' @export
ph2Ess = function( p, n, r) {}

#' Evaluate the probability that the maximum sample size is needed
#'
#' DESCRIPTION HERE
#'
#' @param p p
#' @param n n
#' @param r r
#' @examples
#' # Put exmple code here.
#' @export
ph2mmax = function(p, n, r) {}

#' Probability of rejecting the null hypothesis
#'
#' DESCRIPTION HERE
#'
#' @param p p
#' @param n n
#' @param r r
#' @examples
#' # Put exmple code here.
#' @export
ph2reject = function(p, n, r) {}

#' Probability of rejecting the null hypothesis under curtail ed sampling
#'
#' DESCRIPTION HERE
#'
#' @param p p
#' @param n n
#' @param r r
#' @examples
#' # Put exmple code here.
#' @export
ph2rejcs = function(p, n, r) {}

#' Truncated negative binomial distribution mass function
#'
#' DESCRIPTION HERE
#'
#' @param p p
#' @param n n
#' @param r r
#' @examples
#' # Put exmple code here.
#' @export
ph2tnb = function(p, n, r) {}

#' Check validity of parameter values
#'
#' DESCRIPTION HERE
#'
#' @param p p
#' @param n n
#' @param r r
#' @examples
#' # Put exmple code here.
#' @export
ph2valid = function(p,n,r) {}

