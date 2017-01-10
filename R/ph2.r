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
ph2crit = function(n, p, pearly = .1, alpha = .1) {
  
  if (!ph2valid(p, n, r = c(0, 0)) ||
      pearly < 0 || pearly > 1 || alpha < 0 || alpha > 1)
  {
    warning("Invalid parameter values (PhII)")
    return (NaN)
  }
  
  # Binomial critical value, r2
  r2 <- 1 + qbinom(1 - alpha, sum(n), p[2])     
  
  for (j in 0 : n[1]){
    
    # r1 based on pearly
    if (ph2early(p, n, r = c(j, 0)) > pearly)
      
      return(c(max(j - 1, 0), r2))
  }
  
  return (c(n[1], r2))         # Hmmm. Didn't add up to 1
}

#' Probability of stopping early 
#'
#' Calculates the probability of stopping the trial after Stage 1
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
ph2early = function( p, n, r) {
  
  if(length(p) == 2) p1 <- p[1] # p can be vector or scalar
  
  else p1 <- p
  
  if(length(n) == 2) n1 <- n[1] # n can be vector or scalar
  
  else n1 <- n
  
  if(length(r) == 2) r1 <- r[1] # r can be vector or scalar
  
  else r1 <- r
  
  
  if(n1 <= 0) return(0)    # no patients -> no event
  
  if(r1 <= 0) return(0)    # avoids some rounding errors
  
  if(r1 > n1) return(1)    # always stop early
  
  j <- min(r1, n1) : n1    # range of X1
  
  ph2early <- 1 - sum(dbinom(j, n1, p1))
  
  return(ph2early)
}

#' Expected sample size and SD for decision to continue to Stage 2
#'
#' Mean and standard deviation of the minimum number of Stage 1 patients necessary to be able
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
ph2Eearly = function(p, n1, r1) {
  if(length(p) > 1){ p1 <- p[1] # p can be vector or scalar
  
  }else p1 <- p
  
  if(length(n) > 1){ n1 <- n[1] # n can be vector or scalar
  
  }else n1 <- n
  
  if(length(r) > 1){ r1 <- r[1] # r can be vector or scalar
  
  }else r1 <- r
  
  
  
  tt <- n1 - r1 + 1  # Number of failures to stop early
  
  e <- sum((1 : n1) * ph2dspb(p1, r1, tt))  # Expected value
  
  v <- sum((1 : n1) ^ 2 * ph2dspb(p1, r1, tt)) - e ^ 2   # variance
  
  return(c(e, sqrt(v)))  # Expected value and SD
  
}

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
ph2Ess = function(p, n, r) {
  
  # check parameter values
  if( ! ph2valid(p, n, r)){
    warning("Invalid parameter values (PhII)")
    return(NaN)
  }
  
  pearly <- ph2early( p, n, r)  # Probability of stopping early
  
  #Pr(Y1=y1|continue to stage 2)
  # truncated negative binomial distribution
  # number of trials until we reach r1 successes 
  # r1 <= y1 <= n1
  prY1Continue <- function(y1, r1, n1, p1){
    if(y1<r1 || y1>n1){    
      warning("Invalid parameter values (y1)")
      return(NaN)
    }
    num <- factorial(y1-1)/factorial(y1-r1)
    denom <- 0 
    for(j in r1:n1){
      a <- factorial(j-1)/factorial(j-r1) *(1-p1)^(j-y1)
      denom <- denom + a
    }
    return(num/denom)
  }
  # Pr(Y1=y1|stop early)
  # truncated negative binomial distribution
  # number of trials until we reach n1-r1+1 failures
  # n1-r1+1 <= y1 <= n1
  prY1Stop <- function(y1, r1, n1, p1){
    if(y1 < n1-r1+1 || y1 > n1 ){ 
      warning("Invalid parameter values (y1)")
      return(NaN)
    }
    num <- factorial(y1-1)/factorial(y1-1-n1+r1)
    denom <- 0
    for(j in (n1-r1+1):n1){
      a <- factorial(j-1)/factorial(j-1-n1+r1) * p1^(j-y1)
      denom <- denom+a
    }
    if(is.na(num/denom)){ return(0)}
    return(num/denom)
  }
  
  
  ### expected value of Y1 | continue to Stage 2
  EY1Continue <- 0
  for(i in r[1]:n[1]){
    EY1Continue <- EY1Continue + i*prY1Continue(i, r[1], n[1], p[1])
  }
  
  
  ### expected value of Y1|stop early
  EY1Stop <- 0
  for(i in (n[1]-r[1]+1):n[1]){
    EY1Stop <- EY1Stop + i*prY1Stop(i, r[1], n[1], p[1])
  }
  ### expected value of y2|continue to Stage 2
  a <- 0
  b<- NULL
  
  for(i in r[1]:n[1]){
    py1 <- prY1Continue(i, r[1], n[1], p[1])
    a <- 0
    for(j in 0:r[1]){
      px12 <- dbinom(j, r[1], p[2]/p[1])
      s <- r[2]-j
      t <- n[2]+n[1]-i-r[2]+j+1
      if(s < 0) s <- 0
      if(t < 0) t <- 0
      ey2 <- esnb(p[2], s, t)
      
      a <- a+ px12*ey2
    }
    b <- c(b, (py1*a))
  }
  EY2Continue <- sum(b)
  
  return((EY1Stop*pearly) + ((EY1Continue+EY2Continue)*(1-pearly)))

}

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
ph2mmax = function(p, n, r) {
  if (!ph2valid(p, n, r)) return(NaN) # Check validity of parameters
  if (r[1]==0){
    return(choose(n[1]+n[2]-1, r[2]-1)*p[2]^(r[2]-1)*(1-p[2])^(n[1]+n[2]-r[2]))
  }
  else{
    pearly <- ph2early(p, n, r) # Probability of stopping early
    ph2mmax <- 0
    
    ck <- rep(0, n[1])     # Distribution of Y1  |  Don't stop early
    k <- 1 : n[1]
    
    ck <- p[1] ^ r[1] * (1 - p[1]) ^ (k - r[1]) * choose(k - 1, r[1] - 1)
    ck <- ck / sum(ck[r[1]:n[1]]) 
    
    for (j in r[1] : n[1])
    {
      y1j <- ck[j]           #  Pr[ Y1 = j | don't stop early]
      conv <- 0
      for (i in 0 : r[1])     # Y2 is convoluion sum of two binomials
      {
        x12 <- dbinom(i, r[1], p[2] / p[1])
        xp2 <- dbinom(r[2] - i - 1, sum(n) - j - 1, p[2])
        conv <- conv + x12 * xp2
      }
      ph2mmax <- ph2mmax + conv * y1j
    }
    ph2mmax * (1 - pearly)
  }
}

#' Probability of rejecting the null hypothesis under traditional sampling
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
ph2reject = function(p, n, r) {
  # check validity of parameter values
  
  if ( ! ph2valid(p, n, r)){
  
    warning("Invalid parameter values (PhII)")
    
    return (NaN)
    
  }
  
  
  
  reject <- 0
  
  
  
  # Loop on X1 = number of Stage 1 successes
  
  lowx1 <- max(r[1], r[2] - n[2])    # lower limit for X1
  
  for (x1 in lowx1 : n[1])
    
  {
    
    t1 <- dbinom(x1, n[1], p[1]) # binomial probability of X1
    
    
    
    # Loop on X12 = Stage 1 successes who become Stage 2 successes
    
    lowx12 <- max(0, r[2] - n[2])  # lower limit for X12           
    
    for (x12 in lowx12 : x1)       # X12 conditional on X1
      
    {
      
      t2 <- dbinom(x12, x1, p[2] / p[1])  #  Prob of X12 given X1
      
      
      
      # Loop on X2 = Stage 2 successes
      
      lowx2 <- max(0, r[2] - x12) # lower limit for X2
      
      for (x2 in lowx2 : n[2])
        
      {
        
        reject <- reject + t1 * t2 * dbinom(x2, n[2], p[2])   
        
      }
      
    }
    
  }
  
  return(reject)
}

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
ph2rejcs = function(p, n, r) {
  # check validity of parameter values
  if ( ! ph2valid(p, n, r)){
  
    warning("Invalid parameter values (PhII)")
    return (NaN)
  }
  
  reject <- 0 
  
  # Loop on Y1 = curtailed sample size in Stage 1
  
  
  for (y1 in r[1] : n[1])     #  First summation 
  {
    t1 <- ph2tnb(p[1], n[1], r[1])[y1]
    
    # Loop on X12 = Stage 1 successes who become Stage 2 successes
    lowx12 <- max(0, r[2] - n[2])  # lower limit for X12            
    for (x12 in lowx12 : r[1])       # X12 conditional on X1
    {
      t2 <- dbinom(x12, y1, p[2] / p[1])  #  Prob of X12 given Y1
      
      # Loop on X2 = Stage 2 successes - Third summation
      lowx2 <- max(0, r[2] - x12) # lower limit for X2
      for (x2 in lowx2 : (sum(n) - y1))
      {
        reject <- reject + t1 * t2 * dbinom(x2, n[2], p[2])    
      }
    }
  }
  return(reject)
}


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
ph2dspb = function(p, s, t) {
  
  if( length(p) != 1 || p < 0 || p > 1 || s < 1 || t < 1 ||
      !is.wholenumber(s) || !is.wholenumber(t)){
    
    warning("Invalid parameter values (PhII)")
    
    return(NaN)
    
  }
  
  tmnb <- rep(0, s + t - 1)    # zero out, over the range
  
  for(j in 0 : (t - 1)){      # Last event is success:  eq # (6)
    
    tmnb[j + s] <- choose(s + j - 1, s - 1) * p ^ s * (1 - p) ^ j
    
  }
  
  
  for(j in 0 : (s - 1)){       # Last event is failure:  eq # (7)
    
    tmnb[j + t] <- tmnb[j + t] +
      
      choose(t + j - 1, t - 1) * p ^ j * (1 - p) ^ t
    
  }
  
  return(tmnb / sum(tmnb))    # Normalize
}


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
ph2tnb = function(p, n, r) {
  if(length(p) > 1){ p1 <- p[1] # p can be vector or scalar
  
  }else p1 <- p
  
  if(length(n) > 1){ n1 <- n[1] # n can be vector or scalar
  
  }else n1 <- n
  
  if(length(r) > 1){ r1 <- r[1] # r can be vector or scalar
  
  }else r1 <- r
  
  if(p1 < 0 || p1 > 1 ||
     
     r1 < 1 || n1 < 1 ||
     
     !is.wholenumber(r1))
    
  {
    
    warning("Invalid parameter values (PhII)")
    
    return(NaN)
    
  }
  tnb <- NULL
  for(y1 in 1:n1){
    if(y1-r1 >= 0){
      num <- factorial(y1-1)/factorial(y1-r1)
    } else num <- 0
    denom <- 0 
    for(j in r1:n1){
      a <- factorial(j-1)/factorial(j-r1) *(1-p1)^(j-y1)
      denom <- denom + a
    }
    tnb <- c(tnb, (num/denom))
    
  }
  return(tnb)
}

#' Check validity of parameter values p, n, and r.
#' 
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2). Must have 0 <= p2 <= p1 <= 1.
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2).  Both n1 and n2 should be non-negative integers.
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2). Must have 0 <= r1 <= r2,
#' 0 <= r1 <= n1, and 0 <= r2 <= n1 + n2.
#' @examples
#' ph2valid(p = c(.2, .1), n = c(10, 10), r = c(5, 5))
#' ph2valid(p = c( .2, .3), n = c(10, 10), r = c(5, 5))
#' @export
ph2valid = function(p,n,r) {
  #  Must have:  0 <= p2 <= p1 <= 1
  if( length(p) != 2 ) return(FALSE)
  
  if( p[1] > 1 | p[2] < 0 | p[1] < p[2] )return(FALSE)
  
  
  
  #  The two n's must be non-negative integers
  if( length(n) != 2) return(FALSE)
  
  if(min(n) < 0) return(FALSE)
  
  if(!is.wholenumber(n[1])
     
     || !is.wholenumber(n[2])) return(FALSE)
  
  
  
  # Must have: 0 <= r1 <= n1
  if(length(r) != 2) return(FALSE)
  
  if( r[1] < 0) return(FALSE)
  
  if( r[1] > n[1]) return(FALSE)
  
  
  
  # Must have: 0 <= r2 <= n1+n2
  if(r[2] < 0) return(FALSE)
  
  if(r[2] > sum(n)) return(FALSE)
  
  if(!is.wholenumber(r[1]) || !is.wholenumber(r[2])) return(FALSE) 
  
  return(TRUE)               # Valid parameter values
  
}

#' Test for integer
#' 
#' @param x 
#' @param tot 
#' @examples
#' is.wholenumber(4)
#' is.wholenumber(3.5)
is.wholenumber = function(x, tol = .Machine$double.eps ^ 0.5){
    abs(x - round(x)) < tol
}
