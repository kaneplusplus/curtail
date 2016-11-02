############################################################

############################################################

#                                                          #

#   Two-stage, Phase II clinical trial design with         #

#   nested criteria for early stopping and efficacy        #

#                                                          #

############################################################

############################################################





#  List of functions defined here:



# ph2crit( n, p, pearly = .1, alpha = .1)  # Find critical values (r1, r2)

# ph2dspb(p, s, t)     # Stopped negative binomial distribution mass function

# ph2early( p, n, r)   # Probability of stopping early 

# ph2Eearly(p, n1, r1) # Expected sample size for decision to continue to Stage 2

# ph2Ess( p, n, r)     # Expected curtailed sample size

# ph2mmax(p, n, r)     # Evaluate the probability that the maximum sample size is needed

# ph2reject(p, n, r)   #  Probability of rejecting the null hypothesis

# ph2rejcs(p, n, r) # Probability of rejecting the null hypothesis under curtailed sampling

# ph2tnb(p, n, r)      # Truncated negative binomial distribution mass function

# ph2valid(p,n,r)      # Check validity of parameter values


##################################################################


# Dependencies:



library(clinfun)       # For Simon, Two-stage Phase II designs

library(snb)           # For stopped negative binomial distribution



#############################################################

ph2tnb <- function(p, n, r){
  
  # Truncated negative binomial distribution mass function  
  
  # Returns the probability density function for the number of events to
  
  # reach r = r1 successes (p = probability of success = p1) in Stage 1.
  
  # Can also be the number of events to reach r failures (p = probability of failure = 1-p1).
  
  
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





# Test for integer: From the R Help file



is.wholenumber <-
  
  function(x, tol = .Machine$double.eps ^ 0.5)
    
    abs(x - round(x)) < tol



# # # # # # # # # # # # # # # # # # # # # # # # # # #



is.wholenumber(5.0000)          # Not a test of integer

is.wholenumber(5 + 1.0e-7)      # Small but not whole number

is.wholenumber(5 + 1.0e-8)      # Small enough to be whole number



#############################################################





ph2valid <- function(p,n,r)  # Check validity of parameter values
  
{
  
  
  
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
  
  if(!is.wholenumber(r[1])
     
     || !is.wholenumber(r[2])) return(FALSE) 
  
  
  
  TRUE               # Valid parameter values
  
}



# # # # # # # # Test Validity of parameters



# Test ph2valid in the middle of the ranges:

ph2valid(p = c(.2, .1), n = c(10, 10), r = c(5, 5))  # TRUE



# Test validity at the boundaries:

ph2valid(p = c(.2, .2), n = c(10, 10), r = c(5, 5))  # TRUE: p1=p1

ph2valid(p = c(.0, .0), n = c(10, 10), r = c(5, 5))  # TRUE: p1=p2=0

ph2valid(p = c( 1,  1), n = c(10, 10), r = c(5, 5))  # TRUE: p1=p2=1

ph2valid(p = c(.2, .1), n = c(0,  10), r = c(0, 5))  # TRUE: r1=n1=0

ph2valid(p = c(.2, .1), n = c(10, 0),  r = c(5, 6))  # TRUE: n2=0

ph2valid(p = c(.2, .1), n = c(10, 10), r = c(5, 20)) # TRUE: r2=n1+n2

ph2valid(p = c(.2, .1), n = c(10, 10), r = c(0, 5))  # TRUE: r1=0

ph2valid(p = c(.2, .1), n = c(10, 10), r = c(5, 0))  # TRUE: r2=0



# Test invalid cases:

ph2valid(p =    .2,      n = c(10, 10), r = c(5, 5))  # FALSE: wrong dimensions

ph2valid(p = c( .2, .3), n =   10,      r = c(5, 5))  # FALSE: wrong dimensions

ph2valid(p = c( .2, .3), n = c(10, 10), r =   5    )  # FALSE: wrong dimensions

ph2valid(p = c( .2, .3), n = c(10, 10), r = c(5, 5))  # FALSE: p2 > p1

ph2valid(p = c(-.1, .3), n = c(10, 10), r = c(5, 5))  # FALSE: p1 < 0

ph2valid(p = c(1.2, .3), n = c(10, 10), r = c(5, 5))  # FALSE: p1 > 1

ph2valid(p = c(.2, -.1), n = c(10, 10), r = c(5, 5))  # FALSE: p2 < 0

ph2valid(p = c(.2, 1.3), n = c(10, 10), r = c(5, 5))  # FALSE: p2 > 1

ph2valid(p = c(.2,  .1), n = c(-1, 10), r = c(5, 5))  # FALSE: n1 < 0

ph2valid(p = c(.2,  .1), n = c(10, -1), r = c(5, 5))  # FALSE: n1 < 0

ph2valid(p = c(.2,  .1), n = c(1.2, 1), r = c(5, 5))  # FALSE: n1 not integer

ph2valid(p = c(.2,  .1), n = c(10,1.2), r = c(5, 5))  # FALSE: n2 not integer

ph2valid(p = c(.2,  .1), n = c(10, 10), r = c(-1, 5)) # FALSE: r1 < 0

ph2valid(p = c(.2,  .1), n = c(10, 10), r = c(15, 5)) # FALSE: r1 > n1

ph2valid(p = c(.2,  .1), n = c(10, 10), r = c(5, 35)) # FALSE: r2 > n1+n2

ph2valid(p = c(.2,  .1), n = c(10, 10), r = c(5.1, 5))# FALSE: r1 not integer

ph2valid(p = c(.2,  .1), n = c(10, 1.2),r = c(5, 5.1))# FALSE: r2 not integer





######################################################





ph2reject <- function(p, n, r)
  
  # Probability of rejecting the null hypothesis
  
  # Expression (3) in the paper
  
  # Input: (p1,p2), (n1,n2), (r1,r2)

{
  
  # check validity of parameter values
  
  if ( ! ph2valid(p, n, r))
    
  {
    
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
  
  reject
  
}



# # # # Test ph2reject



ph2reject(p = c( .8, .2), n = c(8, 28), r = c(5, 11)) # Test: 0.08623555

ph2reject(p = c( .8, .2), n = c(12, 24), r = c(8, 11)) # Test: 0.08577967

ph2reject(p = c( .8, .2), n = c(12, 24), r = c(0, 11)) # Test: binomial design

1 - pbinom(10, 36, .2)             #  Confirm binomial: 0.08891278



################################################################

ph2rejcs <- function(p, n, r)
  
  # Probability of rejecting the null hypothesis under curtailed sampling 
  # Input: (p1,p2), (n1,n2), (r1,r2)

{ 
  # check validity of parameter values
  if ( ! ph2valid(p, n, r))
  {
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
  reject
}

##########
ph2rejcs(c(.8, .2), c(6, 30), c(4, 11))
################################################################


ph2early <- function( p, n, r)
  
  
  
  #  Probability of stopping early
  
  #  Returns Pr[ X_1 < r_1 ]
  
  
  
{
  
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



# # # #



ph2early(p = .8, n = 25, r = 18) # .1091228, Test: one p, one n

ph2early(p = c(.8, .2), n = c(25, 9), r = 18)  # .1091228, test 2 n's, 2 p's

ph2early(p = .8, n = 5, r = 3)  # Design A: .05792

ph2early(p = .8, n = 8, r = 5)  # Design B: .0562816

ph2early(p = .8, n = 11, r = 7)  # Design C: .05040957

ph2early(p = .8, n = 12, r = 8)  # Design D: .0725555





########################################################


ph2crit <- function( n, p, pearly = .1, alpha = .1)
  
  
  
  # Find critical values (r1, r2)
  
  # Find r1 so that the probability of early stopping
  
  # is less than or equal to pearly. 
  
  # Find r2 from Binomial model with no early stopping,
  
  # and significance level alpha.



{                    # check validity of parameter values
  
  if (!ph2valid(p, n, r = c(0, 0)) ||
        
        pearly < 0 || pearly > 1 || alpha < 0 || alpha > 1)
    
  {
    
    warning("Invalid parameter values (PhII)")
    
    return (NaN)
    
  }
  
  
  
  r2 <- 1 + qbinom(1 - alpha,
                   
                   sum(n), p[2])          # binomial critical value
  
  
  
  for (j in 0 : n[1])
    
  {
    
    if (ph2early(p, n, r = c(j, 0)) > pearly)
      
      return(c(max(j - 1, 0), r2))
    
  }
  
  return (c(n[1], r2))         # Hmmm. Didn't add up to 1
  
}



# # # # # # Test validity of input parameters



ph2crit( n = 5, p = c(.8, .2), pearly = .1, alpha = .1)  # Fail: one n

ph2crit( n = c(5, 15), p = .5, pearly = .1, alpha = .1)  # Fail: one p

ph2crit( n = c(5, 15), p = c(.8, .2), pearly = -.5, alpha = .1) # Fail: pearly < 0

ph2crit( n = c(5, 15), p = c(.8, .2), pearly = 1.5, alpha = .1) # Fail: pearly > 1

ph2crit( n = c(5, 15), p = c(.8, .2), pearly = .1, alpha = -.2) # Fail: alpha < 0

ph2crit( n = c(5, 15), p = c(.8, .2), pearly = .1, alpha = 2.1) # Fail: pearly > 1



# Test with valid input parameters:

ph2crit( n = c( 5, 31), p = c(.8,.2), pearly = .1, alpha = .1) # Design A: 3 11

ph2crit( n = c( 8, 28), p = c(.8,.2), pearly = .1, alpha = .1) # Design B: 5 11

ph2crit( n = c(11, 25), p = c(.8,.2), pearly = .1, alpha = .1) # Design C: 7 11

ph2crit( n = c(12, 24), p = c(.8,.2), pearly = .1, alpha = .1) # Design D: 8 11








###############################################################





ph2dspb <- function(p, s, t)
  
  
  
  # Stopped negative binomial distribution mass function
  
  # Returns the probability density function for the
  
  #   minimum number of events for either s successes
  
  #   or else t failures
  


{
  
  if( length(p) != 1 || p < 0 || p > 1 ||
        
        s < 1 || t < 1 ||
        
        !is.wholenumber(s) || !is.wholenumber(t))
    
  {
    
    warning("Invalid parameter values (PhII)")
    
    return(NaN)
    
  }
  
  
  
  tmnb <- rep(0, s + t - 1)    # zero out, over the range
  
  for(j in 0 : (t - 1))        # Last event is success:  eq # (6)
    
  {
    
    tmnb[j + s] <- choose(s + j - 1, s - 1) * p ^ s * (1 - p) ^ j
    
  }
  
  
  
  for(j in 0 : (s - 1))       # Last event is failure:  eq # (7)
    
  {
    
    tmnb[j + t] <- tmnb[j + t] +
      
      choose(t + j - 1, t - 1) * p ^ j * (1 - p) ^ t
    
  }
  
  
  
  return(tmnb / sum(tmnb))    # Normalize
  
}



# # # # # # # # # #  Test error checking:



ph2dspb(p=c(.8, .2), s = 3, t = 4) # Invalid length(p) != 1

ph2dspb(p = -.2, s = 1, t = 3)     # Invalid: p < 0

ph2dspb(p = 1.2, s = 1, t = 3)     # Invalid: p > 1

ph2dspb(p = .8, s = 0, t = 3)      # Invalid: s < 1

ph2dspb(p = .8, s = 1.5, t = 3)    # Invalid: s not integer

ph2dspb(p = .8, s = 1, t = -3)     # Invalid: t < 1

ph2dspb(p = .8, s = 1, t = 3.001)  # Invalid: t not integer



ph2dspb(p = .8, s = 3, t = 5)      # Valid example





###################################################################


# # # # # # # # # #  Test error checking:


ph2tnb(p=-.2, n = 6, r = 4)    # Invalid: p < 0

ph2tnb(p=1.2, n = 6, r = 4)     # Invalid: p > 1

ph2tnb(p=.8, n = 6, r = 0)      # Invalid: r < 1

ph2tnb(p=.8, n = 6, r = 1.5)    # Invalid: r not integer

ph2tnb(p = .8, n=0, r = 4)  # Invalid: n < 1


ph2tnb(p = .8, n = 6, r = 4)  # Valid example P(Y1 = y1 | continue)
ph2tnb(p = .2, n = 6, r = 3)  # Valid example P(Y1 = y1 | stop)



########################################################################


ph2Eearly <- function(p, n, r)
  
{
  
  # Mean and SD of the minimum number of Stage 1 patients necessary to be able
  
  # to decide whether to either continue to Stage 2 or else terminate early
  
  
  
  if(length(p) > 1){ p1 <- p[1] # p can be vector or scalar
  
  }else p1 <- p
  
  if(length(n) > 1){ n1 <- n[1] # n can be vector or scalar
  
  }else n1 <- n
  
  if(length(r) > 1){ r1 <- r[1] # r can be vector or scalar
  
  }else r1 <- r
  
  
  
  tt <- n1 - r1 + 1         # Number of failures to stop early
  
  e <- sum((1 : n1) * ph2dspb(p1, r1, tt))              # expected value
  
  v <- sum((1 : n1) ^ 2 * ph2dspb(p1, r1, tt)) - e ^ 2  # variance
  
  c(e, sqrt(v))                                         # Expected value and SD
  
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # #



ph2Eearly( p = .8, n = 5, r = 3)   # Design A:   3.6336000 0.7344052

ph2Eearly( p = .8, n = 8, r = 5)   # Design B:   6.1063680 0.9998149

ph2Eearly( p = .8, n = 11, r = 7)  # Design C:   8.600271 1.230367

ph2Eearly( p = .8, n = 12, r = 8)  # Design D:   9.759577 1.261569

ph2Eearly( p = .2, n = 36, r = 11) # Binomial design X (p=.2): 31.943859  2.414899




#######################################################





ph2Ess <- function( p, n, r)
  
  # Expected curtailed sample size
  
{
  
  if( ! ph2valid(p, n, r))
    
  {                    # check parameter values
    
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
  

ph2Ess(p = c(.8, .2), n=c(6, 30), r=c(4, 11))       # Design A

ph2Ess(p = c(.8, .2), n = c( 8, 28), r = c(5, 11))  # Design B

ph2Ess(p = c(.8, .2), n = c(11, 25), r = c(7, 11))  # Design C

ph2Ess(p = c(.8, .2), n = c(12, 24), r = c(8, 11))  # Design D



ph2Ess(p = c(.8, .2), n = c(18, 18), r = c(12, 11)) # Balanced Design

ph2Ess(p = c(.8, .2), n = c(27,  9), r = c(19, 11))

#####################################################################



ph2mmax <- function(p, n, r)
  #  Evaluate the probability that the maximum sample size is needed
{
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

ph2mmax( p = c(.8, .2), n = c(3,33), r = ph2crit(n=c(3,33), p=c(.8, .2), pearly = .1, alpha =.1))  # one example






########################################################################
########################################################################
#### Build Tables and Figures
########################################################################
########################################################################




#  Build Table 1



p0 <- c(.8, .2)      # Null nypothesis

pA <- c(.8, .4)      # Alternative hypothesis

Des <- matrix( c(6, 30,   9, 27,   12, 24,   15, 21,  18, 18,   0, 36), 2,6)

(Des <- t(Des))      # Table of n values

Table2 <- NULL

ndes <- dim(Des)[1]

for (j in 1 : ndes)           #  loop over all designs (n)
  
{
  
  n <- Des[j, ]
  
  r <- ph2crit(pearly = .1, alpha = .1,  n, p = p0 )
  
  dr <- ph2reject(p = p0, n, r)
  
  de <- ph2early(p = p0, n, r)
  
  dss <- n[1] + (1 - de) * n[2]
  
  dp <- ph2reject(p = pA, n, r)
  
  line <- c(n, r, dr, dss, de, dp) # One line in Table 2
  
  Table2 <- rbind(Table2, line)
  
}



colnames(Table2) <- c("n1", "n2", "r1", "r2", "alpha",
                      
                      "E(SS)", "PrEarly", "Power")

rownames(Table2) <- c("A: Minimax/Optimal", "B", "C", "D", "E: Balanced", "X:Bin")

print(Table2, digits = 2)



# Traditional Simon 2-stage designs: for p1 = p2 = .2



library(clinfun)

Simon <- ph2simon(.2, .4, .1, .1)$out[2 : 1, ]      # Design parameters

Sa1 <- oc.twostage.bdry(.2, .4, Simon[1, 1], Simon[1, 2],
                        
                        Simon[1, 3], Simon[1, 4])              # probs, power, exp'd sample sizes

nS1 <- c(Simon[1, 2], Simon[1, 4] - Simon[1, 2])    # sample sizes: n's

rS1 <- c(Simon[1, 1], Simon[1, 3])                  # critical values: r's

line <- as.vector(c(nS1, rS1, Sa1[c(1, 4, 3, 2)]))  # line in Table 2

Table2 <- rbind(Table2,line)



Sa2 <- oc.twostage.bdry(.2, .4, Simon[2, 1], Simon[2, 2],
                        
                        Simon[2, 3], Simon[2, 4])   # Simon minimax design

nS2 <- c(Simon[2, 2], Simon[2, 4] - Simon[2, 2])

rS2 <- c(Simon[2, 1], Simon[2, 3])

line <- as.vector(c(nS2, rS2, Sa2[c(1, 4, 3, 2)]))

Table2 <- rbind(Table2,line)

rownames(Table2)[6 : 7] <- c("Y:So", "Z:Sm")

print(Table2, digits = 3)


######################################################



## Table 2




nA <- c( 6, 30) ; rA <- c(4, 11)    # Design parameters

nB <- c( 9, 27) ; rB <- c(6, 11)

nC <- c(12, 24) ; rC <- c(8, 11)

nD <- c(15, 21) ; rD <- c(10, 11)

nE <- c(18, 18) ; rE <- c(12, 11)

nBin <- c(0, 36) ; rBin <- c(0, 11) # Binomial design


pe <- ph2Eearly( p = .8, n = nA, r = rA)   # Design A

pc <- ph2Ess(p = c(.8, .2), n = nA, r = rA )

pa <- ph2Ess(p = c(.8, .4), n = nA, r = rA )

print(c(nA, rA, pe, pc, pa), digits = 2)

Table3 <- c(nA, rA, pe, pc, pa)



pe <- ph2Eearly( p = .8, n = nB, r = rB)   # Design B

pc <- ph2Ess(p = c(.8, .2), n = nB, r = rB)   # Design B

pa <- ph2Ess(p = c(.8, .4), n = nB, r = rB)   # Design B

print(c(pe, pc, pa), digits = 2)

Table3 <- rbind(Table3, c(nB, rB, pe, pc, pa))



pe <- ph2Eearly( p = .8, n = nC, r = rC)  # Design C

pc <- ph2Ess(p = c(.8, .2), n = nC, r = rC)  # Design C

pa <- ph2Ess(p = c(.8, .4), n = nC, r = rC)  # Design C

print(c(pe, pc, pa), digits = 2)

Table3 <- rbind(Table3, c(nC, rC, pe, pc, pa))



pe <- ph2Eearly( p = .8, n = nD, r = rD)  # Design D

pc <- ph2Ess(p = c(.8, .2), n = nD, r = rD) 

pa <- ph2Ess(p = c(.8, .4), n = nD, r = rD) 

print(c(pe, pc, pa), digits = 2)

Table3 <- rbind(Table3, c(nD, rD, pe, pc, pa))



pe <- ph2Eearly( p = .8, n = nE, r = rE)  # Design E

pc <- ph2Ess(p = c(.8, .2), n = nE, r = rE) 

pa <- ph2Ess(p = c(.8, .4), n = nE, r = rE) 

print(c(pe, pc, pa), digits = 2)

Table3 <- rbind(Table3, c(nE, rE, pe, pc, pa))




pc <- ph2Ess(p = c(.8, .2), n = nBin, r = rBin)  # Binomial

pa <- ph2Ess(p = c(.8, .4), n = nBin, r = rBin)  # Binomial

print(c(pc, pa), digits = 3)

Table3 <- rbind(Table3, c(nBin, rBin, 0, 0, pc, pa))



row.names(Table3) <- c("A", "B", "C", "D", "E", "Bin")

colnames(Table3) <- c("n1", "n2", "r1", "r2",
                      
                      "ES1", "SDS1", "EssH0", "EssHA")

Table3 <- data.frame(Table3)

print(Table3, digits=3)







#######################################################



#  Power comparisons of all designs for Figure 2



p <- c(.8, .2)

alpha <- .1

pearly <- .1

Des <- matrix( c(nA,  nB,  nC,  nD, nE,  nBin), 2,6)

(Des <- t(Des))      # Table of n values



library(clinfun)

Simon <- ph2simon(.2, .4, .1, .1)$out[2 : 1 , ]



ph2power <- function(design)
  
{
  
  nn <- dim(design)[1]         # number of designs
  
  altlist <- seq(p[2], .4, length.out = 12) # alternative list
  
  malt <- matrix(altlist, 1, length(altlist))
  
  power <- data.frame(malt)    # first row are alternatives
  
  
  
  for (i in 1 : (nn + 2))          # different designs + 2 Simon
    
  {
    
    if(i <= nn)n <- design[i, ]
    
    r <- ph2crit(n, p, alpha, pearly)
    
    ###        print(c(n, r))
    
    thisp <- NULL            # Start with empty output row
    
    for (alt in altlist)     # different alternatives
      
    {
      
      ppow <- c(p[1], alt)
      
      lp <- ph2reject(ppow, n, r) # power here
      
      if(i > nn)          # Simon 2-stage designs
        
      {
        
        sp <- as.vector(Simon[i - nn, 1 : 4])
        
        lp <- oc.twostage.bdry(ppow[1], ppow[2],
                               
                               sp[1], sp[2], sp[3], sp[4])[2]
        
      }
      
      thisp <- c(thisp, lp) # add to the row
      
    }
    
    
    
    ###        print(thisp)
    
    power <- rbind(power, thisp) # add row to data.frame
    
  }
  
  return(power)
  
}



power <- ph2power(Des)

print(power, digits = 2)





############################################################





#  Figure 3: Power curves



### pdf(file = "power.pdf")



plot(range(power[1, ]), c(.05, .95), type = "n",
     
     ylab = "Power", xlab = expression(p[2]), cex.lab = .8, cex.axis=.8)

for (i in 2 : (dim(power)[1] - 1))    # plot all designs
  
{
  
  lines(power[1, ], power[i, ], type = "l")
  
}



#text(.38, c( .69, .94),
     
#     labels = c("Designs A-E", "Designs X-Z"), cex=.9)

text(.37, .55, "Designs from top to bottom:", cex=.9)
#text(.4, c(.5, .45, .4, .3, .25, .2, .15, .1), c("X", "Y", "Z", "A", "B", "C", "D", "E"), cex=.9)
text(.4, seq(.5, .44, length.out=3), c("X", "Y", "Z"))
text(.4, seq(.38, .26, length.out=5), c("A", "B", "C", "D", "E"))
text(.42, .9, "X", cex=.9)
### dev.off()





###############################################################



# Figure 4: Contour plot of power in p1 and p2



mo <- function(x,y)   #  my "outer()" function
  
{
  
  mo <- matrix(0, length(x), length(y))
  
  for (i in 1 : length(x))
    
  {
    
    for (j in 1 : length(y))
      
    {
      
      mo[i,j] <- zfun( x[i], y[j])
      
    }}
  
  return(mo)
  
}



nF <- c(6, 30)  # Suggested design

rF <- c(4, 11)



zfun <- function(p1, p2)
  
{ ph2reject(c(p2, p1), nF, rF) }



p1 <- seq(from = .5, to = .99, length.out = 30)

p2 <- seq(from = .18, to = .5, length.out = 30)



### pdf(file="powerF.pdf")



z <- mo (p2, p1)

contour(p2, p1, z, xlab = expression(p[2]),
        
        ylab = expression(p[1]), cex.axis=.8, cex.lab=.8, ylim=c(.5, .9))

points(.2, .8, pch = 16, cex = 1.3)
text(.26, .82, label = "Null hypothesis", pos = 1, cex=.9)



### dev.off()



##############################################################



ph2all <- function( p, ns)
  
{
  
  #  systematically try ALL possible designs with a sample of size ns
  
  Des <- NULL
  
  for (i in c(1 : (ns-1)))
    
  {
    
    n <- c( i, ns - i)
    
    r <- ph2crit( n, p, pearly = .1, alpha = .1)
    
    ess <- ph2Ess( p, n, r)
    
    
    
    Des <- rbind( Des, c(n, r, ess))
    
  }
  
  colnames(Des) <- c("n1", "n2", "r1", "r2", "Ess")
  
  Des
  
}



ph2all(p = c( .8, .2), ns = 36)


#######################################################




## Figure 5:  Minimax Design
#  Look at all designs with 36 patients.
p <- c(.8, .2)

mmax <- rep(NA, 35)
for (i in 1:35)
{
  nn <- c(i, 36 - i)
  r <- ph2crit(n=nn, p=p, pearly = .1, alpha =.1)
  prob <- ph2mmax(n = nn, p=p, r)
  if(length(prob)!=0){
    mmax[i] <- prob
  }
}
mmax <- cbind(1:35, mmax)

plot(mmax[,1], mmax[,2], type='l', ylab="Probability of maximum sample size", 
     xlab=expression(n[1]),cex.axis=.8, cex.lab=.8 )





##############################################################

#######################  end of this file ####################

##############################################################
