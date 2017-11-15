#' Generic function for creating single-stage curtail trials
#' @export
setGeneric("two_stage_curtail_trial", function(p_null, p_alt, n, n_total, r,
                                               prob_early, alpha) {
  standardGeneric("two_stage_curtail_trial")
})

#' @examples
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30), r=c(4, 11))
#' @export
setMethod("two_stage_curtail_trial",
          signature(p_null="numeric", p_alt="numeric", n="numeric", 
                    n_total="missing",r="numeric", alpha="missing", 
                    prob_early="missing"),
          function(p_null, p_alt, n, r) {
            two_stage_valid_parameters_n_r(p_null, p_alt, n, r)
            Power <- Significance <- ess <- NULL
            ret <- data.frame(p_null.1=p_null[1], p_null.2=p_null[2],
                              p_alt.1=p_alt[1],p_alt.2=p_alt[2], n.1=n[1], 
                              n.2=n[2], r.1=r[1], r.2=r[2])
            class(ret) <- c("two_stage_curtail_trial", class(ret))
            ret
          })


#' @examples
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30), alpha=0.1, prob_early=0.1)
#' @export
setMethod("two_stage_curtail_trial",
          signature(p_null="numeric", p_alt="numeric", n="numeric", 
                    n_total="missing", r="missing", alpha="numeric", 
                    prob_early="numeric"),
          function(p_null, p_alt, n, prob_early, alpha) {
            two_stage_valid_parameters_n_alpha(p_null, p_alt, n,
                                                  prob_early, alpha)
            
            r <- Power <- Significance <- ess <- NULL
            r <- critical_values(n, p_null, prob_early, alpha)
            ret <- data.frame(p_null.1=p_null[1], p_null.2=p_null[2],
                              p_alt.1=p_alt[1],p_alt.2=p_alt[2], n.1=n[1], 
                              n.2=n[2], r.1=r[1], r.2=r[2])
            class(ret) <- c("two_stage_curtail_trial", class(ret))
            ret
          })

#' @examples
#' trials <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n_total=36, alpha=0.1, prob_early=0.1)
#' @export
setMethod("two_stage_curtail_trial",
          signature(p_null="numeric", p_alt="numeric", n="missing", 
                    n_total="numeric", r="missing", alpha="numeric", 
                    prob_early="numeric"),
          function(p_null, p_alt, prob_early, alpha, n_total) {
            two_stage_valid_parameters_n_total_alpha(p_null, p_alt, n_total,
                                                  prob_early, alpha)
            
            n1 <- Power <- Significance <- ess <- n <- r <- NULL
            n1 <- seq_len(n_total-1)
            n2 <- n_total - n1
            n <- cbind(n1, n2)
            r <- matrix(apply(n, 1, function(x){
              critical_values(x, p_null,prob_early, alpha)
            }), nrow=length(n1), ncol=2, byrow=TRUE)

            ret <- data.frame(p_null=matrix(rep(p_null,length(n1)), ncol=2, 
                  byrow=TRUE), p_alt=matrix(rep(p_alt,length(n1)), ncol=2, 
                  byrow=TRUE), n.1 = n1, n.2 = n2, r=r)
            ret <- ret[-which(ret$r.1==0),]
            class(ret) <- c("two_stage_curtail_trial_sawtooth", 
                            "two_stage_curtail_trial", class(ret))
            ret
})

#' @import ggplot2
#' @importFrom reshape melt
#' @export
#' 
plot.two_stage_curtail_trial_sawtooth <- function(x, ...) {
  x["Optimal Criteria"] <- sample_size(x)
  x["Minimax Criteria"] <- minimax_probability(x)

  m <- melt(subset(x, select=c("n.1", "Optimal Criteria", "Minimax Criteria")), id.var="n.1")
  ggplot(m, aes(x = n.1, y = value)) + geom_line() + 
    facet_grid(variable ~ ., scales = "free_y") + xlab(expression(n[1]))+ 
    ylab("")
}

#' @export
significance.two_stage_curtail_trial <- function(x) {
  apply(x, 1, function(x){
    p <- x[c("p_null.1", "p_null.2")]
    n <- x[c("n.1", "n.2")]
    r <- x[c("r.1", "r.2")]
    
    reject <- 0 
    
    # Loop on Y1 = curtailed sample size in Stage 1
    
    for (y1 in r[1] : n[1]) {    #  First summation 
      t1 <- (1-prob_early_stop(p, n, r))*ph2tnb(p, n, r)[y1]
      
      # Loop on X12 = Stage 1 successes who become Stage 2 successes
      lowx12 <- max(0, r[2] - (sum(n)-y1))  # lower limit for X12
      
      for (x12 in lowx12 : r[1]) {      
        # X12 conditional on X1
        
        t2 <- dbinom(x12, r[1], p[2]/p[1])  #  Prob of X12 given Y1
        # Loop on X2 <- Stage 2 successes - Third summation
        lowx2 <- max(0, r[2] - x12) # lower limit for X2
        
        for (x2 in lowx2:(sum(n) - y1)) {
          reject <- reject + t1 * t2 * dbinom(x2, sum(n)-y1, p[2])    
        }
      }
    }
    return(reject)
    })
  
}

#' @export
power.two_stage_curtail_trial <- function(x) {
  apply(x, 1,
        function(x) {
          p <- x[c("p_alt.1", "p_alt.2")]
          n <- x[c("n.1", "n.2")]
          r <- x[c("r.1", "r.2")]
          
          reject <- 0 
          
          # Loop on Y1 = curtailed sample size in Stage 1
          
          for (y1 in r[1] : n[1]) {    #  First summation 
            t1 <- (1-prob_early_stop(p, n, r))*ph2tnb(p, n, r)[y1]
            
            # Loop on X12 = Stage 1 successes who become Stage 2 successes
            lowx12 <- max(0, r[2] - (sum(n)-y1))  # lower limit for X12
            
            for (x12 in lowx12 : r[1]) {      
              # X12 conditional on X1
              
              t2 <- dbinom(x12, r[1], p[2]/p[1])  #  Prob of X12 given Y1
              # Loop on X2 <- Stage 2 successes - Third summation
              lowx2 <- max(0, r[2] - x12) # lower limit for X2
              
              for (x2 in lowx2:(sum(n) - y1)) {
                reject <- reject + t1 * t2 * dbinom(x2, sum(n)-y1, p[2])    
              }
            }
          }
          return(reject)
        })
}

#' @export
sample_size.two_stage_curtail_trial <- function(x) {
  
  apply(x, 1,
        function(x) {
          expected_total_sample_size(p = c(x["p_null.1"], x["p_null.2"]), 
                                     n = c(x["n.1"], x["n.2"]), r=c(x["r.1"], x["r.2"]))
        })
}

#' @export
stage1_sample_size <- function(x) {
  UseMethod("stage1_sample_size")
}

#' @export
stage1_sample_size.two_stage_curtail_trial <- function(x) {
  
  apply(x, 1,
        function(x) {
          
          n1 <- x["n.1"]
          r1 <- x["r.1"]
          p1 <- x["p_null.1"]
          
          tt <- n1 - r1 + 1  # Number of failures to stop early
          
          e <- esnb(p1, r1, tt)  # Expected value
          
          v <- vsnb(p1, r1, tt)   # variance
          
          eAndSD<-list("Expectation"=e, "Standard_Deviation"=sqrt(v))
          
          eAndSD  # Expected value and SD
          
        })
}

#' @export
PET <- function(x) {
  UseMethod("PET")
}

PET.two_stage_curtail_trial <- function(x){
  
  apply(x, 1,
        function(x) {
          j <- min(x[c("r.1", "n.1")]) : x["n.1"] # range of X1
          
          1 - sum(dbinom(j, x["n.1"], x["p_null.1"]))
        }) 
}

#' @export
minimax_probability <- function(x) {
  UseMethod("minimax_probability")
}

minimax_probability.two_stage_curtail_trial <- function(x){
  
  apply(x, 1,
        function(x) {
          n <- c(x["n.1"], x["n.2"])
          p <- c(x["p_null.1"], x["p_null.2"])
          r <- c(x["r.1"], x["r.2"])
          if (r[1]==0) {
            return(choose(n[1]+n[2]-1, r[2]-1)*p[2]^(r[2]-1)*
                     (1-p[2])^(n[1]+n[2]-r[2]))
          }
          else{
            # Probability of stopping early
            pearly <- prob_early_stop(p, n, r) 
            minimaxDesign<- 0
            
            # Distribution of Y1  |  Don't stop early
            ck <- rep(0, n[1])     
            k <- 1 : n[1]
            
            ck <- p[1] ^ r[1] * (1 - p[1]) ^ (k - r[1])* 
              choose(k - 1, r[1] - 1)
            ck <- ck / sum(ck[r[1]:n[1]]) 
            
            for (j in r[1] : n[1]) {
              #  Pr[ Y1 = j | don't stop early]
              y1j <- ck[j]           
              conv <- 0
              # Y2 is convolution sum of two binomials
              for (i in 0 : r[1]) {    
                x12 <- dbinom(i, r[1], p[2] / p[1])
                xp2 <- dbinom(r[2] - i - 1, sum(n) - j - 1, p[2])
                conv <- conv + x12 * xp2
              }
              minimaxDesign<- minimaxDesign+ conv * y1j
            }
            minimaxDesign* (1 - pearly)
          }
        })
}

#' @export
minimax_design <- function(x) {
  UseMethod("minimax_design")
}


#' @export
minimax_design.two_stage_curtail_trial_sawtooth <- function(x){
  x$minimax_probability <- minimax_probability(x)
  return(x[which.min(x$minimax_probability),])
  
}


#' @export
optimal_design <- function(x) {
  UseMethod("optimal_design")
}

#' @export
optimal_design.two_stage_curtail_trial_sawtooth <- function(x){
  x$sample_size <- sample_size(x)
  x[which.min(x$sample_size),]
}

summary.two_stage_curtail_trial_sawtooth <- function(x){
  opt_design <- optimal_design(x)
  class(opt_design) <- c("two_stage_curtail_trial", "data.frame")
  min_design <- minimax_design(x)
  class(min_design) <- c("two_stage_curtail_trial", "data.frame")
  list("Optimal_Design" = opt_design, "Minimax_Design" = min_design)
}

summary.two_stage_curtail_trial <- function(x){
  x$Power <- power(x)
  x$Significance <- significance(x)
  x$PET <- PET(x)
  x$Stage1_ESS <- stage1_sample_size(x)[[1]]$Expectation
  x$ESS <- sample_size(x)
  x$Minimax_Probability <- minimax_probability(x)
  x
}

##################################################################
two_stage_valid_parameters_n_r <- function(p_null, p_alt, n, r){
  
  if(any(p_null > 1) | any(p_null < 0)) stop("p_null should be between 
                                             0 and 1") 
  if(any(p_alt > 1) | any(p_alt < 0)) stop("p_alt should be between 
                                           0 and 1") 
  if(any(!is.wholenumber(n))) stop("n1 and n2 should be whole numbers") 
  if(any(n < 1)) stop("n1 and n2 should be at least 1") 
  if(any(r < 0)) stop("r1 and r2 should not be less than 0")
  if( r[1] > n[1]) stop("r1 should be less than or equal to n1")
  # Must have: 0 <= r2 <= n1+n2
  if(r[2] > sum(n)) stop("r2 should be less than or equal to n1+n2")
  if(any(!is.wholenumber(r))) stop("r1 and r2 should be whole numbers")
}

two_stage_valid_parameters_n_alpha <- function(p_null, p_alt, n, 
                                          prob_early, alpha){
                                                       
  if(any(p_null > 1) | any(p_null < 0)) stop("p_null should be between 
                                             0 and 1") 
  if(any(p_alt > 1) | any(p_alt < 0)) stop("p_alt should be between 
                                           0 and 1") 
  if(any(!is.wholenumber(n))) stop("n1 and n2 should be whole numbers") 
  if(any(n < 1)) stop("n1 and n2 should be at least 1") 
  if(alpha > 1 | alpha < 0) stop("alpha should be between 
                                 0 and 1")
  if(prob_early > 1 | prob_early < 0) stop("prob_early should be between 
                                           0 and 1")
}

two_stage_valid_parameters_n_total_alpha <- function(p_null, p_alt, n_total, 
                                                     prob_early, alpha){
  
  if( any(p_null > 1) | any(p_null < 0)) stop("p_null should be between 
                                              0 and 1") 
  if( any(p_alt > 1) | any(p_alt < 0)) stop("p_alt should be between 
                                            0 and 1") 
  if(!is.wholenumber(n_total)) stop("n_total should be a whole number") 
  if(n_total < 2) stop("n_total should be at least 2")
  if(alpha > 1 | alpha < 0) stop("alpha should be between 0 and 1")
  if(prob_early > 1 | prob_early < 0) stop("prob_early should be between
                                           0 and 1")
}

prob_early_stop <- function(p, n, r){
  j <- min(r[1], n[1]) : n[1]    # range of X1
  
  return(1 - sum(dbinom(j, n[1], p[1])))
}


#' Test for integer
#' 
#' @param x the number to test.
#' @param tol the tolerance for x being a whole number (default 
#' \code{.Machine$double.eps ^ 0.5})
is.wholenumber = function(x, tol = .Machine$double.eps ^ 0.5){
  abs(x - round(x)) < tol
}

#' Truncated negative binomial distribution mass function
#'
#' Returns the probability density function for the number of events
#' to reach r= r1 successes in Stage 1 with p = p1 probability of success.
#' Can also be the number of events to reach r failures with p <- (1-p1) 
#' probability of failure.
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) or a scalar containing p1
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) or a scalar containing n1
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2) or a scalar containing r1
#' @examples
#' # ph2tnb(p=0.8, n=6, r=4) 
#' # ph2tnb(p=0.2, n=6, r=3) 
ph2tnb = function(p, n, r) {
  if(length(p) > 1) { 
    p1 <- p[1] # p can be vector or scalar
  } else {
    p1 <- p
  }
  
  if(length(n) > 1) { 
    n1 <- n[1] # n can be vector or scalar
  } else {
    n1 <- n
  }
  
  if(length(r) > 1) { 
    r1 <- r[1] # r can be vector or scalar
  } else {
    r1 <- r
  }
  
  if(p1 < 0 || p1 > 1 ||
     r1 < 1 || n1 < 1 ||
     !is.wholenumber(r1)) {
    warning("Invalid parameter values")
    return(NaN)
  }
  tnb <- NULL
  for (y1 in 1:n1) {
    if (y1-r1 >= 0) {
      num <- factorial(y1-1)/factorial(y1-r1)
    } else {
      num <- 0
    }
    denom <- 0 
    for (j in r1:n1) {
      a <- factorial(j-1)/factorial(j-r1) *(1-p1)^(j-y1)
      denom <- denom + a
    }
    tnb <- c(tnb, (num/denom))
  }
  return(tnb)
}

