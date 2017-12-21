
#' Calculate significance for the single stage trial
#' 
#' Calculates the probability of rejecting the null hypothesis assuming the 
#' null probability of success
#'
#' @param pNull probability of success under the null hypothesis
#' @param s number of successes to stop the trial
#' @param t number of failures to stop the trial
#' @examples
#' single_stage_significance(0.2, 7, 11)
single_stage_significance <- function(pNull, s, t) {
  sum(dsnb_stacked(min(s, t):(s+t-1), p=pNull, s=s, t=t)[,'s'])
}

#' Calculate power for the single stage trial
#' 
#' Calculates the probability of rejecting the null hypothesis assuming an 
#' alternative probability of success
#'
#' @param pAlt probability of success under the alternative hypothesis
#' @param s number of successes to stop the trial
#' @param t number of failures to stop the trial
#' @examples
#' single_stage_power(0.5, 7, 11)
single_stage_power <- function(pAlt, s, t) {
  sum(dsnb_stacked(min(s, t):(s+t-1), p=pAlt, s=s, t=t)[,'s'])
}


#' Expected sample size and SD for the one-stage design
#'
#' Mean and standard deviation of the minimum number of patients to make a 
#' decision about rejecting the null for the one-stage design.
#'
#' @param p scalar for the one-stage design containing the probability of 
#' successful outcome (p)
#' @param s number of successes to stop the trial
#' @param t number of failures to stop the trial
#' @examples
#' single_stage_expected_sample_size(p = 0.2, s = 7, t = 11)
single_stage_expected_sample_size <- function(p, s, t) {
  if(length(p) != length(s) | length(s) != length(t) | length(p) != length(t))
    stop(paste("Parameters p, s, and t must all be the same length",
               "(all length 1 for the one-stage design, or all length 2 for",
               "the two-stage design)"))
  e <- sum((1:(s+t-1)) * ph2snb(p, s, t))  # Expected value
  v <- sum((1:(s+t-1))^2 * ph2snb(p, s, t)) - e^2   # variance
  eAndSD <- list("expectation"=e, "standardDeviation"=sqrt(v))
  return(eAndSD)  # Expected value and SD
  
}

#' Plot Power and Significance for One-Stage Design
#'
#' Plot power and significance across all trial designs with a maximum number of patients, n,
#' with a varying number of responses, s, required to reach the success endpoint
#'
#' @param n maximum number of patients in the trial
#' @param pNull probability of success under the null hypothesis
#' @param pAlt probability of success under the alternative hypothesis
#' @import ggplot2
#' @importFrom foreach foreach %do%
#' @importFrom tidyr gather
#' @examples
#' power_significance_plot(17, 0.2, 0.4)
power_significance_plot <- function(n, pNull, pAlt){
  value <- Significance <- Power <- si <- s <- `Design Feature` <- NULL
  
  designs <- foreach (si=seq_len(n-1), .combine=rbind) %do% {
    ti <- n + 1 - si
    c(si, ti, single_stage_significance(pNull, si, ti), 
      single_stage_power(pAlt, si, ti), esnb(pNull, si, ti))
  }
  rownames(designs) <- NULL
  colnames(designs) <- c("s", "t", "Significance", "Power", "ess")
  designs <- as.data.frame(designs)
  
  designs_long <- gather(designs, `Design Feature`, value, Significance:Power)
  
  ggplot(designs_long, aes(x=s, y=value, color=`Design Feature`)) + 
    geom_line() + ylab("Probability") + 
    xlab("Number of Responses to Stop the Trial (s)") + theme_bw() +
    scale_x_continuous(breaks=seq_len(n-1))
  
  
}
#' ROC curve for Power vs. 1-Significance
#'
#' ROC curve of all trial designs with a maximum number of patients, n,
#' with a varying number of responses, s, required to reach the success endpoint
#'
#' @param n maximum number of patients in the trial
#' @param pNull probability of success under the null hypothesis
#' @param pAlt probability of success under the alternative hypothesis
#' @param all_labels controls how many labels of s will be included in the plot.  if set to TRUE,  all labels of s will appear 
#' on the plot. otherwise, if FALSE (default), the labels for the extreme values of s will be removed to improve the readability
#' @import ggplot2
#' @importFrom foreach foreach %do%
#' @examples
#' power_significance_ROC(17, 0.2, 0.4)
power_significance_ROC <- function(n, pNull, pAlt, all_labels=FALSE){
  si <- Power <- Significance <- s <- NULL
  designs <- foreach (si=seq_len(n-1), .combine=rbind) %do% {
    ti <- n + 1 - si
    c(si, ti, single_stage_significance(pNull, si, ti), single_stage_power(pAlt, si, ti),
      esnb(pNull, si, ti))
  }
  rownames(designs) <- NULL
  colnames(designs) <- c("s", "t", "Significance", "Power", "ess")
  designs <- as.data.frame(designs)
  
  data_pos <- designs[1:(n-1), c("Power", "Significance", "s")]
  data_pos$Power <- data_pos$Power + 0.02
  data_pos$Significance <- data_pos$Significance - 0.02
  # stop labeling when power < .05 or significance>.95 when all_labels=FALSE
  if(all_labels==FALSE)  
    data_pos[which(data_pos$Power<0.05 | data_pos$Significance>.95), "s"] <- ""
  ggplot(designs, aes(x=Power, y=1-Significance)) + 
    geom_line() +   
    geom_text(data=data_pos, aes(x=Power, y=1-Significance, label=s)) +
    theme_bw()
}



#' Find critical values for decision making during the one-stage or two-stage trial 
#' 
#' Finds the critical number of successes necessary to reject 
#' the null hypothesis in the one- and two-stage designs, 
#' as well as the number of Stage 1 successes neccesary to 
#' continue to Stage 2 in the two-stage design.
#' 
#' For the one-stage design:  
#' Finds the value of s, the minimum number of successes to reject the 
#' null hypothesis, from the Binomial model with significance level alpha.
#' 
#' For the two-stage design:
#' Finds r1, the minimum number of Stage 1 successes to continue to 
#' Stage 2, and r2, the minimum number of Stage 2 successes to 
#' reject the null hypothesis.  Find r1 so that the probability 
#' of early stopping is less than or equal to pearly. 
#' Find r2 from Binomial model with no early stopping
#' and significance level alpha.
#'
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'          and Stage 2 (n2).
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) under the null hypothesis.
#' @param pearly desired probability of early stopping (default is 0.1).  
#' Not necessary for determining critical values for the one-stage design.
#' @param alpha desired significance level (default is 0.1).
#' @importFrom stats qbinom
#' @examples
#' critical_values(n=36, p=0.2, alpha=.1)
#' critical_values(n=c( 5, 31), p=c(.8,.2), pearly=.1, alpha=.1)
critical_values <- function(n, p, pearly=.1, alpha=.1) {
  if(length(p) != length(n))
    stop('Parameters p and n must be the same length (both length 1 for the one-stage design, or both length 2 for the two-stage design)')
  
  # one-stage design
  if(length(p) == 1 && length(n) == 1){
    
    if (!ph2valid(p, n, r=0) ||
        pearly < 0 || pearly > 1 || alpha < 0 || alpha > 1)
    {
      warning("Invalid parameter values")
      return (NaN)
    }
    return(1 + qbinom(1 - alpha, n, p))  
  }
  
  
  ## two-stage design
  else{
    
    if (!ph2valid(p, n, r=c(0, 0)) ||
        pearly < 0 || pearly > 1 || alpha < 0 || alpha > 1)
    {
      warning("Invalid parameter values, pearly and alpha must be in [0, 1]")
      return (NaN)
    }
    # Binomial critical value, r2
    r2 <- 1 + qbinom(1 - alpha, sum(n), p[2])     
    
    for (j in 0 : n[1]){
      
      # r1 based on pearly
      if (prob_early_stop(p, n, r=c(j, 0)) > pearly)
        
        return(c(max(j - 1, 0), r2))
    }
    
    return (c(n[1], r2))        
  }
}

#' Probability of stopping early in the two-stage design
#'
#' Calculates the probability of stopping the trial after Stage 1
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'          and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2) 
#' @examples
#' prob_early_stop(p=c(0.8, 0.2), n=c(25, 9), r=c(17, 11))
prob_early_stop <- function(p, n, r) {
  if( ! ph2valid(p, n, r)){
    warning("Invalid parameter values")
    return(NaN)
  }
  
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
  
  1 - sum(dbinom(j, n1, p1))
}

#' Expected sample size and SD for Stage 1 of the two-stage design 
#'
#' Mean and standard deviation of the minimum number of patients to make a decision to either 
#' continue to Stage 2 or else terminate early. 
#'
#' @param p a vector containing the probability of successful outcomes in 
#'          in Stage 1 (p1) and Stage 2 (p2)
#' @param n a vector containing maximum sample sizes planned for Stage 1 (n1) and Stage 2 (n2)
#' @param r a vector containing the minimum number of Stage 1 successes to continue to 
#' Stage 2 (r1) and the minimum number of Stage 2 successes to reject the null hypothesis (r2) 
#' @examples
#' expected_stage1_sample_size(p=c(.8, .2), n=c(5, 31), r=c(3, 11))
expected_stage1_sample_size <- function(p, n, r) {
  if(length(p) != length(n) | length(n) != length(r) | length(p) != length(r))
    stop('Parameters p, n, and r must all be the same length (all length 1 for the one-stage design, or all length 2 for the two-stage design)')
  if(length(p) > 1){ p1 <- p[1] # p can be vector or scalar
  
  }else p1 <- p
  
  if(length(n) > 1){ n1 <- n[1] # n can be vector or scalar
  
  }else n1 <- n
  
  if(length(r) > 1){ r1 <- r[1] # r can be vector or scalar
  
  }else r1 <- r
  
  
  
  tt <- n1 - r1 + 1  # Number of failures to stop early
  
  e <- sum((1 : n1) * ph2snb(p1, r1, tt))  # Expected value
  
  v <- sum((1 : n1) ^ 2 * ph2snb(p1, r1, tt)) - e ^ 2   # variance
  
  eAndSD<-list("expectation"=e, "standardDeviation"=sqrt(v))
  
  return(eAndSD)  # Expected value and SD
  
}

#' Expected curtailed sample size for the two-stage design
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
#' expected_total_sample_size(p=c(0.8, 0.2), n=c(6, 30), r=c(4, 11))
#' expected_total_sample_size(p=c(0.8, 0.2), n=c(18, 18), r=c(12, 11))
#' @importFrom stats dbinom
expected_total_sample_size <- function(p, n, r) {
  
  if (length(p) != 2 || length(n) != 2 || length(r) != 2) {
    stop("Parameters must be vectors of length 2")
  }
  # check parameter values
  if (!ph2valid(p, n, r)) {
    warning("Invalid parameter values")
    return(NaN)
  }
  
  pearly <- prob_early_stop(p, n, r)  # Probability of stopping early
  
  #Pr(Y1=y1|continue to stage 2)
  # truncated negative binomial distribution
  # number of trials until we reach r1 successes 
  # r1 <= y1 <= n1
  prY1Continue <- function(y1, r1, n1, p1) {
    if( y1 < r1 || y1 > n1){    
      warning("Invalid parameter values")
      return(NaN)
    }
    num <- factorial(y1-1) / factorial(y1-r1)
    denom <- sum(sapply(r1:n1, 
                        function(j) {
                          factorial(j-1)/factorial(j-r1) *(1-p1)^(j-y1)
                        }))
    
    if (is.na(num/denom)) { 
      return(0)
    } else {
      return(num/denom)
    }
  }
  
  # Pr(Y1=y1|stop early)
  # truncated negative binomial distribution
  # number of trials until we reach n1-r1+1 failures
  # n1-r1+1 <= y1 <= n1
  
  prY1Stop <- function(y1, r1, n1, p1) {
    if (y1 < n1-r1+1 || y1 > n1) { 
      warning("Invalid parameter values")
      return(NaN)
    }
    num <- factorial(y1-1) / factorial(y1-1-n1+r1)
    
    denom <- sum(sapply((n1-r1+1):n1, 
                        function(j) {
                          factorial(j-1)/factorial(j-1-n1+r1) * p1^(j-y1)
                        }))
    
    if(is.na(num/denom)){ 
      return(0)
    } else {
      return(num/denom)
    }
  }
  
  ### expected value of Y1 | continue to Stage 2
  EY1Continue <- sum(sapply(r[1]:n[1], 
                            function(i) {
                              i * prY1Continue(i, r[1], n[1], p[1])
                            }))
  
  ### expected value of Y1|stop early
  EY1Stop <- sum(sapply((n[1]-r[1]+1):n[1], 
                        function(i) {
                          i*prY1Stop(i, r[1], n[1], p[1])
                        }))
  
  ### expected value of y2|continue to Stage 2
  # i = possible values of y1
  # j = possible values of x12
  py1 <- sapply(r[1]:n[1], 
                function(i) {
                  prY1Continue(i, r[1], n[1], p[1])
                })
  
  px12 <- sapply(0:r[1], 
                 function(j) {
                   dbinom(j, r[1], p[2]/p[1])
                 })
  
  #for the snb distribution,
  #s <- r[2]-j
  #t <- n[2]+n[1]-i-r[2]+j+1
  ey2 <- function(y1) {
    sapply(0:r[1], 
           function(j) {
             esnb(p[2], s=ifelse((r[2]-j)>0, (r[2]-j), 0), 
                  t=ifelse((n[2]+n[1]-y1-r[2]+j+1)>0, n[2]+n[1]-y1-r[2]+j+1, 0))
           })
  }
  
  a <- sapply(r[1]:n[1], function(y1) sum(px12*ey2(y1)))
  EY2Continue <- sum(py1*a)
  
  return((EY1Stop*pearly) + ((EY1Continue+EY2Continue)*(1-pearly))) #Eq. 7
}


#' Find the minimax design for a two-stage trial
#'
#' Computes the minimax probability (probability of requiring all n1 + n2 
#' patients to reach a statistical endpoint) for each combination of n1 and 
#' n2 for a given total n. The minimax design is the design in the first row 
#' of output with the smallest minimax probability.
#'
#' @param p vector containing the probability of successful outcomes
#'        in Stage 1 (p1) and Stage 2 (p2) 
#' @param ntot scalar containing the total sample size planned for the 
#'        trial (n1+n2) 
#' @param pearly desired probability of early stopping (default = .1).
#' @param alpha desired significance level (default = .1).
#' @examples 
#' all_minimax_designs(c(.8, .2), 36)
#' all_minimax_designs(c(.7, .3), 40, pearly = .08, alpha=.1)
all_minimax_designs <- function(p, ntot, pearly = .1, alpha = .1) {
  
  if (any(p > 1) || any(p < 0) || (ntot < 0) || (pearly < 0) || (pearly > 1) || 
      (alpha < 0) || (alpha > 1)) {
    warning("Invalid parameter values")
    return (NaN)
  }
  
  rcrit <- sapply(1:(ntot-1), 
                  function(n1) {
                    critical_values(n=c(n1, ntot-n1), p=p, pearly=pearly, alpha = alpha)
                  })
  
  n1 <- (1:(ntot-1))[which(rcrit[1,]>0)]
  prob <- sapply(n1, 
                 function(n1) {
                   minimax_design_old(p=p, n=c(n1, ntot-n1), 
                                      r=critical_values(n=c(n1, ntot-n1), p=p, pearly=pearly, 
                                                        alpha = alpha))
                 })
  
  mmax <- data.frame(n1, (ntot-n1), prob)
  names(mmax) <- c("n1", "n2", "Probability of Maximum Sample Size")
  
  return(mmax[order(mmax[,3]),])
}

#' Find the optimal design for a two-stage trial
#'
#' Computes the expected sample size under curtailed sampling for each combination of 
#' n1 and n2 for a given total n.  The optimal design is the design in the first row 
#' of output with the smallest expected sample size under curtailed sampling.
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param ntot scalar containing the total sample size planned for the trial (n1+n2) 
#' @param pearly desired probability of early stopping (default = .1).
#' @param alpha desired significance level (default = .1).
#' @examples 
#' all_optimal_designs(c(.8, .2), 36)
#' all_optimal_designs(c(.7, .3), 40, pearly = .08, alpha=.1)
all_optimal_designs <- function(p, ntot, pearly = .1, alpha = .1) {
  if (any(p>1) || any(p < 0) || (ntot<0) || (pearly < 0) || (pearly > 1) || 
      (alpha < 0) || (alpha > 1)) {
    warning("Invalid parameter values")
    return (NaN)
  }
  
  rcrit <- sapply(1:(ntot-1), 
                  function(n1) {
                    critical_values(n=c(n1, ntot-n1), p=p, pearly=pearly, alpha = alpha)
                  })
  n1 <- (1:(ntot-1))[which(rcrit[1,]>0)]
  ss <- sapply(n1, 
               function(n1) {
                 expected_total_sample_size(p=p, n=c(n1, ntot-n1), 
                                            r=critical_values(n=c(n1, ntot-n1), p=p, 
                                                              pearly=pearly, alpha=alpha))
               })
  
  ess <- data.frame(n1, (ntot-n1), ss)
  names(ess) <- c("n1", "n2", "Expected Sample Size")
  return(ess[order(ess[,3]),])
}


#' Finds the minimax and optimal designs for a two-stage trial
#' 
#' Returns the minimax and optimal designs and their associated statistical 
#' properties.
#'
#' @param p vector containing the probability of successful outcomes
#' in Stage 1 (p1) and Stage 2 (p2) 
#' @param ntot scalar containing the total sample size planned for the 
#' trial (n1+n2) 
#' @param pearly desired probability of early stopping (default = .1).
#' @param alpha desired significance level (default = .1).
#' @examples 
#' best_designs(c(0.8, 0.2), 36)
#' best_designs(c(0.7, 0.3), 40, pearly = 0.08, alpha=0.1)
#' best_designs(c(0.7, 0.3), 40, pearly = 0.08, alpha=0.1)$designs
best_designs <- function(p, ntot, pearly = 0.1, alpha = 0.1) {
  if (any(!is.numeric(p)) || any(p > 1) || any(p < 0) || !is.numeric(ntot) || 
      (ntot < 0) || !is.numeric(pearly) || (pearly < 0) || (pearly > 1) || 
      !is.numeric(alpha) || (alpha < 0) || (alpha > 1)) {
    warning("Invalid parameter values")
    return (NaN)
  }
  
  opt <- all_optimal_designs(p, ntot, pearly, alpha)
  mini <- all_minimax_designs(p, ntot, pearly, alpha)
  optimalDes <- opt[which.min(opt[,3]),]
  minimaxDes <- mini[which.min(mini[,3]),]
  nOpt <- as.vector(c(optimalDes[1,1], optimalDes[1,2]))
  nMini <- as.vector(c(minimaxDes[1,1], minimaxDes[1,2]))
  rOpt <- critical_values(c(optimalDes[1,1], optimalDes[1,2]), p, pearly, alpha)
  rMini <- critical_values(c(minimaxDes[1,1], minimaxDes[1,2]), p, pearly, alpha)
  
  optimalDesign <- cbind(p[1], nOpt[1], rOpt[1], p[2], nOpt[2], rOpt[2], 
                         alpha, prob_early_stop(p, nOpt, rOpt), 
                         optimalDes[1,3], minimax_design_old(p, nOpt, rOpt))
  minimaxDesign <- cbind(p[1], nMini[1], rMini[1], p[2], nMini[2], rMini[2], 
                         alpha, prob_early_stop(p, nMini, rMini), 
                         expected_total_sample_size(p, nMini, rMini), minimaxDes[1,3])
  designs <- data.frame(rbind(optimalDesign, minimaxDesign))
  rownames(designs) <- c("Optimal", "Minimax")
  colnames(designs) <- c("p1", "n1", "r1", "p2", "n2", "r2", "Alpha", "PET", 
                         "ECSS", "P(MaxSS)")
  ret <- list(designs=designs, p=p, ntot=ntot, pearly=pearly, alpha=alpha)
  class(ret) = c("ph2_design", class(ret))
  return(ret)
}

#' @importFrom ggplot2 ggplot aes geom_line facet_grid xlab
plot.ph2_design = function(x, ...) {
  n1 <- values <- NULL
  p <- x$p
  ntot <- x$ntot
  pearly <- x$pearly
  alpha <- x$alpha
  opt <- all_optimal_designs(p, ntot, pearly, alpha)
  opt <- opt[order(opt[,1]),]
  mini <- all_minimax_designs(p, ntot, pearly, alpha)
  mini <- mini[order(mini[,1]),]
  
  df <- data.frame(n1 = c(opt[,1], mini[,1]), values = c(opt[,3], mini[,3]), 
                   type=c(rep("Optimal Criteria", dim(opt)[1]), rep("Minimax Criteria", 
                                                                    dim(mini)[1])))
  ggplot(data=df, aes(x=n1, y=values)) + geom_line() + 
    facet_grid(type ~ ., scales="free")+ xlab(expression(n[1])) + ylab("")
}

#' Evaluate the probability that the maximum sample size is needed for the two-stage design
#'
#' Evaluate the minimax probability that the maximum sample size n1+n2 is
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
#' minimax_design_old(p=c(0.8, 0.2), n=c(3,33), r=critical_values(n=c(3,33), 
#'   p=c(0.8, 0.2), pearly=0.1, alpha=0.1))
minimax_design_old <- function(p, n, r) {
  if (!ph2valid(p, n, r)){
    warning("Invalid parameter values")
    return(NaN) # Check validity of parameters
  } 
  if (r[1]==0) {
    return(choose(n[1]+n[2]-1, r[2]-1)*p[2]^(r[2]-1)*(1-p[2])^(n[1]+n[2]-r[2]))
  }
  else{
    pearly <- prob_early_stop(p, n, r) # Probability of stopping early
    minimaxDesign<- 0
    
    ck <- rep(0, n[1])     # Distribution of Y1  |  Don't stop early
    k <- 1 : n[1]
    
    ck <- p[1] ^ r[1] * (1 - p[1]) ^ (k - r[1]) * choose(k - 1, r[1] - 1)
    ck <- ck / sum(ck[r[1]:n[1]]) 
    
    for (j in r[1] : n[1]) {
      y1j <- ck[j]           #  Pr[ Y1 = j | don't stop early]
      conv <- 0
      for (i in 0 : r[1]) {    # Y2 is convolution sum of two binomials
        x12 <- dbinom(i, r[1], p[2] / p[1])
        xp2 <- dbinom(r[2] - i - 1, sum(n) - j - 1, p[2])
        conv <- conv + x12 * xp2
      }
      minimaxDesign<- minimaxDesign+ conv * y1j
    }
    minimaxDesign* (1 - pearly)
  }
}

#' Significance of the two-stage design under curtailed sampling
#' 
#' Calculates the probability of rejecting the null hypothesis in the 
#' two-stage design under curtailed sampling, assuming null probabililties of 
#' success for Stage 1 and Stage 2
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) under the null hypothesis
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
#' @examples
#' two_stage_significance(p = c( 0.8, 0.2), n = c(12, 24), r = c(8, 11))
#' two_stage_significance(p = c( 0.8, 0.2), n = c(6, 30), r = c(4, 11))
two_stage_significance <- function(p, n, r) {
  return(prob_reject(p=p, n=n, r=r))
}

#' Power of the two-stage design under curtailed sampling
#' 
#' Calculates the probability of rejecting the null hypothesis in the 
#' two-stage design under curtailed sampling, assuming alternative 
#' probabililties of success for Stage 1 and Stage 2.
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) under the alternative hypothesis
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
#' @examples
#' two_stage_power(p = c(0.8, 0.4), n = c(12, 24), r = c(8, 11))
#' two_stage_power(p = c(0.8, 0.4), n = c(6, 30), r = c(4, 11))
two_stage_power <- function(p, n, r) {
  return(prob_reject(p=p, n=n, r=r))
}

#' Probability of rejecting the null hypothesis in the two-stage design
#' under curtailed sampling
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
prob_reject <- function(p, n, r) {
  # check validity of parameter values
  if (!ph2valid(p, n, r)) {
    warning("Invalid parameter values")
    return (NaN)
  }
  
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
}


#' Probability of rejecting the null hypothesis in the two-stage design
#' under traditional sampling
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param n vector containing sample sizes planned for Stage 1 (n1) 
#'           and Stage 2 (n2) 
#' @param r vector containing the minimum number of Stage 1 successes 
#' to continue to Stage 2 (r1) and the minimum number of Stage 2 
#' successes to reject the null hypothesis (r2)
prob_reject_traditional = function(p, n, r) {
  # check validity of parameter values
  
  if ( !ph2valid(p, n, r)) {
    warning("Invalid parameter values")
    return (NaN)
  }
  
  reject <- 0
  
  # Loop on X1 = number of Stage 1 successes
  
  lowx1 <- max(r[1], r[2] - n[2])    # lower limit for X1
  
  for (x1 in lowx1 : n[1]) {
    
    t1 <- dbinom(x1, n[1], p[1]) # binomial probability of X1
    
    # Loop on X12 = Stage 1 successes who become Stage 2 successes
    
    lowx12 <- max(0, r[2] - n[2])  # lower limit for X12           
    for (x12 in lowx12 : x1) {
      # X12 conditional on X1
      
      t2 <- dbinom(x12, x1, p[2] / p[1])  #  Prob of X12 given X1
      
      # Loop on X2 = Stage 2 successes
      lowx2 <- max(0, r[2] - x12) # lower limit for X2
      for (x2 in lowx2 : n[2]) {
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
#' # ph2snb(p = .8, s = 3, t = 5)
ph2snb = function(p, s, t) {
  if( length(p) != 1 || p < 0 || p > 1 || s < 1 || t < 1 ||
      !is.wholenumber(s) || !is.wholenumber(t)){
    warning("Invalid parameter values")
    return(NaN)
  }
  
  tmnb <- rep(0, s + t - 1)    # zero out, over the range
  for(j in 0 : (t - 1)) {      # Last event is success:  eq # (6)
    tmnb[j + s] <- choose(s+j-1, s-1) * p^s * (1-p)^j
  }
  
  for(j in 0 : (s - 1)) {       # Last event is failure:  eq # (7)
    tmnb[j+t] <- tmnb[j+t] +
      choose(t+j-1, t-1) * p^j * (1-p)^t
  }
  return(tmnb / sum(tmnb))    # Normalize
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
#' # ph2valid(p = c(.2, .1), n = c(10, 10), r = c(5, 5))
#' # ph2valid(p = c( .2, .3), n = c(10, 10), r = c(5, 5))
ph2valid <- function(p,n,r) {
  if ((length(p) == 1) && (length(n) == 1) && (length(r) == 1)) {
    if( p > 1 | p < 0) return(FALSE) 
    if(n < 0) return(FALSE)
    if(!is.wholenumber(n)) return(FALSE) 
    if( r < 0) return(FALSE)
    if( r > n) return(FALSE)
    if(!is.wholenumber(r)) return(FALSE) 
    return(TRUE)               # Valid parameter values for one-stage
  } else {
    #  Must have:  0 <= p2 <= p1 <= 1
    if( length(p) != 2 ) return(FALSE)
    if( p[1] > 1 | p[2] < 0 | p[1] < p[2] )return(FALSE)
    #  The two n's must be non-negative integers
    if( length(n) != 2) return(FALSE)
    if(min(n) < 0) return(FALSE)
    if(!is.wholenumber(n[1]) || !is.wholenumber(n[2])) return(FALSE)
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
}



#' Test for integer
#' 
#' @param x the number to test.
#' @param tol the tolerance for x being a whole number (default 
#' \code{.Machine$double.eps ^ 0.5})
is.wholenumber = function(x, tol = .Machine$double.eps ^ 0.5){
  abs(x - round(x)) < tol
}

