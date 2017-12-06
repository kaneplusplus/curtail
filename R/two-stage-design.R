#' Create a two stage curtailed trial design
#'
#' @description WHAT DOES THE FUNCTION DO?
#' @param p1_null TODO:write this
#' @param p2_null TODO:write this
#' @param p1_alt TODO:write this
#' @param p2_alt TODO:write this
#' @param n1 TODO:write this
#' @param n2 TODO:write this
#' @param n_total TODO:write this
#' @param r1 TODO:write this
#' @param r2 TODO:write this
#' @param prob_early TODO:write this
#' @param alpha TODO:write this
#' @details
#' WHAT ARE THE DIFFERENT WAYS OF SPECIFYING THIS TRIAL
#' @examples
#' # Case 1
#' trial <- two_stage_curtail_trial(p1_null = 0.8, p2_null = 0.2, 
#' p1_alt = 0.8, p2_alt = 0.4, n1 = 6, n2 = 30, r1 = 4, r2 = 11)
#' 
#' # Case 2
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30))
#'
#' # Case 3
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30), prob_early=0.1, alpha=0.1)
#' 
#' # Case 4
#' trials <- two_stage_curtail_trial(p1_null = 0.8, p2_null=0.2, 
#' p1_alt = 0.8, p2_alt = 0.4, n_total=36, prob_early=0.1, alpha=0.1)
#' 
#' # Case 5
#' trials <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n_total=36)
#' @export
setGeneric("two_stage_curtail_trial", 
           function(p1_null, p2_null, p1_alt, p2_alt, n1, n2, n_total, 
                    r1, r2, prob_early, alpha) {
  standardGeneric("two_stage_curtail_trial")
})


#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="numeric", n2="numeric", n_total="missing",
                    r1="numeric", r2="numeric", 
                    prob_early="missing", alpha="missing"),
          function(p1_null, p2_null, p1_alt, p2_alt, 
                   n1, n2, r1, r2) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null), 
                          p_alt=c(p1_alt, p2_alt), n=c(n1, n2), 
                          n_total=NULL, r=c(r1, r2), prob_early=NULL, 
                          alpha=NULL)
            ret <- data.frame(p1_null=p1_null, p2_null=p2_null,
                              p1_alt=p1_alt, p2_alt=p2_alt, n1=n1, 
                              n2=n2, r1=r1, r2=r2)
            class(ret) <- c("two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
          })

# Create a two_stage_curtail_trial
# Case 2:  User inputs p, n
# Using default values of alpha and prob_early
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="numeric", n2="numeric", n_total="missing",
                    r1="missing", r2="missing", 
                    prob_early="missing", alpha="missing"),
          function(p1_null, p2_null, p1_alt, p2_alt, n1, n2) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null), 
                                       p_alt=c(p1_alt, p2_alt), 
                                       n=c(n1, n2), n_total=NULL, r=NULL, 
                                       prob_early=NULL, alpha=NULL)
            r <- two_stage_critical_values(c(n1, n2), 
                                           c(p1_null, p2_null))
            ret <- data.frame(p1_null=p1_null, p2_null=p2_null,
                              p1_alt=p1_alt,p2_alt=p2_alt, n1=n1, 
                              n2=n2, r1=r[1], r2=r[2])
            class(ret) <- c("two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
          })


# Create a two_stage_curtail_trial
# Case 3:  User inputs p, n, prob_early, alpha
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="numeric", n2="numeric", n_total="missing",
                    r1="missing", r2="missing", 
                    prob_early="numeric", alpha="numeric"),
          function(p1_null, p2_null, p1_alt, p2_alt, n1, n2,
                   prob_early, alpha) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null), 
                                       p_alt=c(p1_alt, p2_alt), 
                                       n=c(n1, n2), n_total=NULL, r=NULL,
                                       prob_early=prob_early, alpha=alpha)
            
            r <- two_stage_critical_values(c(n1, n2), 
                                           c(p1_null, p2_null), 
                                           prob_early, alpha)
            
            ret <- data.frame(p1_null=p1_null, p2_null=p2_null,
                              p1_alt=p1_alt,p2_alt=p2_alt, n1=n1, 
                              n2=n2, r1=r[1], r2=r[2])
            class(ret) <- c("two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
          })
#' Create a two_stage_curtail_trial
#' Case 3a:  User inputs p, n, prob_early
#' @examples
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30), prob_early=0.1, alpha=0.1)
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="numeric", n2="numeric", n_total="missing",
                    r1="missing", r2="missing", 
                    prob_early="numeric", alpha="missing"),
          function(p1_null, p2_null, p1_alt, p2_alt, n1, n2,
                   prob_early, alpha) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null), 
                                       p_alt=c(p1_alt, p2_alt), 
                                       n=c(n1, n2), n_total=NULL, r=NULL,
                                       prob_early=prob_early, alpha=NULL)
            
            r <- two_stage_critical_values(c(n1, n2), 
                                           c(p1_null, p2_null), 
                                           pearly=prob_early)
            
            ret <- data.frame(p1_null=p1_null, p2_null=p2_null,
                              p1_alt=p1_alt,p2_alt=p2_alt, n1=n1, 
                              n2=n2, r1=r[1], r2=r[2])
            class(ret) <- c("two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
          })

#' Create a two_stage_curtail_trial
#' Case 3b:  User inputs p, n, alpha
#' @examples
#' trial <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n=c(6, 30), prob_early=0.1, alpha=0.1)
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="numeric", n2="numeric", n_total="missing",
                    r1="missing", r2="missing", 
                    prob_early="missing", alpha="numeric"),
          function(p1_null, p2_null, p1_alt, p2_alt, n1, n2,
                   prob_early, alpha) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null), 
                                       p_alt=c(p1_alt, p2_alt), 
                                       n=c(n1, n2), n_total=NULL, r=NULL,
                                       prob_early=NULL, alpha=alpha)
            
            r <- two_stage_critical_values(c(n1, n2), 
                                           c(p1_null, p2_null), 
                                           alpha=alpha)
            
            ret <- data.frame(p1_null=p1_null, p2_null=p2_null,
                              p1_alt=p1_alt,p2_alt=p2_alt, n1=n1, 
                              n2=n2, r1=r[1], r2=r[2])
            class(ret) <- c("two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
          })

# Create a two_stage_curtail_trial
# Case 4:  User inputs p, n_total, prob_early, alpha
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="missing", n2="missing", n_total="numeric",
                    r1="missing", r2="missing", 
                    prob_early="numeric", alpha="numeric"),
          function(p1_null, p2_null, p1_alt, p2_alt, n_total, 
                   prob_early, alpha) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null),
                                       p_alt=c(p1_alt, p2_alt),
                                       n=NULL, n_total=n_total, r = NULL, 
                                       prob_early=prob_early, 
                                       alpha=alpha)
            
            n1 <- seq_len(n_total-1)
            n2 <- n_total - n1
            n <- cbind(n1, n2)
            r <- matrix(apply(n, 1, function(x){
              two_stage_critical_values(x, c(p1_null, p2_null),
                                        prob_early, alpha)
            }), nrow=length(n1), ncol=2, byrow=TRUE)

            ret <- data.frame(p1_null=rep(p1_null, length(n1)),
                              p2_null=rep(p2_null, length(n1)),
                              p1_alt=rep(p1_alt, length(n1)),
                              p2_alt=rep(p2_alt, length(n1)),
                              n1=n1, n2=n2, r1=r[,1], r2=r[,2])
            if(min(ret$r1)==0){
              ret <- subset(ret, r1>0)
              rownames(ret) <- NULL
            }
            class(ret) <- c("two_stage_curtail_trial_sawtooth", 
                            "two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
})

#' Create a two_stage_curtail_trial
#' Case 4a:  User inputs p, n_total, prob_early
#' @examples
#' trials <- two_stage_curtail_trial(p1_null = 0.8, p2_null=0.2, 
#' p1_alt = 0.8, p2_alt = 0.4, n_total=36, prob_early=0.1, alpha=0.1)
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="missing", n2="missing", n_total="numeric",
                    r1="missing", r2="missing", 
                    prob_early="numeric", alpha="missing"),
          function(p1_null, p2_null, p1_alt, p2_alt, n_total, 
                   prob_early, alpha) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null),
                                       p_alt=c(p1_alt, p2_alt),
                                       n=NULL, n_total=n_total, r = NULL, 
                                       prob_early=prob_early, 
                                       alpha=NULL)
            
            n1 <- seq_len(n_total-1)
            n2 <- n_total - n1
            n <- cbind(n1, n2)
            r <- matrix(apply(n, 1, function(x){
              two_stage_critical_values(x, c(p1_null, p2_null),
                                        pearly=prob_early)
            }), nrow=length(n1), ncol=2, byrow=TRUE)
            
            ret <- data.frame(p1_null=rep(p1_null, length(n1)),
                              p2_null=rep(p2_null, length(n1)),
                              p1_alt=rep(p1_alt, length(n1)),
                              p2_alt=rep(p2_alt, length(n1)),
                              n1=n1, n2=n2, r1=r[,1], r2=r[,2])
            if(min(ret$r1)==0){
              ret <- subset(ret, r1>0)
              rownames(ret) <- NULL
            }
            class(ret) <- c("two_stage_curtail_trial_sawtooth", 
                            "two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
          })
#' Create a two_stage_curtail_trial
#' Case 4b:  User inputs p, n_total, alpha
#' @examples
#' trials <- two_stage_curtail_trial(p1_null = 0.8, p2_null=0.2, 
#' p1_alt = 0.8, p2_alt = 0.4, n_total=36, prob_early=0.1, alpha=0.1)
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="missing", n2="missing", n_total="numeric",
                    r1="missing", r2="missing", 
                    prob_early="missing", alpha="numeric"),
          function(p1_null, p2_null, p1_alt, p2_alt, n_total, 
                   prob_early, alpha) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null),
                                       p_alt=c(p1_alt, p2_alt),
                                       n=NULL, n_total=n_total, r = NULL, 
                                       prob_early=NULL, alpha=alpha)
            
            n1 <- seq_len(n_total-1)
            n2 <- n_total - n1
            n <- cbind(n1, n2)
            r <- matrix(apply(n, 1, function(x){
              two_stage_critical_values(x, c(p1_null, p2_null),
                                        alpha=alpha)
            }), nrow=length(n1), ncol=2, byrow=TRUE)
            
            ret <- data.frame(p1_null=rep(p1_null, length(n1)),
                              p2_null=rep(p2_null, length(n1)),
                              p1_alt=rep(p1_alt, length(n1)),
                              p2_alt=rep(p2_alt, length(n1)),
                              n1=n1, n2=n2, r1=r[,1], r2=r[,2])
            if(min(ret$r1)==0){
              ret <- subset(ret, r1>0)
              rownames(ret) <- NULL
            }
            class(ret) <- c("two_stage_curtail_trial_sawtooth", 
                            "two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
          })
# Create a two_stage_curtail_trial
# Case 5:  User inputs p, n_total
# Using default values of prob_early and alpha, with n_total
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="numeric", p2_null="numeric", 
                    p1_alt="numeric",  p2_alt="numeric", 
                    n1="missing", n2="missing", n_total="numeric",
                    r1="missing", r2="missing", 
                    prob_early="missing", alpha="missing"),
          function(p1_null, p2_null, p1_alt, p2_alt, n_total, 
                   prob_early, alpha) {
            two_stage_valid_parameters(p_null=c(p1_null, p2_null), 
                                       p_alt=c(p1_alt, p2_alt), n=NULL,
                                       n_total=n_total, r=NULL, 
                                       prob_early=NULL, alpha=NULL)
            
            n1 <- seq_len(n_total-1)
            n2 <- n_total - n1
            n <- cbind(n1, n2)
            r <- matrix(apply(n, 1, function(x){
              two_stage_critical_values(x, c(p1_null, p2_null))
            }), nrow=length(n1), ncol=2, byrow=TRUE)
            
            ret <- data.frame(p1_null=rep(p1_null, length(n1)),
                              p2_null=rep(p2_null, length(n1)),
                              p1_alt=rep(p1_alt, length(n1)),
                              p2_alt=rep(p2_alt, length(n1)),
                              n1=n1, n2=n2, r1=r[,1], r2=r[,2])
            if(min(ret$r1)==0){
              ret <- subset(ret, r1>0)
              rownames(ret) <- NULL
            }
            class(ret) <- c("two_stage_curtail_trial_sawtooth", 
                            "two_stage_curtail_trial", class(ret))
            ret$power <- power(ret)
            ret$significance <- significance(ret)
            ret$stage1_mean_ss <- stage1_sample_size(ret)
            ret$mean_ss_null <- expected_sample_size(ret)
            ret$mean_ss_alt <- expected_sample_size_alt(ret)
            ret$PET <- PET(ret)
            ret$minimax_prob <- minimax_probability(ret)
            ret
          })

#' Create a two_stage_curtail_trial
#' Case 6:  Anything else
#' @examples
#' trials <- two_stage_curtail_trial(p_null=c(0.8, 0.2), p_alt=c(0.8, 0.4), 
#' n_total=36)
#' @export
setMethod("two_stage_curtail_trial",
          signature(p1_null="ANY", p2_null="ANY", 
                    p1_alt="ANY",  p2_alt="ANY", 
                    n1="ANY", n2="ANY", n_total="ANY",
                    r1="ANY", r2="ANY", 
                    prob_early="ANY", alpha="ANY"),
          function(p1_null, p2_null, p1_alt, p2_alt, n1, n2, n_total, 
                   r1, r2, prob_early, alpha) {
            stop("Design parameters not entered correctly. Refer to documentation")
          })

#' @import ggplot2
#' @importFrom tidyr gather
#' @export
#' 
plot.two_stage_curtail_trial_sawtooth <- function(x, ...) {
  y <- subset(x, select=c("n1", "mean_ss_null", "minimax_prob"))
  names(y) <- c("n1", "Optimal Criteria", "Minimax Criteria")
  m <- gather(y, variable, value, -n1)
  ggplot(m, aes(x = n1, y = value)) + geom_line() + 
    facet_grid(variable ~ ., scales = "free_y") + xlab(expression(n[1]))+ 
    ylab("")
}

#' @export
significance.two_stage_curtail_trial <- function(x) {
  apply(x, 1, function(x){
    p <- as.numeric(x[c("p1_null", "p2_null")])
    n <- as.numeric(x[c("n1", "n2")])
    r <- as.numeric(x[c("r1", "r2")])
    
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
          p <- as.numeric(x[c("p1_alt", "p2_alt")])
          n <- as.numeric(x[c("n1", "n2")])
          r <- as.numeric(x[c("r1", "r2")])
          
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
expected_sample_size.two_stage_curtail_trial <- function(x) {
  
  apply(x, 1,
        function(x) {
          p <- as.numeric(x[c("p1_null", "p2_null")])
          n <- as.numeric(x[c("n1", "n2")])
          r <- as.numeric(x[c("r1", "r2")])
          
          # Probability of stopping early
          pearly <- prob_early_stop(p, n, r)  
          
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
            denom <- sum(sapply(r1:n1, function(j) {
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
            
            denom <- sum(sapply((n1-r1+1):n1, function(j) {
                            factorial(j-1)/factorial(j-1-n1+r1) * p1^(j-y1)
                            }))
            
            if(is.na(num/denom)){ 
              return(0)
            } else {
              return(num/denom)
            }
          }
          
          ### expected value of Y1 | continue to Stage 2
          EY1Continue <- sum(sapply(r[1]:n[1], function(i) {
                                      i * prY1Continue(i, r[1], n[1], p[1])
                                    }))
          
          ### expected value of Y1|stop early
          EY1Stop <- sum(sapply((n[1]-r[1]+1):n[1], function(i) {
                                  i*prY1Stop(i, r[1], n[1], p[1])
                                }))
          
          ### expected value of y2|continue to Stage 2
          # i = possible values of y1
          # j = possible values of x12
          py1 <- sapply(r[1]:n[1], function(i) {
                          prY1Continue(i, r[1], n[1], p[1])
                        })
          
          px12 <- sapply(0:r[1], function(j) {
                           dbinom(j, r[1], p[2]/p[1])
                         })
          
          #for the snb distribution,
          #s <- r[2]-j
          #t <- n[2]+n[1]-i-r[2]+j+1
          ey2 <- function(y1) {
            sapply(0:r[1], function(j) {
                     esnb(p[2], s=ifelse((r[2]-j)>0, (r[2]-j), 0), 
                          t=ifelse((n[2]+n[1]-y1-r[2]+j+1)>0, 
                                   n[2]+n[1]-y1-r[2]+j+1, 0))
                   })
          }
          
          a <- sapply(r[1]:n[1], function(y1) sum(px12*ey2(y1)))
          EY2Continue <- sum(py1*a)
          
          return((EY1Stop*pearly) + ((EY1Continue+EY2Continue)*(1-pearly))) 
        })
}

#' @export
expected_sample_size_alt <- function(x) {
  UseMethod("expected_sample_size_alt")
}

expected_sample_size_alt.two_stage_curtail_trial <- function(x) {
  apply(x, 1,
        function(x) {
          p <- as.numeric(x[c("p1_alt", "p2_alt")])
          n <- as.numeric(x[c("n1", "n2")])
          r <- as.numeric(x[c("r1", "r2")])
          
          # Probability of stopping early
          pearly <- prob_early_stop(p, n, r)  
          
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
            denom <- sum(sapply(r1:n1, function(j) {
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
            
            denom <- sum(sapply((n1-r1+1):n1, function(j) {
              factorial(j-1)/factorial(j-1-n1+r1) * p1^(j-y1)
            }))
            
            if(is.na(num/denom)){ 
              return(0)
            } else {
              return(num/denom)
            }
          }
          
          ### expected value of Y1 | continue to Stage 2
          EY1Continue <- sum(sapply(r[1]:n[1], function(i) {
            i * prY1Continue(i, r[1], n[1], p[1])
          }))
          
          ### expected value of Y1|stop early
          EY1Stop <- sum(sapply((n[1]-r[1]+1):n[1], function(i) {
            i*prY1Stop(i, r[1], n[1], p[1])
          }))
          
          ### expected value of y2|continue to Stage 2
          # i = possible values of y1
          # j = possible values of x12
          py1 <- sapply(r[1]:n[1], function(i) {
            prY1Continue(i, r[1], n[1], p[1])
          })
          
          px12 <- sapply(0:r[1], function(j) {
            dbinom(j, r[1], p[2]/p[1])
          })
          
          #for the snb distribution,
          #s <- r[2]-j
          #t <- n[2]+n[1]-i-r[2]+j+1
          ey2 <- function(y1) {
            sapply(0:r[1], function(j) {
              esnb(p[2], s=ifelse((r[2]-j)>0, (r[2]-j), 0), 
                   t=ifelse((n[2]+n[1]-y1-r[2]+j+1)>0, 
                            n[2]+n[1]-y1-r[2]+j+1, 0))
            })
          }
          
          a <- sapply(r[1]:n[1], function(y1) sum(px12*ey2(y1)))
          EY2Continue <- sum(py1*a)
          
          return((EY1Stop*pearly) + ((EY1Continue+EY2Continue)*(1-pearly))) 
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
          
          n1 <- as.numeric(x["n1"])
          r1 <- as.numeric(x["r1"])
          p1 <- as.numeric(x["p1_null"])
          
          tt <- n1 - r1 + 1  # Number of failures to stop early
          
          e <- esnb(p1, r1, tt)  # Expected value
          
          v <- vsnb(p1, r1, tt)   # variance
          sd <- sqrt(v)
          e
          
        })
}

#' @export
PET <- function(x) {
  UseMethod("PET")
}

PET.two_stage_curtail_trial <- function(x){
  apply(x, 1,
        function(x) {
          r1 <- as.numeric(x["r1"])
          n1 <- as.numeric(x["n1"])
          p1 <- as.numeric(x["p1_null"])
          
          j <- min(r1, n1) : n1 # range of X1
          
          1 - sum(dbinom(j, n1, p1))
        }) 
}

#' @export
minimax_probability <- function(x) {
  UseMethod("minimax_probability")
}

minimax_probability.two_stage_curtail_trial <- function(x){
  
  apply(x, 1,
        function(x) {
          n <- as.numeric(c(x["n1"], x["n2"]))
          p <- as.numeric(c(x["p1_null"], x["p2_null"]))
          r <- as.numeric(c(x["r1"], x["r2"]))
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
  x[which.min(x$minimax_prob),]
  
}


#' @export
optimal_design <- function(x) {
  UseMethod("optimal_design")
}

#' @export
optimal_design.two_stage_curtail_trial_sawtooth <- function(x){
  x[which.min(x$mean_ss_null),]
}

summary.two_stage_curtail_trial_sawtooth <- function(x){
  opt_design <- round(optimal_design(x), 3)
  class(opt_design) <- c("two_stage_curtail_trial", "data.frame")
  min_design <- round(minimax_design(x), 3)
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
  round(x, 3)
}

##################################################################
two_stage_valid_parameters <- function(p_null, p_alt, n, n_total, r, 
                                       prob_early, alpha){
  
  if(!is.null(p_null) && (any(p_null > 1) || any(p_null < 0))){
    stop("p1_null and p2_null should be between 0 and 1") 
  } 
  if(!is.null(p_alt) && (any(p_alt > 1) || any(p_alt < 0))){
    stop("p1_alt and p2_alt should be between 0 and 1") 
  } 
  if(!is.null(n) && any(!is.wholenumber(n))){
    stop("n1 and n2 must be be integers") 
  } 
  if(!is.null(n) && any(n < 1)){
    stop("n1 and n2 should be at least 1") 
  } 
  if(!is.null(r) && any(r < 0)){
    stop("r1 and r2 should not be less than 0")
  } 
  if(!is.null(r) && !is.null(n) && r[1] > n[1]){
    stop("r1 should be less than or equal to n1")
  } 
  # Must have: 0 <= r2 <= n1+n2
  if(!is.null(r) && !is.null(n) && r[2] > sum(n)){
    stop("r2 should be less than or equal to n1+n2")
  } 
  if(!is.null(r) && any(!is.wholenumber(r))){
    stop("r1 and r2 should be integers")
  } 
  if(!is.null(n_total) && !is.wholenumber(n_total)){
    stop("n_total should be an integer") 
  } 
  if(!is.null(n_total) && n_total < 2){
    stop("n_total should be at least 2")
  } 
  if(!is.null(alpha) && (alpha > 1 | alpha < 0)){
    stop("alpha should be between 0 and 1")
  } 
  if(!is.null(prob_early) && (prob_early > 1 || prob_early < 0)){
    stop("prob_early should be between 0 and 1")
  }
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
#' two_stage_critical_values(n=c( 5, 31), p=c(.8,.2), pearly=.1, alpha=.1)
two_stage_critical_values <- function(n, p, pearly=.1, alpha=.1) {

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
