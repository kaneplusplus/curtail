S = function(k, p, s) {
  if (p < 0 || p > 1) stop("p must be between zero and one.")
  if (s < 0) stop("s must be non-negative")
  ret = choose(k-1, s-1) * p^s * (1-p)^(k-s)
  ret[k < s] = 0
  ret
}

#' Stack the distribution by responders and non-responders.
#'
#' Stacked distribution function for the stopped negative binomial distribution.
#' @param x quantile
#' @param p success probability
#' @param s number of successes
#' @param t number of failures
#' @export
dsnb_stacked = function(x, p, s, t) {
  ret = cbind(x, S(x, p, s), S(x, 1-p, t))
  colnames(ret) = c("x", "s", "t")
  ret
}

#R = function(k, p, t) {
#  choose(k-1, k-t) * p^(k-t) * (1-p)^t
#}

#N = function(k, p, s) {
#  choose(k-1, s-1) * p^s * (1-p)^(k-s)
#}

#Rc = function(k, s, t, shape1, shape2) {
#  suppressWarnings({ret = choose(k-1, k-t) * beta(shape1 + k - t, t + shape2) / 
#    beta(shape1, shape2)})
#  ret[!is.finite(ret)] = 0
#  ret
#}

#Nc = function(k, s, t, shape1, shape2) {
#  suppressWarnings({ret = choose(k-1, s-1) * beta(shape1 + s, k - s + shape2) / 
#    beta(shape1, shape2)})
#  ret[!is.finite(ret)] = 0
#  ret
#}

dsnb_private = function(x, p, s, t) {
  k=NULL
  a = foreach(k=1:(s+t-1), .combine=c) %do% N(k, p, s)
  b = foreach(k=1:(s+t-1), .combine=c) %do% R(k, p, t)
  d = a + b
  inds = which(x %in% 1:length(d))
  ret = rep(0, length(x))
  ret[inds] = d[x[inds]]
  ret
}

## Remember, shape1 and shape2 are data plus priors.
#dsnbc_private = function(x, s, t, shape1, shape2) {
#  k = NULL
#  a = foreach(k=1:(s+t-1), .combine=c) %do% Nc(k, s, t, shape1, shape2)
#  b = foreach(k=1:(s+t-1), .combine=c) %do% Rc(k, s, t, shape1, shape2)
#  d = a + b
#  inds = which(x %in% 1:length(d))
#  ret = rep(0, length(x))
#  ret[inds] = d[x[inds]]
#  ret
#}

#dsnb_private_stacked = dsnb_stacked

#dsnbc_private_stacked = function(x, shape1, shape2, s, t) {
#  k=i=NULL
#  a = foreach(k=1:(s+t-1), .combine=c) %do% Nc(k, s, t, shape1, shape2)
#  b = foreach(k=1:(s+t-1), .combine=c) %do% Rc(k, s, t, shape1, shape2)
#  d = a + b
#  u = a/sum(d)
#  r = b/sum(d)
#  ret = foreach (i=x, .combine=rbind) %do% {
#    v = c(i, 0, 0)
#    if (i %in% 1:length(d)) {
#      v[2] = a[i] #u[i] 
#      v[3] = b[i] #r[i] 
#    }
#    v
#  }
#  colnames(ret) = c("x", "s", "t")
#  rownames(ret) = NULL
#  ret
#}

#' The Stacked Plot
#'
#' The stacked plot of the probability mass function for the snb showing
#' the contributions from N (the top barrier) and R (the right barrier).
#' @param p the probability of a success on each coin flip. 
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @param offset an offset on the domain of the distribution. This is 
#' used when getting the conditional distribution where the domain does 
#' not start at 1.
#' @import ggplot2
#' @return a plot of the probability mass function.
#' @export
stacked_plot = function(x, s, t) {
  if (missing(s) && missing(t) && all(names(x) %in% c("x", "s", "t"))) {
    s = x$s
    t = x$t
    x = x$x
  }
  d = data.frame(list(x=x, s=s, t=t))
  d = melt(data=d, id.vars="x") 
  names(d)[names(d) == "variable"] = "Outcome"
  ggplot(data=d, aes(x=factor(x), y=value, fill=Outcome)) +
    geom_bar(position="stack", stat="identity") + xlab("k") +
    ylab("f(k,p,s,t)")
}

#' The Stopped Negative Binomial P.M.F. Plot
#'
#' The plot of the probability mass function for the SNB showing
#' the contributions from N (the top barrier) and R (the right barrier).
#' @param p the probability of a success on each coin flip. 
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @param offset an offset on the domain of the distribution. This is 
#' used when getting the conditional distribution where the domain does 
#' not start at 1.
#' @import ggplot2
#' @return a plot of the probability mass function.
#' @export
dsnb_plot = function(p, s, t, x, offset) {
  value = Outcome = k = NULL
  if (missing(x))
    x = min(s,t):(t+s-1)
  d = as.data.frame(
    dsnb_stacked(x, p=p, s=s, t=t))
  if (!missing(offset))
    d$x = d$x+offset
  d$y = apply(d[,2:3], 1, sum)
  ggplot(data=d, aes(x=factor(x), y=y)) + 
    geom_bar(position="stack", stat="identity") + xlab("k") +
    ylab("f(k,p,s,t)")
}

#' The Stopped Negative Binomial P.M.F. Stack-Plot
#'
#' The stacked plot of the probability mass function for the SNB showing
#' the contributions from N (the top barrier) and R (the right barrier) by
#' color.
#' @param p the probability of a success on each coin flip. 
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @param offset an offset on the domain of the distribution. This is 
#' used when getting the conditional distribution where the domain does 
#' not start at 1.
#' @import ggplot2
#' @return a plot of the probability mass function.
#' @export
dsnb_stack_plot = function(p, s, t, x, offset) {
  value = Outcome = k = NULL
  if (missing(x))
    x = min(s,t):(t+s-1)
  d = as.data.frame(
    dsnb_stacked(x, p=p, s=s, t=t))
  if (!missing(offset))
    d$x = d$x+offset
  d = melt(data=d, id.vars="x") 
  names(d)[names(d) == "variable"] = "Outcome"
  ggplot(data=d, aes(x=factor(x), y=value, fill=Outcome)) +
    geom_bar(position="stack", stat="identity") + xlab("k") +
    ylab("f(k,p,s,t)")
}

#' The Conditional Stopped Negative Binomial Density
#'
#' The conditional stacked snb density function. This function gets
#' the distribution of the stopping time when the binomial process has not 
#' reached one of its endpoints. The success probability is fitted using 
#' fit_flips with specified prior.
#' @param x quantile
#' @param shape the shape parameters of the beta prior.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @export
cdsnb_stacked = function(x, shape, s, t) {
  ret=foreach(k=x, .combine=rbind) %do% {
    rets=0
    rett=0
    normalizer = beta(shape[1], shape[2])
    if (s <= k && k <= s+t-1)
      rets = rets + choose(k-1, s-1) * beta(shape[1]+s, k-s+shape[2])/normalizer
    if (t <= k && k <= s+t-1)
      rett = rett + choose(k-1, t-1) * beta(shape[1]+k-t, t+shape[2])/normalizer
    c(k, rets, rett)
  }
  rownames(ret) = NULL
  colnames(ret) = c("k", "s", "t")
  as.data.frame(ret)
}

#' The Conditional Stopped Negative Binomial Density Plot
#' 
#' A plot of the stacked snb density function. The plot shows the distribution
#' of the stopping time when the binomial process has not reached one of 
#' its endpoints. The success probabilty is fitted using the fit_flips function
#' with specefied prior.
#' @param d a sequence of 1's and 0's corresponding to the binomial process.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param prior the shape parameters of the prior on the success probability.
#' @importFrom reshape2 melt
#' @export
cdsnb_stack_plot = function(d, shape, s, t) {
  value = Outcome = NULL
  x = cdsnb(d, shape, s, t)
  x = melt(data=x, id.vars="x") 
  names(x)[names(x) == "variable"] = "Outcome"
  qplot(x=factor(x), y=value, data=x, fill=Outcome, geom="bar", 
        position="stack", stat="identity", ylab="f(k,p,s,t)", xlab="k")
}

#' Stacked Plot of the Compound Stopped Negative Binomial Density 
#'
#' The stacked plot of the probability mass function for the snb showing
#' the contributions from N (the top barrier) and R (the right barrier).
#' @param d the data, a vector of 0 and 1 values.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param shape1 the value of the first shape parameter on the prior
#' @param shape2 the value of the second shape parameter on the prior
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @return a plot of the probability mass function.
#' @export
dsnbc_stack_plot = function(d, s, t, shape1=0.5, shape2=0.5,
                            x=min(s,t):(t+s-1)) {
  value = Outcome = NULL
  d = dsnbc_stack(d, s, t, shape1, shape2, x)
  d = melt(data=d, id.vars="x") 
  names(d)[names(d) == "variable"] = "Outcome"
  qplot(x=factor(x), y=value, data=d, fill=Outcome, geom="bar", 
        position="stack", stat="identity", ylab="f(k|x,s,t,alpha,beta)", xlab="k")
}

#' The "Stacked" Compound Negative Binomial Density Function
#' 
#' This function returns the "stacked" density function showing the 
#' contribution from each of the end points to the total mass.
#' @param d the data, a vector of 0 and 1 values.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param shape1 the value of the first shape parameter on the prior
#' @param shape2 the value of the second shape parameter on the prior
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @return the "stacked" density.
#' @export
dsnbc_stack = function(d, s, t, shape1=0.5, shape2=0.5,
                       x=min(s,t):(t+s-1)) {
  num_heads = sum(d)
  num_flips = length(d)
  as.data.frame(
    dsnbc_private_stacked(x, shape1+num_heads, shape2+num_flips-num_heads,s,t))
}


#' The Stopped Negative Binomial Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the stopped negative binomial distribution with parameters, 'prob',
#' 's', and 't'.
#' 
#' @rdname snb
#' @aliases psnb dsnb qsnb rsnb
#' @usage dsnb(x, prob, s, t)
#' psnb(q, prob, s, t)
#' qsnb(p, prob, s, t)
#' rsnb(n, prob, s, t)
#' @param x vector of quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param prob probility of success on each coin flip.
#' @param s the ceiling for the snb process
#' @param t the length of the the process can run for.
#' @return 'dsnb' give the density, 'psnb' give the distribution function
#' function, 'qsnb' gives the quantile function, 'rsnb' generates random
#' deviates.
#' @importFrom foreach foreach %do%
#' @export
dsnb = function(x, prob, s, t) {
  if (s < 1) stop("dsnb s-parameter must be at least 1")
  if (t < 1) stop("dsnb t-parameter must be at least 1")
  if (any(prob > 1) || any(prob < 0))
    stop("dsnb prob-parameter must be between zero and one inclusive")
  apply(dsnb_stacked(x, prob, s, t)[,2:3], 1, sum)
}

#' @export
rsnb = function(n, prob, s, t) {
  # Get the distribution function.
  support = min(s,t):(t+s-1)
  ps = dsnb(support, prob, s, t)
  sample(support, n, replace=TRUE, prob=ps)
}

##' Simulate the binomial process
##'
##' Generate coin flip trajectories that stop after either s heads or t tails.
##' @param n the number of trajectories to simulate.
##' @param prob the probability of a head.
##' @param s the number of heads to stop at.
##' @param t the number of tails to stop at.
##' @param drop if TRUE then return the results as a matrix. Otherwise return as a list.
## @export
#snb_flips = function(n, prob, s, t, drop=TRUE) {
#  if (length(prob) > 1)
#    stop("rsnb prob-parameter must have length 1")
#  flips = foreach(i=1:n) %do% {
#    flip = rbinom(s+t-1, 1, prob=prob)
#    path = c(cumsum(flip), sum(flip))
#    m = which(path >= s)
#    if (length(m) == 0) {
#      m = s+t-1
#    } else {
#      m = min(m)
#    }
#    r = which(path < 0:(s+t-1)-(t-s+1))
#    if (length(r) == 0) {
#      r = s+t-1
#    } else {
#      r = min(r)
#    }
#    flip[1:min(m, r)]
#  }
#  if (n == 1 && drop) {
#    flips = unlist(flips)
#  }
#  flips
#}

#flips_to_zplot_df = function(flips) {
#  d = data.frame(k=0:length(flips))
#  d$head = c(0, cumsum(flips))
#  d$tail= c(0, cumsum(!(flips)))
#  d$headEnd = c(d$head[-1], NA)
#  d$tailEnd = c(d$tail[-1], NA)
#  d
#}

#' The Z-Plot for the Binomial Process
#'
#' Visualize the stopped Bernoulli process with horizontal axis counting 
#' successes and vertical axis counting failure.
#'
#' @param flips the sequence of coing flips (1's and 0's) to visualize.
#' Note that this can be a list in which case multiple processes will be 
#' shown.
#' @param s the top barrier for the Bernoulli process.
#' @param t the right barrier for the Bernoulli process.
#' @param show_arrows should arrows be shown in the Bernoullis process path?
#' @param unif_jitter for multiple flip paths, how much jitter to add 
#' (default is 0.2).
#' @param xlab the name of the x axis.
#' @param ylab the name of the y axis.
#' @importFrom grid arrow
#' @examples
#' flips = c(0, 0, 1)
#' zplot(flips, 2, 3)
#' @export
zplot = function(flips, s, t, show_arrows=TRUE, unif_jitter=0.2, xlab=NULL,
                 ylab=NULL) {
  p = tailEnd = headEnd = num = NULL
  if (!is.list(flips)) {
    d =flips_to_zplot_df(flips)
    if (show_arrows) {
      p = ggplot(data=na.omit(d)) + 
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd), arrow=arrow()) 
    } else {
      p = ggplot(data=na.omit(d)) + 
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd)) 
    }
  } else {
    flip_set = lapply(flips, flips_to_zplot_df)
    for (i in 1:length(flip_set)) {
      flip_set[[i]] = na.omit(flip_set[[i]])
      flip_set[[i]]$num = as.factor(i)
      if (tail(flip_set[[i]]$headEnd, 1) == s) {
        # We hit the top barrier. Jitter on the x values
        flip_set[[i]]$tail= flip_set[[i]]$tail + runif(nrow(flip_set[[i]]),
                                                       -unif_jitter, unif_jitter) 
        flip_set[[i]]$tailEnd = flip_set[[i]]$tailEnd + 
          runif(nrow(flip_set[[i]]), -unif_jitter, unif_jitter)
      } else {
        flip_set[[i]]$head = flip_set[[i]]$head + runif(nrow(flip_set[[i]]), 
                                                        -unif_jitter, unif_jitter)
        flip_set[[i]]$headEnd = flip_set[[i]]$headEnd + 
          runif(nrow(flip_set[[i]]), -unif_jitter, unif_jitter)
      }
      # Make sure that the paths "connect".
      for (j in nrow(flip_set[[i]]):2) {
        flip_set[[i]][j, c("head", "tail")] = 
          flip_set[[i]][j-1, c("headEnd", "tailEnd")]
      }
    }
    d = Reduce(rbind, flip_set)
    if (show_arrows) {
      p = ggplot(data=na.omit(d)) + 
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd, group=num), arrow=arrow())
    } else {
      p = ggplot(data=na.omit(d)) + 
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd, group=num))
    }
  }
  p = p+scale_x_continuous(breaks=0:t, limits=c(-unif_jitter, t)) +
    scale_y_continuous(breaks=0:s, limits=c(-unif_jitter, s)) +
    geom_segment(x=0, y=s, xend=t-1, yend=s, color="red") +
    geom_segment(x=t, y=0, xend=t, yend=s-1, color="green")
  if (!is.null(xlab))
    p = p + xlab(xlab)
  if (!is.null(ylab))
    p = p + ylab(ylab)
  p
}

stairs = function(p, xstart, xend) {
  x = c(xstart, rep((xstart+1):xend, each=2))
  y = rep(0:(xend-xstart), each=2)
  y = y[-length(y)]
  for (i in 1:(length(x)-1)) {
    p = p + geom_segment(x=x[i], y=y[i], xend=x[i+1], yend=y[i+1],
                         color="green")
  }
  p
}

flips_to_kplot_df = function(flips) {
  d = data.frame(k=0:length(flips))
  d$head = c(0, cumsum(flips))
  d$tail= c(0, cumsum(1-flips))
  d$headEnd = c(d$head[-1], NA)
  d$tailEnd = c(d$tail[-1], NA)
  d$path = c(0, cumsum(flips))
  d$k = 0:(nrow(d)-1)
  d
}

#' The K-Plot for the Binomial Process
#'
#' Visualize the stopped Bernoulli process with a horizontal step axis and a 
#' vertical axis counting the number of successes.
#'
#' @param flips the sequence of coin flips (1's and 0's) to visualize.
#' @param s the top barrier for the Bernoulli process.
#' @param t the right barrier for the Bernoulli process.
#' @param bw should the plot be in black and white?
#' @examples
#' flips = c(0, 0, 1)
#' kplot(flips, 2, 3)
#' @export
kplot = function(flips, s, t, bw=FALSE) {
  if (!is.list(flips)) {
    d = flips_to_kplot_df(flips)
    if (bw) {
      p = qplot(k, path, data = d, geom = "line") +
        scale_x_continuous(breaks = 0:(t + s), limits = c(0, t + s)) +
        scale_y_continuous(breaks = 0:s, limits=c(0, s+0.15)) +
        geom_segment(x=s, y=s, xend=(t+s-1), yend=s, linetype=2) +
        geom_segment(x=t, y=0, xend=(s+t-1), yend=s-1, linetype=2)
    } else {
      p = qplot(k, path, data = d, geom = "line") +
        scale_x_continuous(breaks = 0:(t + s), limits = c(0, t + s)) +
        scale_y_continuous(breaks = 0:s, limits=c(0, s+0.15)) +
        geom_segment(x=s,y=s,xend=(t+s-1),yend=s, color="green", linetype=1) +
        geom_segment(x=t, y=0, xend=(s+t-1), yend=s-1, col="red")
    }
  } else {
    flip_set = lapply(flips, flips_to_kplot_df)
    for (i in 1:length(flip_set)) {
      flip_set[[i]]$num = as.factor(i)
      flip_set[[i]]$k = jitter(flip_set[[i]]$k)
      flip_set[[i]]$k[flip_set[[i]]$k < 0] = 0
    }
    d = Reduce(rbind, flip_set)[, -(4:5)]
    if (bw) {
      p = qplot(k, path, data = d, geom = "path", group = num) +
        scale_x_continuous(breaks=0:(t+s), limits = c(0, t+s)) +
        geom_segment(x = s, y = s, xend = (t + s - 1), yend = s,
                     linetype=2) +
        geom_segment(x=t, y=0, xend=(s+t-1), yend=s-1, linetype=2)
    } else {
      p = qplot(k, path, data = d, geom = "path", group = num) +
        scale_x_continuous(breaks=0:(t+s), limits = c(0, t+s)) +
        geom_segment(x = s, y = s, xend = (t + s - 1), yend = s,
                     color = "green") +
        geom_segment(x=t, y=0, xend=(s+t-1), yend=s-1, col="red")
    }
  }
  p
}


#' @export
psnb = function(q, prob, s, t) {
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support = min(s, t):(t+s-1)
  cdf = c(rep(0, support[1]-1), cumsum(dsnb(support, prob, s, t)))
  qs = floor(q)
  qs[qs < support[1]] = support[1]-1
  qs[qs > support[length(support)]] = support[length(support)]
  cdf[qs]
}

#' @export
qsnb = function(p, prob, s, t) {
  pr = NULL
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support = min(s, t):t
  cdf = c(rep(0, support[1]-1), cumsum(dsnb(support, prob, s, t)))
  ret = foreach(pr=p, .combine=c) %do% {
    r = NA
    if (!is.na(pr)) {
      r = which(pr < cdf)[1]
      if (is.na(r))
        r = support[length(support)]
    }
    if (pr > 1 || pr < 0)
      r = NaN
    r
  }
  ret[ret < support[1]-1] = support[1] - 1
  ret
}

#' Expected Value of the SNB Distribution
#' 
#' Find the expected size of an SNB distribution with specified parameters.
#' @param p success probability
#' @param s number of successes 
#' @param t number of failures
#' @export
esnb = function(p, s, t) {
  ds = dsnb_stacked(min(s,t):(s+t-1), p, s, t)
  ds[,2:3] = ds[,1] * ds[,2:3]
  sum(as.vector(ds[,2:3]))
}



#' Variance of the SNB Distribution
#' 
#' Find the variance of the SNB distribution with specified parameters.
#' @param p success probability
#' @param s number of successes 
#' @param t number of failures
#' @export
vsnb = function(p, s, t) {
  ds = dsnb_stacked(min(s,t):(s+t-1), p, s, t)
  ds[,2:3] = ds[,1]^2 * ds[,2:3]
  sum(as.vector(ds[,2:3])) - esnb(p, s, t)^2
}

#' Expected Value of the Conditional SNB Distribution
#' 
#' Find the expected size of the conditional SNB distribution with specified 
#' parameters.
#' @param shape the shape parameters of the beta prior.
#' @param s number of successes 
#' @param t number of failures
#' @export
ecsnb = function(shape, s, t) {
  ds = cdsnb_stacked(min(s,t):(s+t-1), shape, s, t)
  ds[,2:3] = ds[,1] * ds[,2:3]
  sum(as.vector(ds[,2:3]))
}

#' Variance of the Conditional SNB Distribution
#' 
#' Find the variance of the conditional SNB distribution with specified 
#' parameters.
#' @param shape the shape parameters of the beta prior.
#' @param s number of successes.
#' @param t number of failures.
#' @export
vcsnb = function(shape, s, t) {
  ds = cdsnb_stacked(min(s,t):(s+t-1), shape, s, t)
  ds[,2:3] = ds[,1]^2 * ds[,2:3]
  sum(as.vector(ds[,2:3])) - ecsnb(shape, s, t)^2 
}

#' Expected size for the DKZ 2-stage trial
#' 
#' Find the expected size of the DKZ trial with specified parameters.
#' @param n1 maximum number of enrollees in the first stage.
#' @param r1 number of successes to move to stage-2.
#' @param p1 success probability in stage-1.
#' @param n2 maximum number of enrollees in stage-2.
#' @param r2 number of successes in stage-2 for success endpoint.
#' @param p2 success probability in stage-2.
#' @export
edkz = function(n1, r1, p1, n2, r2, p2) {
  EY1 = esnb(p1, r1, n1-r1+1)
  X12 = cbind(0:r1, dbinom(0:r1, r1, p2/p1))
  stage1_success = cbind((n1-r1):n1, S((n1-r1):n1, p1, r1))
  EY2 = 0
  for (i in 1:nrow(stage1_success)) {
    for (j in 1:nrow(X12)) {
      EY2 = EY2 + stage1_success[i,2] * X12[j,2] * 
        esnb(p2, r2-X12[j,1], n2+n1-stage1_success[i,1]-r2-r1+X12[j,1]+1)
    }
  }
  #  for (i in 1:nrow(X12)) {
  #    EY2 = EY2 + X12[i,2] * esnb(p2, r2-X12[i,1], n2-r2-r1+X12[i,1]+1)
  #  }  
  EY1 + EY2
}


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
#' successes to reject the null hypothesis (r2) or a scalar containing r1
#' @examples
#' ph2early(p = .8, n = 25, r = 18)
#' ph2early(p = c(.8, .2), n = c(25, 9), r = 18)
#' @export
ph2early = function( p, n, r) {
  if( any(p>1) || any(p<0)) stop('p must be in [0, 1]')

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
#' ph2Eearly(p = c(.8, .2), n = c(5, 31), r = 3)
#' @export
ph2Eearly = function(p, n, r) {
  if(length(p) > 1){ p1 <- p[1] # p can be vector or scalar
  
  }else p1 <- p
  
  if(length(n) > 1){ n1 <- n[1] # n can be vector or scalar
  
  }else n1 <- n
  
  if(length(r) > 1){ r1 <- r[1] # r can be vector or scalar
  
  }else r1 <- r
  
  
  
  tt <- n1 - r1 + 1  # Number of failures to stop early
  
  e <- sum((1 : n1) * ph2snb(p1, r1, tt))  # Expected value
  
  v <- sum((1 : n1) ^ 2 * ph2snb(p1, r1, tt)) - e ^ 2   # variance
  
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


#' Find the minimax design
#'
#'Calculates the minimax probability for each combination of n1 and n2 for a given total n
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param ntot scalar containing the total sample size planned for the trial (n1+n2) 
#' @param pearly desired probability of early stopping (default = .1).
#' @param alpha desired significance level (default = .1).
#' @examples 
#' minimax(c(.8, .2), 36)
#' minimax(c(.7, .3), 40, pearly = .08, alpha=.1)
minimax = function(p, ntot, pearly = .1, alpha = .1){
  
  mmax <- NULL
  nOne <- NULL
  nTwo <- NULL
  # for each combination of n1 and n2, calculate probability of using maximum sample size
  for (i in 1:(ntot-1)){
    nn <- c(i, ntot - i)
    r <- ph2crit(n=nn, p=p, pearly = pearly, alpha = alpha)
    if(r[1]>0){
      prob <- ph2mmax(p=p, n = nn, r)
      mmax <- c(mmax, prob)
      nOne <- c(nOne, nn[1])
      nTwo <- c(nTwo, nn[2])
    }
    
  }
  mmax <- data.frame(nOne, nTwo, mmax)
  names(mmax) <- c("n1", "n2", "Probability of Maximum Sample Size")
  
  return(mmax[order(mmax[,3]),])
}

#' Find the optimal design
#'
#'Calculates the expected sample size under curtailed sampling for each combination of 
#'n1 and n2 for a given total n
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param ntot scalar containing the total sample size planned for the trial (n1+n2) 
#' @param pearly desired probability of early stopping (default = .1).
#' @param alpha desired significance level (default = .1).
#' @examples 
#' optimal(c(.8, .2), 36)
#' optimal(c(.7, .3), 40, pearly = .08, alpha=.1)
#' @export
optimal = function(p, ntot, pearly = .1, alpha = .1){
  
  ess <- NULL
  nOne <- NULL
  nTwo <- NULL
  for (i in 1:(ntot))
  {
    nn <- c(i, ntot - i)
    r <- ph2crit(n=nn, p=p, pearly = pearly, alpha = alpha)
    if(r[1]>0){
      ss <- ph2Ess(p=p, n=nn, r)
      ess <- c(ess, ss)
      nOne <-c(nOne, nn[1])
      nTwo <-c(nTwo, nn[2])
    }
  }
  ess <- data.frame(nOne, nTwo, ess)
  names(ess) <- c("n1", "n2", "Expected Sample Size")
  return(ess[order(ess[,3]),])
}

#' Finds the minimax and optimals designs
#'
#' @param p vector containing the probability of successful outcomes
#'          in Stage 1 (p1) and Stage 2 (p2) 
#' @param ntot scalar containing the total sample size planned for the trial (n1+n2) 
#' @param pearly desired probability of early stopping (default = .1).
#' @param alpha desired significance level (default = .1).
#' @examples 
#' ph2designs(c(.8, .2), 36)
#' ph2designs(c(.7, .3), 40, pearly = .08, alpha=.1)
#' @export
ph2designs = function(p, ntot, pearly = .1, alpha = .1){
  opt <- optimal(p, ntot, pearly, alpha)
  mini <- minimax(p, ntot, pearly, alpha)
  optimalDes <- opt[which.min(opt[,3]),]
  minimaxDes <- mini[which.min(mini[,3]),]
  nOpt <- as.vector(c(optimalDes[1,1], optimalDes[1,2]))
  nMini <- as.vector(c(minimaxDes[1,1], minimaxDes[1,2]))
  rOpt <- ph2crit(c(optimalDes[1,1], optimalDes[1,2]), p, pearly, alpha)
  rMini <- ph2crit(c(minimaxDes[1,1], minimaxDes[1,2]), p, pearly, alpha)
  
  optimalDesign <- cbind(p[1], nOpt[1], rOpt[1], p[2], nOpt[2], rOpt[2], alpha, ph2early(p, nOpt, rOpt), 
                         optimalDes[1,3], ph2mmax(p, nOpt, rOpt))
  minimaxDesign <- cbind(p[1], nMini[1], rMini[1], p[2], nMini[2], rMini[2], alpha, ph2early(p, nMini, rMini), 
                         ph2Ess(p, nMini, rMini), minimaxDes[1,3])
  designs <- data.frame(rbind(optimalDesign, minimaxDesign))
  rownames(designs) <- c("Optimal", "Minimax")
  colnames(designs) <- c("p1", "n1", "r1", "p2", "n2", "r2", "Alpha", "PET", "ECSS", "P(MaxSS)")
  class(designs) = c("ph2_design", "data.frame")
  return(designs)
}

#' @importFrom ggplot2 ggplot
#' @export
plot.ph2_design = function(x, ...) {
  # Get the arguments from dots
  dot_args = list(...)
  
  # Fill in the value for pearly.
  if ("pearly" %in% names(dot_args)) {
    pearly = dot_args$pearly
  } else {
    pearly = .1
  }
  
  # Fill in the value for alpha.
  if ("alpha" %in% names(dot_args)) {
    alpha = dot_args$alpha
  } else {
    alpha = .1
  }
  
  p <- c(x[1, "p1"], x[1, "p2"])
  ntot <- x[1, "n1"]+x[1, "n2"]
  
  opt <- optimal(p, ntot, pearly, alpha)
  opt <- opt[order(opt[,1]),]
  mini <- minimax(p, ntot, pearly, alpha)
  mini <- mini[order(mini[,1]),]
  
  
  df <- data.frame(n1 = c(opt[,1], mini[,1]), values = c(opt[,3], mini[,3]), type=c(rep("Optimal Criteria", dim(opt)[1]), rep("Minimax Criteria", dim(mini)[1])))
  ggplot(data=df, aes(x=n1, y=values)) + geom_line() + facet_grid(type ~ ., scales="free")+ 
    xlab(expression(n[1]))+ylab("")
  
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
      for (i in 0 : r[1])     # Y2 is convolution sum of two binomials
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
#' ph2snb(p = .8, s = 3, t = 5)
ph2snb = function(p, s, t) {
  
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
#' @param x the number to test.
#' @param tot the tolerance for x being a whole number (default \code{.Machine$double.eps ^ 0.5})
#' @examples
#' is.wholenumber(4)
#' is.wholenumber(3.5)
is.wholenumber = function(x, tol = .Machine$double.eps ^ 0.5){
    abs(x - round(x)) < tol
}
