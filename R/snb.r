
# The truncated negative binomial function.
S <- function(k, p, s) {
  if (p < 0 || p > 1) stop("p must be between zero and one.")
  if (s < 0) stop("s must be non-negative")
  ret <- choose(k-1, s-1) * p^s * (1-p)^(k-s)
  ret[k < s] <- 0
  ret
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
dsnb <- function(x, prob, s, t) {
  if (s < 1) stop("dsnb s-parameter must be at least 1")
  if (t < 1) stop("dsnb t-parameter must be at least 1")
  if (any(prob > 1) || any(prob < 0))
    stop("dsnb prob-parameter must be between zero and one inclusive")
  apply(dsnb_stacked(x, prob, s, t)[, 2:3, drop=FALSE], 1, sum)
}

#' @export
rsnb <- function(n, prob, s, t) {
  # Get the distribution function.
  support <- min(s,t):(t+s-1)
  ps <- dsnb(support, prob, s, t)
  sample(support, n, replace=TRUE, prob=ps)
}

#' @export
psnb <- function(q, prob, s, t) {
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support <- min(s, t):(t+s-1)
  cdf <- c(rep(0, support[1]-1), cumsum(dsnb(support, prob, s, t)))
  qs <- floor(q)
  qs[qs < support[1]] <- support[1]-1
  qs[qs > support[length(support)]] <- support[length(support)]
  cdf[qs]
}

#' @export
qsnb <- function(p, prob, s, t) {
  pr <- NULL
  if (length(prob) > 1)
    stop("psnb prob-parameter may only have length 1")
  support <- min(s, t):t
  cdf <- c(rep(0, support[1]-1), cumsum(dsnb(support, prob, s, t)))
  ret <- foreach(pr=p, .combine=c) %do% {
    r <- NA
    if (!is.na(pr)) {
      r <- which(pr < cdf)[1]
      if (is.na(r))
        r <- support[length(support)]
    }
    if (pr > 1 || pr < 0)
      r <- NaN
    r
  }
  ret[ret < support[1]-1] <- support[1] - 1
  ret
}

#' Expected Value of the SNB Distribution
#' 
#' Find the expected size of an SNB distribution with specified parameters.
#' @param p success probability
#' @param s number of successes 
#' @param t number of failures
#' @examples 
#' esnb(p=0.2, s=7, t=11)
#' @export
esnb <- function(p, s, t) {
  ds <- dsnb_stacked(min(s,t):(s+t-1), p, s, t)
  ds[,2:3] <- ds[,1] * ds[,2:3]
  sum(as.vector(ds[,2:3]))
}

#' Variance of the SNB Distribution
#' 
#' Find the variance of the SNB distribution with specified parameters.
#' @param p success probability
#' @param s number of successes 
#' @param t number of failures
#' @examples
#' vsnb(p=0.2, s=7, t=11)
#' @export
vsnb <- function(p, s, t) {
  ds <- dsnb_stacked(min(s,t):(s+t-1), p, s, t)
  ds[,2:3] <- ds[,1]^2 * ds[,2:3]
  sum(as.vector(ds[,2:3])) - esnb(p, s, t)^2
}

#' Expected Value of the Conditional SNB Distribution
#' 
#' Find the expected size of the conditional SNB distribution with specified 
#' parameters.
#' @param shape the shape parameters of the beta prior.
#' @param s number of successes 
#' @param t number of failures
#' @examples 
#' ecsnb(c(0.5, 0.5), s=7, t=11)
#' @export
ecsnb <- function(shape, s, t) {
  ds <- cdsnb_stacked(min(s,t):(s+t-1), shape, s, t)
  ds[,2:3] <- ds[,1] * ds[,2:3]
  sum(as.vector(ds[,2:3]))
}

#' Variance of the Conditional SNB Distribution
#' 
#' Find the variance of the conditional SNB distribution with specified 
#' parameters.
#' @param shape the shape parameters of the beta prior.
#' @param s number of successes.
#' @param t number of failures.
#' @importFrom foreach foreach %do%
#' @examples 
#' vcsnb(c(0.5, 0.5), s=7, t=11)
#' @export
vcsnb <- function(shape, s, t) {
  ds <- cdsnb_stacked(min(s,t):(s+t-1), shape, s, t)
  ds[,2:3] <- ds[,1]^2 * ds[,2:3]
  sum(as.vector(ds[,2:3])) - ecsnb(shape, s, t)^2 
}

#' Stack the Stopped Negative Binomial distribution by responders and 
#' non-responders.
#'
#' Stacked distribution function for the Stopped Negative Binomial distribution.
#' @param x quantile
#' @param p success probability
#' @param s number of successes
#' @param t number of failures
#' @examples
#' dsnb_stacked(14, 0.2, 7, 11) 
#' @export
dsnb_stacked <- function(x, p, s, t) {
  ret <- cbind(x, S(x, p, s), S(x, 1-p, t))
  colnames(ret) <- c("x", "s", "t")
  ret
}

#' The Stacked Plot of the Stopped Negative Binomial Distribution
#'
#' The stacked plot of the probability mass function for the snb showing
#' the contributions from N (the top barrier) and R (the right barrier).
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @import ggplot2
#' @importFrom tidyr gather
#' @return a plot of the probability mass function.
#' @examples 
#' stacked_plot(x=7:17, s=7, t=11)
#' @export
stacked_plot <- function(x, s, t) {
  Outcome <- value <- NULL
  if (missing(s) && missing(t) && all(names(x) %in% c("x", "s", "t"))) {
    s <- x$s
    t <- x$t
    x <- x$x
  }
  d <- data.frame(list(x=x, s=s, t=t))
  d <- gather(d, Outcome, value, -x)
  names(d)[names(d) == "variable"] <- "Outcome"
  ggplot(data=d, aes(x=factor(x), y=value, fill=Outcome)) +
    geom_bar(position="stack", stat="identity") + xlab("k") +
    ylab("f(k,p,s,t)")
}

#' The Stopped Negative Binomial P.M.F. Plot
#'
#' The plot of the probability mass function for the Stopped Negative Binomial.
#' @param p the probability of a success on each coin flip. 
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @param offset an offset on the domain of the distribution. This is 
#' used when getting the conditional distribution where the domain does 
#' not start at 1.
#' @import ggplot2
#' @return a plot of the probability mass function.
#' @examples 
#' dsnb_plot(p=0.2, s=7, t=11)
#' @export
dsnb_plot <- function(p, s, t, x, offset) {
  y <- NULL
  value <- Outcome <- k <- NULL
  if (missing(x))
    x <- min(s,t):(t+s-1)
  d <- as.data.frame(
    dsnb_stacked(x, p=p, s=s, t=t))
  if (!missing(offset))
    d$x <- d$x+offset
  d$y <- apply(d[,2:3], 1, sum)
  ggplot(data=d, aes(x=factor(x), y=y)) +
    geom_bar(position="stack", stat="identity") + xlab("k") +
    ylab("f(k,p,s,t)")
}

#' The Stopped Negative Binomial P.M.F. Stack-Plot
#'
#' The stacked plot of the probability mass function for the Stopped Negative 
#' Binomial showing the contributions from N (the top barrier) and R (the 
#' right barrier) by color.
#' @param p the probability of a success on each coin flip. 
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @param x the range of the distribution (defaults to min(s,t):(t+s-1)).
#' @param offset an offset on the domain of the distribution. This is 
#' used when getting the conditional distribution where the domain does 
#' not start at 1.
#' @import ggplot2
#' @importFrom tidyr gather
#' @return a plot of the probability mass function.
#' @examples 
#' dsnb_stack_plot(p=0.2, s=7, t=11)
#' @export
dsnb_stack_plot <- function(p, s, t, x, offset) {
  Outcome <- value <- NULL
  if (missing(x))
    x <- min(s,t):(t+s-1)
  d <- as.data.frame(
    dsnb_stacked(x, p=p, s=s, t=t))
  if (!missing(offset))
    d$x <- d$x+offset
  d <- gather(d, Outcome, value, -x)
  ggplot(data=d, aes(x=factor(x), y=value, fill=Outcome)) +
    geom_bar(position="stack", stat="identity") + xlab("k") +
    ylab("f(k,p,s,t)")
}

#' The Conditional Stopped Negative Binomial Density
#'
#' The conditional stacked snb density function. This function gets
#' the distribution of the stopping time when the binomial process has not 
#' reached one of its endpoints. 
#' @param x quantile
#' @param shape the shape parameters of the beta prior.
#' @param s the top barrier for the snb process.
#' @param t the right barrier for the snb process.
#' @importFrom foreach foreach %do%
#' @examples 
#' cdsnb_stacked(8:11, c(0.5, 0.5), s=7, t=11)
#' @export
cdsnb_stacked <- function(x, shape, s, t) {
  k <- NULL
  ret <- foreach(k=x, .combine=rbind) %do% {
    rets <- 0
    rett <- 0
    normalizer <- beta(shape[1], shape[2])
    if (s <= k && k <= s+t-1) {
      rets <- rets + choose(k-1, s-1) *
        beta(shape[1]+s, k-s+shape[2])/normalizer
    }
    if (t <= k && k <= s+t-1) {
      rett <- rett + choose(k-1, t-1) *
        beta(shape[1]+k-t, t+shape[2])/normalizer
    }
    c(k, rets, rett)
  }
  rownames(ret) <- NULL
  colnames(ret) <- c("k", "s", "t")
  as.data.frame(ret)
}

#' Format a data.frame for plotting by the zplot function
#' @param flips the seqence of ones and zeros denoting responses and 
#' non-responses respecitively.
flips_to_zplot_df <- function(flips) {
  d <- data.frame(k=0:length(flips))
  d$head <- c(0, cumsum(flips))
  d$tail <- c(0, cumsum(!(flips)))
  d$headEnd <- c(d$head[-1], NA)
  d$tailEnd <- c(d$tail[-1], NA)
  d
}

#' The Z-Plot for the Stopped Negative Binomial Process
#'
#' Visualize the Stopped Negative Binomial process with horizontal axis 
#' counting successes and vertical axis counting failure.
#'
#' @param flips the sequence of coin flips (1's and 0's) to visualize.
#' Note that this can be a list in which case multiple processes will be 
#' shown.
#' @param s the top barrier for the Bernoulli process.
#' @param t the right barrier for the Bernoulli process.
#' @param show_arrows should arrows be shown in the Bernoullis process path?
#' @param unif_jitter for multiple flip paths, how much jitter to add 
#' (default is 0.2).
#' @param xlab the name of the x axis.
#' @param ylab the name of the y axis.
#' @import ggplot2
#' @importFrom grid arrow
#' @importFrom stats na.omit runif
#' @importFrom utils tail head
#' @examples
#' flips <- c(0, 0, 1)
#' zplot(flips, 2, 3)
#' @export
zplot <- function(flips, s, t, show_arrows=TRUE, unif_jitter=0.2, xlab=NULL,
                  ylab=NULL) {
  p <- tailEnd <- headEnd <- num <- NULL
  if (!is.list(flips)) {
    d <- flips_to_zplot_df(flips)
    if (show_arrows) {
      p <- ggplot(data=na.omit(d)) +
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd), arrow=arrow())
    } else {
      p <- ggplot(data=na.omit(d)) +
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd))
    }
  } else {
    flip_set <- lapply(flips, flips_to_zplot_df)
    for (i in seq_along(flip_set)) {
      flip_set[[i]] <- na.omit(flip_set[[i]])
      flip_set[[i]]$num <- as.factor(i)
      if (tail(flip_set[[i]]$headEnd, 1) == s) {
        # We hit the top barrier. Jitter on the x values
        flip_set[[i]]$tail <- flip_set[[i]]$tail +
          runif(nrow(flip_set[[i]]), -unif_jitter, unif_jitter)
        flip_set[[i]]$tailEnd <- flip_set[[i]]$tailEnd +
          runif(nrow(flip_set[[i]]), -unif_jitter, unif_jitter)
      } else {
        flip_set[[i]]$head <- flip_set[[i]]$head +
          runif(nrow(flip_set[[i]]), -unif_jitter, unif_jitter)
        flip_set[[i]]$headEnd <- flip_set[[i]]$headEnd +
          runif(nrow(flip_set[[i]]), -unif_jitter, unif_jitter)
      }
      # Make sure that the paths "connect".
      for (j in nrow(flip_set[[i]]):2) {
        flip_set[[i]][j, c("head", "tail")] <-
          flip_set[[i]][j-1, c("headEnd", "tailEnd")]
      }
    }
    d <- Reduce(rbind, flip_set)
    if (show_arrows) {
      p <- ggplot(data=na.omit(d)) +
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd, group=num), arrow=arrow())
    } else {
      p <- ggplot(data=na.omit(d)) +
        geom_segment(mapping=aes(x=tail, y=head, xend=tailEnd,
                                 yend=headEnd, group=num))
    }
  }
  p <- p+scale_x_continuous(breaks=0:t, limits=c(-unif_jitter, t)) +
    scale_y_continuous(breaks=0:s, limits=c(-unif_jitter, s)) +
    geom_segment(x=0, y=s, xend=t-1, yend=s, color="red") +
    geom_segment(x=t, y=0, xend=t, yend=s-1, color="green")
  if (!is.null(xlab))
    p <- p + xlab(xlab)
  if (!is.null(ylab))
    p <- p + ylab(ylab)
  p
}

stairs <- function(p, xstart, xend) {
  x <- c(xstart, rep((xstart+1):xend, each=2))
  y <- rep(0:(xend-xstart), each=2)
  y <- y[-length(y)]
  for (i in 1:(length(x)-1)) {
    p <- p + geom_segment(x=x[i], y=y[i], xend=x[i+1], yend=y[i+1],
                         color="green")
  }
  p
}

flips_to_kplot_df <- function(flips) {
  d <- data.frame(k=0:length(flips))
  d$head <- c(0, cumsum(flips))
  d$tail <- c(0, cumsum(1-flips))
  d$headEnd <- c(d$head[-1], NA)
  d$tailEnd <- c(d$tail[-1], NA)
  d$path <- c(0, cumsum(flips))
  d$k <- 0:(nrow(d)-1)
  d
}


#' The K-Plot for the Stopped Negative Binomial Process
#'
#' Visualize the Stopped Negative Binomial process with a horizontal step axis and a 
#' vertical axis counting the number of successes.
#'
#' @param flips the sequence of coin flips (1's and 0's) to visualize.
#' @param s the top barrier for the Bernoulli process.
#' @param t the right barrier for the Bernoulli process.
#' @param bw should the plot be in black and white?
#' @examples
#' flips <- c(0, 0, 1)
#' kplot(flips, 2, 3)
#' @export
kplot <- function(flips, s, t, bw=FALSE) {
  k <- path <- num <- NULL
  if (!is.list(flips)) {
    d <- flips_to_kplot_df(flips)
    if (bw) {
      p <- qplot(k, path, data=d, geom="line") +
        scale_x_continuous(breaks=0:(t + s), limits = c(0, t + s)) +
        scale_y_continuous(breaks=0:s, limits=c(0, s+0.15)) +
        geom_segment(x=s, y=s, xend=(t+s-1), yend=s, linetype=2) +
        geom_segment(x=t, y=0, xend=(s+t-1), yend=s-1, linetype=2)
    } else {
      p <- qplot(k, path, data=d, geom="line") +
        scale_x_continuous(breaks=0:(t + s), limits=c(0, t + s)) +
        scale_y_continuous(breaks=0:s, limits=c(0, s+0.15)) +
        geom_segment(x=s,y=s,xend=(t+s-1),yend=s, color="green", linetype=1) +
        geom_segment(x=t, y=0, xend=(s+t-1), yend=s-1, col="red")
    }                            
  } else { 
    flip_set <- lapply(flips, flips_to_kplot_df)
    for (i in seq_along(flip_set)) {
      flip_set[[i]]$num <- as.factor(i)
      flip_set[[i]]$k <- jitter(flip_set[[i]]$k)
      flip_set[[i]]$k[flip_set[[i]]$k < 0] <- 0
    }
    d <- Reduce(rbind, flip_set)[, -(4:5)]
    if (bw) {
      p <- qplot(k, path, data=d, geom="path", group=num) +
        scale_x_continuous(breaks=0:(t+s), limits=c(0, t+s)) +
        geom_segment(x=s, y=s, xend=(t + s - 1), yend=s,
                     linetype=2) +
        geom_segment(x=t, y=0, xend=(s+t-1), yend=s-1, linetype=2)
    } else {
      p <- qplot(k, path, data=d, geom="path", group=num) +
        scale_x_continuous(breaks=0:(t+s), limits = c(0, t+s)) +
        geom_segment(x=s, y=s, xend=(t + s - 1), yend=s,
                     color="green") +
        geom_segment(x=t, y=0, xend=(s+t-1), yend=s-1, col="red")
    }     
  }   
  p   
}     

