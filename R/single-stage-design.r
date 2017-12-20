#' Create a single-stage curtailed trial design
#' @description The single_stage_curtail_trial function creates a single-stage 
#' trial object containing design parameters and statistical properties of 
#' one or more curtailed designs
#' @param p_null probability of success under the null hypothesis
#' @param p_alt probability of success under the alternative hypothesis
#' @param s number of successes to stop the trial
#' @param t number of failures to stop the trial
#' @param n total maximum planned sample size (n = s+t-1)
#' @details
#' There are two ways that the user can specify the parameters to create
#' a single-stage curtailed trial design, as defined in the two cases below.  
#' In case 1, the user provide sparameters values of s and t, in which case the trial
#' object is created for that one specific design.  In the second case, the user 
#' specifies an n rather than s and t and the trial object created contains all 
#' possible designs (values of s and t) for the given total sample size, n.
#' 
#' Case 1:  User specifies p_null, p_alt, s, t
#' 
#' Case 2:  User specifies p_null, p_alt, n 
#' 
#' @examples
#' Case 1:
#' trial <- single_stage_curtail_trial(p_null=0.2, p_alt=0.4, s=7, t=11)
#' 
#' Case 2:
#' trials <- single_stage_curtail_trial(p_null=0.2, p_alt=0.4, n=17)
#' @export
setGeneric("single_stage_curtail_trial", function(p_null, p_alt, s, t, n) {
  standardGeneric("single_stage_curtail_trial")
})

# trial <- single_stage_curtail_trial(p_null=0.2, p_alt=0.4, s=7, t=11)
#' @export
setMethod("single_stage_curtail_trial",
  signature(p_null="numeric", p_alt="numeric", s="numeric", t="numeric",
            n="missing"),
  function(p_null, p_alt, s, t) {
    ret <- data.frame(p_null=p_null, p_alt=p_alt, s=s, t=t)
    class(ret) <- c("single_stage_curtail_trial", class(ret))
    ret$power <- power(ret)
    ret$significance <- significance(ret)
    ess <- expected_sample_size(ret)
    ret$mean_ss <- ess$mean
    ret$sd <- ess$sd_ss
    ret
  })

#trials <- single_stage_curtail_trial(p_null=0.2, p_alt=0.4, n=17)
#' @export
setMethod("single_stage_curtail_trial",
  signature(p_null="numeric", p_alt="numeric", s="missing", t="missing",
            n="numeric"),
  function(p_null, p_alt, n) {
    si <- Power <- Significance <- s <- NULL
    si <- seq_len(n-1)
    ti <- n + 1 - si
    ret <- data.frame(p_null=rep(p_null, length(si)), 
      p_alt=rep(p_alt, length(si)), s=si, t=ti)
    class(ret) <- c("single_stage_curtail_trial", 
      "single_stage_curtail_trial_roc", class(ret))
    ret$power <- power(ret)
    ret$significance <- significance(ret)
    ess <- expected_sample_size(ret)
    ret$mean_ss <- ess$mean
    ret$sd_ss <- ess$sd
    ret
  })

#' @export
plot.single_stage_curtail_trial_roc <- function(x, y, ...) {
  x$Significance <- significance(x)
  x$Power <- power(x)
  x$ess <- apply(x, 1, function(x) esnb(x["p_null"], x["s"], x["t"]))
  n <- x$s[1] + x$t[1] - 1 
  data_pos <- x[1:(n-1), c("Power", "Significance", "s")]
  data_pos$Power <- data_pos$Power + 0.02
  data_pos$Significance <- data_pos$Significance - 0.02
  # stop labeling when power < .05 or significance>.95 when all_labels=FALSE
  data_pos[which(data_pos$Power < 0.05 | data_pos$Significance > 0.95), "s"] <- 
    ""
  ggplot(x, aes(x=Power, y=1-Significance)) +
    geom_line() +
    geom_text(data=data_pos, aes(x=Power, y=1-Significance, label=s)) +
    theme_bw()

}


significance <- function(x) {
  UseMethod("significance")
}


significance.single_stage_curtail_trial <- function(x) {
  apply(x, 1, 
    function(x) {
      sum(dsnb_stacked(min(x[c("s", "t")]):(sum(x[c("s", "t")])-1),
        p=x["p_null"], s=x["s"], t=x["t"])[,"s"])
    })
}


power <- function(x) {
  UseMethod("power")
}


power.single_stage_curtail_trial <- function(x) {
  apply(x, 1,
    function(x) {
      sum(dsnb_stacked(min(x[c("s", "t")]):(sum(x[c("s", "t")])-1), 
        p=x["p_alt"], s=x["s"], t=x["t"])[,'s'])
    })
}


expected_sample_size <- function(x) {
  UseMethod("expected_sample_size")
}


expected_sample_size.single_stage_curtail_trial <- function(x) {
  foreach(i=1:nrow(x), .combine=rbind) %do% {
    mean <- sum((1:(x$s[i]+x$t[i]-1)) * ph2snb(x$p_null[i], x$s[i], x$t[i]))
    sd <- sqrt(sum((1:(x$s[i]+x$t[i]-1))^2 * 
      ph2snb(x$p_null[i], x$s[i], x$t[i])) - mean^2)
    data.frame(mean=mean, sd=sd)
  }
}


critical_values <- function(x, alpha=0.1, p_early=0.1) {
  UseMethod("critical_values")
}

# The following isn't right. Take the individual parameters.

critical_values.single_stage_curtail_trial <- 
  function(x, alpha=0.1, p_early=NULL) {
  if (!is.null(p_early)) {
    warning("p_early specified but not used.")
  }
  if ((alpha > 1) || (alpha < 0)) {
    stop("alpha must be between zero and one inclusive.")
  }
  1 + qbinom(1 - alpha, x$s + x$t - 1, x$p_null)
}

#'
#print.single_stage_curtail_trial <- function(object, ...) {
#  if (nrow(object) == 1) {
#    ret <- object
#    # It's a single, single stage curtail trial.
#    cat("\n")
#    cat(" Single-Stage Curtail Trial\n")
#    cat("\n")
#    cat("Null response rate: ", object$p_null[1], "\n")
#    cat("Alternative response rate: ", object$p_alt[1], "\n")
#    cat("Responses to stop the trial: ", object$s[1], "\n")
#    cat("Non-responses to stop the trial: ", object$t[1], "\n")
#    cat("Power: ", object$power, "\n")
#    cat("Significance: ", object$significance, "\n")
#    cat("Expected sample size: ", object$mean, "\n")
#    cat("Sample size standard deviation: ", object$sd, "\n\n")
#  }
#}

