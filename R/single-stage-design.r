#' Generic function for creating single-stage curtail trials
#' @export
setGeneric("single_stage_curtail_trial", function(p_null, p_alt, s, t, n) {
  standardGeneric("single_stage_curtail_trial")
})

# @examples
# trial <- single_stage_curtail_trial(p_null=0.2, p_alt=0.4, s=7, t=11)
#' @export
setMethod("single_stage_curtail_trial",
  signature(p_null="numeric", p_alt="numeric", s="numeric", t="numeric",
            n="missing"),
  function(p_null, p_alt, s, t) {
    ret <- data.frame(p_null=p_null, p_alt=p_alt, s=s, t=t)
    class(ret) <- c("single_stage_curtail_trial", class(ret))
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

#' @export
significance <- function(x) {
  UseMethod("significance")
}

#' @export
significance.single_stage_curtail_trial <- function(x) {
  apply(x, 1, 
    function(x) {
      sum(dsnb_stacked(min(x[c("s", "t")]):(sum(x[c("s", "t")])-1),
        p=x["p_null"], s=x["s"], t=x["t"])[,"s"])
    })
}

#' @export
power <- function(x) {
  UseMethod("power")
}

#' @export
power.single_stage_curtail_trial <- function(x) {
  apply(x, 1,
    function(x) {
      sum(dsnb_stacked(min(x[c("s", "t")]):(sum(x[c("s", "t")])-1), 
        p=x["p_alt"], s=x["s"], t=x["t"])[,'s'])
    })
}

#' @export
expected_sample_size <- function(x) {
  UseMethod("expected_sample_size")
}

#' @export
expected_sample_size.single_stage_curtail_trial <- function(x) {
  foreach(i=1:nrow(x), .combine=rbind) %do% {
    mean <- sum((1:(x$s[i]+x$t[i]-1)) * ph2snb(x$p_null[i], x$s[i], x$t[i]))
    sd <- sqrt(sum((1:(x$s[i]+x$t[i]-1))^2 * 
      ph2snb(x$p_null[i], x$s[i], x$t[i])) - mean^2)
    data.frame(mean=mean, sd=sd)
  }
}

#' @export
critical_values <- function(x, alpha=0.1, p_early=0.1) {
  UseMethod("critical_values")
}

#' @export
critical_values.single_stage_curtail_trial <- 
  function(x, alpha=0.1, p_early=NULL) {
  if (!is.null(p_early)) {
    warning("p_early specified but not used.")
  }
  if ((p_early > 1) || (p_early < 0)) {
    stop("p_early must be between zero and one inclusive.")
  }
  if ((alpha > 1) || (alpha < 0)) {
    stop("alpha must be between zero and one inclusive.")
  }
  1 + qbinom(1 - alpha, x$s + x$t - 1, x$p_null)
}

#' @export
print.single_stage_curtail_trial <- function(object, ...) {
  if (nrow(object) == 1) {
    # It's a single, single stage curtail trial.
    cat("\n")
    cat(" Single-Stage Curtail Trial\n")
    cat("\n")
    cat("Null response rate: ", object$p_null[1], "\n")
    cat("Alternative response rate: ", object$p_alt[1], "\n")
    cat("Responses to stop the trial: ", object$s[1], "\n")
    cat("Non-responses to stop the trial: ", object$t[1], "\n")
    cat("Power: ", power(object), "\n")
    cat("Significance: ", significance(object), "\n")
    ess <- expected_sample_size(object)
    cat("Expected sample size: ", ess$mean, "\n")
    cat("Sample size standard deviation: ", ess$sd, "\n\n")
  }
}

