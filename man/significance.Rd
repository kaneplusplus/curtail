% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two-stage-design.R
\name{significance}
\alias{significance}
\title{Significance of a two-stage curtail trial}
\usage{
significance(x)
}
\arguments{
\item{x}{a two-stage curtailed trial object}
}
\description{
This function returns the significance (probability of rejecting the 
null hypothesis when the null hypothesis is true) of all designs contained in a 
two-stage curtail trial object
}
\examples{
trials <- two_stage_curtail_trial(p1_null = 0.8, p2_null=0.2, 
p1_alt = 0.8, p2_alt = 0.4, n_total=36, prob_early=0.2, alpha=0.05)
significance(trials)
}
