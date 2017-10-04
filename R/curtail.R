#' @title Design of One- and Two-Stage Clinical Trials with Curtailed Sampling
#' @name curtail
#' @docType package
#' @description The \pkg{curtail} package is used for the planning of one- and two-
#' stage clinical trials with curtailed sampling.  Under curtailed sampling, an early 
#' decision in the trial is allowed as soon as any predefined statistical endpoint is 
#' reached.  The package provides functions to help select a design, including 
#' visualizations to compare criteria for different choices of parameters, and functions 
#' to calculate power, significance, and expected sample size, among others.
#' 
#' In the one-stage design, patients are assumed to be enrolled into the trial 
#' sequentially up to a maximum number of patients.  A critical value of observed 
#' patient successes needed to deem the therapy superior is set prior to the start of
#' the study.  The \code{power_significance_plot} and \code{power_significance_ROC} 
#' functions provide visualizations to compare power and significance levels for 
#' various choices of critical values.  Also, the \code{critical_values} function 
#' will calculate the critical value to maintain a desired significance level.  
#' Under curtailed sampling, the study ends as soon as enough the observed patient 
#' successes meets the critical value or as soon as too many patient failures have 
#' been observed.  The smallest number of patient enrollees needed to reach a 
#' decision under curtailed sampling is described by the Stopped Negative Binomial 
#' distribution.  This package also provides density and other related functions for 
#' the Stopped Negative Binomial Distribution.  
#' 
#' The two-stage design presented in this package is a modification of Simon's two-
#' stage design with separate, but nested, criteria for early stopping in Stage 1 
#' and efficacy in Stage 2.  The two-stage design has two critical values which can 
#' both be determined with the \code{critical_values} function to maintain a desired 
#' significance level and probability of early stopping.  These critical values are set 
#' prior to the start of the study to determine the number of patient successes in 
#' the first stage needed to continue the trial to the second stage and the critical 
#' number of efficacy successes throughout the trial needed to deem the therapy 
#' superior.  Under curtailed sampling, early decisions can be made in Stage 1 and 
#' Stage 2.  The \code{best_designs} function finds the optimal and minimax design 
#' for a fixed total sample size.  Other functions are provided to calculate the 
#' expected sample size, power, significance, and probability of early stopping in 
#' the trial given a choice of design parameters.
#' 
#' 
#' One-Stage Design Function Calls
#' \itemize{
#' \item{critical_values: }{finds the minimum number of successes to reject the null 
#' hypothesis with significance level alpha}
#' \item{single_stage_significance: }{computes the probability of rejecting the null 
#' hypothesis assuming the null probability of success}
#' \item{single_stage_power: }{computes the probability of rejecting the null 
#' hypothesis assuming an alternative probability of success}
#' \item{expected_stage1_sample_size: }{computes the mean and standard deviation of 
#' the sample size for the one-stage design}
#' \item{zplot: }{visualize the Stopped Negative Binomial process with horizontal 
#' axis counting patient successes and vertical axis counting patient failures}
#' \item{kplot: }{visualize the Stopped Negative Binomial process with horizontal 
#' axis counting patients enrolled and vertical axis counting the number of 
#' successes}
#' \item{power_significance_plot: }{Plot power and significance across all trial 
#' designs with a fixed maximum number of patients}
#' \item{power_significance_ROC: }{ROC curve of Power vs. 1-Significance for all 
#' trial designs with a fixed number of maximum patients}
#' }
#' 
#' #' Stopped Negative Binomial Distribution Function Calls
#' \itemize{
#' \item{dsnb: }{density for the Stopped Negative Binomial Distribution}
#' \item{psnb: }{distribution function for the Stopped Negative Binomial Distribution}
#' \item{qsnb: }{quantile function for the Stopped Negative Binomial Distribution}
#' \item{rsnb: }{randomly generated value from the Stopped Negative Binomial 
#' distribution}
#' \item{esnb: }{expected value of the Stopped Negative Binomial distribution}
#' \item{vsnb: }{variance of the Stopped Negative Binomial distribution}
#' \item{ecsnb: }{expected value of the Stopped Negative Binomial distribution}
#' \item{vcsnb: }{variance of the Stopped Negative Binomial distribution}
#' \item{dsnb_stacked: }{}
#' \item{stacked_plot: }{}
#' \item{dsnb_plot: }{}
#' \item{dsnb_stack_plot: }{}
#' \item{cdsnb_stacked: }{}
#' }
#' 
#' Two-Stage Design Function Calls
#' \itemize{
#' \item{critical_values: }{finds the minimum number of successes to continue the 
#' trial to the Stage 2 and the minimum number of successes to reject the null 
#' hypothesis with a given significance level}
#' \item{two_stage_significance: }{calculates the probability of rejecting the null 
#' hypothesis assuming null probabililties of success for Stage 1 and Stage 2}
#' \item{two_stage_power: }{calculates the probability of rejecting the null 
#' hypothesis assuming alternative probabililties of success for Stage 1 and Stage 
#' 2}
#' \item{prob_early_stop: }{computes the probability of stopping early for futility 
#' in the two-stage design}
#' \item{expected_stage1_sample_size: }{computes the mean and standard deviation of 
#' the sample size for Stage 1 of a two-stage design}
#' \item{expected_total_sample_size: }{computes the expected sample size for a given
#' two-stage design}
#' \item{minimax_design: }{computes the minimax probability for a given two-stage 
#' design}
#' \item{all_minimax_designs: }{computes the minimax probability for each possible 
#' design for a fixed total sample size}
#' \item{all_optimal_designs: }{computes the expected total sample size for each
#' possible design for a fixed total sample size}
#' \item{best_designs: }{finds the optimal and minimax design for a fixed total sample 
#' size }
#' \item{plot.ph2_design: }{plots the optimal and minimax criteria for all possible 
#' designs with a fixed total sample size}
#' 
#' 
#' Acknowledgements: This work was partially supported through a 
#' Patient-Centered Outcomes Research Institute (PCORI) Award (ME-1511-32832).
#'
#' Disclaimer: All statements in this report, including its findings and 
#' conclusions, are solely those of the authors and do not necessarily 
#' represent the views of the Patient-Centered Outcomes Research Institute 
#' (PCORI), its Board of Governors or Methodology Committee.
#' 
#' @references Rucker G (1989). "A two-stage trial design for testing treatment,
#' self-selection and treatment preference effects." \emph{Stat Med}, 
#' \strong{8}(4):477-485. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/2727471}{PubMed})
NULL