#' @title One- and Two-Stage Curtailed Designs Using the Stopped Negative Binomial Distribution
#' @name curtail
#' @docType package
#' @description The \pkg{curtail} package is used for the planning of one- or
#' two-stage designs with curtailed sampling.  These designs make use of the Stopped Negative 
#' Binomial distribution, which characterizes the distribution of the sample size under 
#' curtailed sampling.  Under curtailed samplng, an early decision in the trial is allowed as 
#' soon as any predefined endpoint is reached.  The Stopped Negative Binomial Distribution
#' describes the smallest number of independent and indentically distributed Bernoulli trials
#' needed to observe a given number of successes or a given number of failures.  This package 
#' also contains density and other related functions of the Stopped Negative Binomial 
#' Distribution.
#' 
#' In the one-stage design, patients are assumed to be enrolled into the trial sequentially
#' until a maximum number of patients has been enrolled.  A critical value (threshold) of 
#' observed patient successes to deem the therapy superior is set prior to the start of the 
#' study.  This threshold is set to maintain a desired significance level based on the 
#' Bernoulli probability of patient success under the null hypothesis.  Under curtailed 
#' sampling,the study ends as soon as any predefined endpoint is reach.  In this case, 
#' patients will continue to be enrolled into the study until the observed number of patient 
#' successes exceeds the set critical value, or until the number of patient failures is too 
#' great for us to be able to deem the therapy superior based on maximum number of patients 
#' planned for the study.  The \code{single_stage_significance} and \code{single_stage_power} 
#' functions can be used to aid in the planning of the one-stage design.
#' 
#' A maximum number of patients are allocated to Stage~1 and Stage~2. Patients are
#' 
#' Stopped Negative Binomial Distribution Function Calls
#' \itemize{
#' \item{dsnb_stacked: }{}
#' \item{stacked_plot: }{}
#' \item{dsnb_plot: }{}
#' \item{dsnb_stack_plot: }{}
#' \item{cdsnb_stacked: }{}
#' \item{dsnb: }{}
#' \item{psnb: }{}
#' \item{qsnb: }{}
#' \item{rsnb: }{}
#' \item{esnb: }{}
#' \item{vsnb: }{}
#' \item{ecsnb: }{}
#' \item{vcsnb: }{}
#' \item{zplot: }{}
#' \item{kplot: }{} 
#' 
#' }
#' 
#' One-Stage Design Function Calls
#' \itemize{
#' \item{critical_values: }{}
#' \item{single_stage_significance: }{}
#' \item{single_stage_power: }{}
#' \item{expected_stage1_sample_size: }{}
#' \item{zplot: }{}
#' \item{kplot: }{}
#' \item{power_significance_plot: }{}
#' \item{power_significance_ROC: }{}
#' }
#' 
#' Two-Stage Design Function Calls
#' \itemize{
#' \item{critical_values: }{}
#' \item{prob_early_stop: }{}
#' \item{expected_stage1_sample_size: }{}
#' \item{expected_total_sample_size: }{}
#' \item{all_minimax_designs: }{}
#' \item{all_optimal_designs: }{}
#' \item{best_designs: }{}
#' \item{plot.ph2_design: }{}
#' \item{minimax_design: }{}
NULL