#' Likelihood method for analyzing Bioequivalence
#'
#' This package will calculate and plot the profile likelihoods for the mean difference and standard deviation ratios of a test drug to a reference drug for
#' AUC or Cmax from a full replicate 2x4 (RTRT/TRTR) cross-over, 2x3 partial replicate (TRT/RTR) cross-over, or a 2x2 (RT/TR) cross-over design (where R and T for reference and test drugs, respectively). 
#'
#' @docType package
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats as.formula na.exclude nlm
#' @import ggplot2
"_PACKAGE"
