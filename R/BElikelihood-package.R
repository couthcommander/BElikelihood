#' Likelihood method for analyzing bioequivalence (BE) trial data
#'
#' This package will calculate and plot the profile likelihoods for the mean difference and standard deviation ratios of a test drug 
#' to a reference drug for AUC or Cmax from various crossover designs commonly used in BE studies, such as a fully replicated crossover 
#' design (e.g., 2x4 two-sequence, four-period, RTRT/TRTR), a partially replicated crossover design (e.g., 2x3, two-sequence, three-period, RTR/TRT), and a two-sequence, two-period, crossover design design (2x2, RT/TR), where "R" stands for a reference drug and "T" stands for a test drug.
#'
#' @docType package
#'
#' @importFrom stats as.formula na.exclude nlm
#' @import ggplot2
"_PACKAGE"
