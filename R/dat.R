#' Data example for bioequivalence (BE) study
#'
#' The dataset is a bioequivalence dataset from a fully repicated 2x4 crossover design with RTRT and TRTR as sequences. 
#' It is a subset of Example 4.4 in Chapter 4 of Patterson and Jones's book.
#'
#' @format
#' A data frame with 176 observations (from 44 subjects) on 6 variables:
#' \describe{
#'   \item{subject}{subject ID}
#'   \item{sequence}{RTRT or TRTR, where T and R stand for test and reference drugs, respectively}
#'   \item{period}{1 to 4 for crossover period}
#'   \item{formula}{T or R stand for test and reference drugs, respectively}
#'   \item{AUC}{a pharmacokinetic parameter - the area under the blood/plasma concentration-time curve}
#'   \item{CMAX}{a pharmacokinetic parameter - the peak concentration}
#' }
#'
#' @usage data(dat, package = 'BElikelihood')
#'
#' @keywords dataset
#'
#' @source Patterson S and Jones B (2023). Bioequivalence and Statistics in Clinical Pharmacology. Chapman Hall/CRC Press.
#' 
#' @examples
#' data(dat)
"dat"
