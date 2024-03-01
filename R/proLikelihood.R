#' Calculate profile likelihood for bioequivalence data
#'
#' This is a general function to calculate the profile likelihoods for the mean difference, total standard deviation ratio, 
#' and within-subject standard deviation ratio of the test drug to the reference drug from bioequivalence (BE) study data. 
#' Standardized profile likelihood plots with the 1/8 and 1/32 likelihood intervals can be generated using the plot method. 
#' The within-subject standard deviation ratio can be obtained only for a fully replicated 2x4 or a partially replicated 2x3 design.  
#' 
#' @details This function implements a likelihood method for evaluating BE for pharmacokinetic parameters (AUC and Cmax) (see reference below). It accepts a dataframe collected with various crossover designs commonly used in BE studies such as 
#' a fully replicated crossover design (e.g., 2x4 two-sequence, four-period, RTRT/TRTR), a partially replicated crossover design 
#' (e.g., 2x3, two-sequence, three-period, RTR/TRT), and a two-sequence, two-period, crossover design design (2x2, RT/TR), 
#' where "R" stands for a reference drug and "T" stands for a test drug. 
#' It allows missing data (for example, a subject may miss the period 2 data) and utilizes all available data. It will 
#' calculate the profile likelihoods for the mean difference, total standard deviation ratio, and within-subject standard deviation ratio. 
#' Plots of standardized profile likelihood can be generated and provide evidence for various quantities of interest for evaluating 
#' BE in a unified framework.
#' 
#' @param dat data frame contains BE data (AUC and Cmax) with missing data allowed.
#' @param colSpec a named list that should specify columns in \sQuote{dat}; \sQuote{subject} (subject ID),
#' \sQuote{formula} (must be coded as T or R, where T for test drug and R for reference drug), and \sQuote{y} (either AUC or Cmax) are
#' required. \sQuote{period} and \sQuote{seq} may also be provided.
#' The \sQuote{formula} column should identify a test or a reference drug with R and T.
#' @param theta An optional numeric vector contains initial values of the parameters for use in optimization.
#' For example, in a 2x4 fully replicated design, the vector is [mu, p2, p3, p4, S, phi,log(sbt2), log(sbr2), log(swt2), log(sbr2), rho], where  
#' \sQuote{mu} is the population mean for the reference drug when there are no period or sequence effects; \sQuote{p2} to \sQuote{p4} are fixed 
#' period effects with period 1 as the reference period; \sQuote{S} the fixed sequence effect with seq 1 as the reference sequence; \sQuote{phi} 
#' is the mean difference between the two drugs; \sQuote{sbt2} and \sQuote{sbr2} are between-subject variances for the test and reference drugs, 
#' respectively; \sQuote{swt2} and \sQuote{swr2} are within-subject variances for the test and reference drugs, respectively; \sQuote{rho} is 
#' the correlation coefficient within a subject. When \sQuote{theta} (default is null) is not provided, the function 
#' will choose the starting values automatically based on a linear mixed-effects model. If users want to provide these values, for method 
#' \sQuote{average} (mean difference), user may put any value for \sQuote{phi}. Similarly, for method \sQuote{total}, user can put any value 
#' for \sQuote{log(sbt2)}, and for method \sQuote{within}, user can put any value for \sQuote{log(swt2)}.
#' @param xlow numeric value, the lower limit of x-axis for the profile likelihood plot, at which the profile likelihood is calculated. It is 
#' optional and can be automatically generated using the maximum likelihood estimate (MLE) depending on the \sQuote{method}. We strongly 
#' recommend users trying a better value that would better fit for purpose. 
#' @param xup numeric value, the upper limit of x-axis for the profile likelihood plot, at which the profile likelihood is calculated. It is 
#' optional and can be automatically generated using the MLE depending on the \sQuote{method}. We strongly recommend users trying 
#' a better value that would better fit for purpose. 
#' @param xlength numeric value. Defaults to 100. It is the number of grids between the lower and upper limits, which controls smoothness of 
#' the curve. It will take longer time to run for larger number of grids, but we strongly recommend users using a larger number than the default 
#' value.
#' @param method character value. Should be one of \sQuote{average}, \sQuote{total}, or \sQuote{within}. 
#' \sQuote{average} will provide the profile likelihood for the mean difference between test and reference drugs. 
#' \sQuote{total} will provide the profile likelihood for the total standard deviation ratio of test to reference drug. \sQuote{within} 
#' will provide the profile likelihood for the within-subject standard deviation ratio of test to reference drug when appropriate. 
#' 
#' @return A \sQuote{proLikelihood} object, with elements \sQuote{poi}, \sQuote{maxLik}, \sQuote{MAX}, \sQuote{LI}, and \sQuote{method}. 
#' \sQuote{poi} and \sQuote{maxLik} are the interested parameter (mean difference, total standard deviation ratio 
#' or within-subject standard deviation ratio) values and the corresponding profile likelihood values, respectively. \sQuote{MAX} is the MLE 
#' estimate for that parameter. \sQuote{LI} is the likelihood intervals with the 1/4.5, 1/8 and 1/32 intervals. 
#' \sQuote{method} is one of \sQuote{average},\sQuote{total}, and \sQuote{within}. 
#' 
#' @references Liping Du and Leena Choi, Likelihood approach for evaluating bioequivalence of highly variable drugs, Pharmaceutical Statistics, 14(2): 82-94, 2015
#'
#' @examples
#' \donttest{
#' data(dat)
#' cols <- list(subject = 'subject', formula = 'formula', y = 'AUC')
#' p4a <- proLikelihood(dat, colSpec = cols, xlength = 300, method = 'average')
#' p4t <- proLikelihood(dat, colSpec = cols, xlength = 300, method = 'total')
#' p4w <- proLikelihood(dat, colSpec = cols, xlength = 300, method = 'within')
#' # three period case
#' dd3 <- dat[dat$period < 4,]
#' p3a <- averageBE(dd3, colSpec = cols, xlength = 300)
#' p3t <- totalVarianceBE(dd3, colSpec = cols, xlength = 300)
#' p3w <- withinVarianceBE(dd3, colSpec = cols, xlength = 300)
#' # two period case
#' dd2 <- dat[dat$period < 3,]
#' p2a <- averageBE(dd2, colSpec = cols, xlength = 300)
#' p2t <- totalVarianceBE(dd2, colSpec = cols, xlength = 300)
#' }
#'
#' @export

proLikelihood <- function(dat, colSpec = list(), theta = NULL, xlow, xup, xlength = 100, method) {
  m <- match.arg(method, c('average','total','within'))

  e <- setup_env(dat, colSpec)
  TRname <- e$TRname
  TRnum <- e$TRnum
  nSeq <- max(e$seq)
  nPeriods <- max(e$period)
  if(nSeq != 2) stop('data must contain two sequences')
  if(m == 'within' && (nPeriods < 3 || nPeriods > 4)) stop('data must contain 3-4 periods')
  if(nPeriods < 2 || nPeriods > 4) stop('data must contain 2-4 periods')

  subject1 <- unique(e$subject[e$seq == 1]) ## unique subjects in seq 1
  subject2 <- unique(e$subject[e$seq == 2]) ## unique subjects in seq 2
  n1 <- length(subject1) ## number of subjects in seq 1
  n2 <- length(subject2) ## number of subjects in seq 2

  ###get the starting value for theta if not provided###
  if(is.null(theta)) theta <- select_theta(e, nPeriods)
  expThetaSize <- nPeriods + 4 + 3 # period 4?x, logsigma 4x, S, phi, rho
  if(length(theta) < expThetaSize) {
    stop('the specified theta has too few values')
  }
  if(missing(xlow) && missing(xup)) {
    if(m == 'average') {
      spot <- theta[nPeriods + 2]
      xlow <- min(-0.225, spot - 0.223)
      xup <- max(0.225, spot + 0.223)
    } else {
      seq2 <- seq(nPeriods+3, length.out = 4)
      sigma2 <- exp(theta[seq2])
      bt <- sigma2[1]
      br <- sigma2[2]
      wt <- sigma2[3]
      wr <- sigma2[4]
      tt <- bt + wt
      tr <- br + wr
      if(m == 'total') {
        spot <- tt / tr
      } else if(m == 'within') {
        spot <- wt / wr
      }
      xlow <- min(0.7, spot * 0.7)
      xup <- max(1.3, spot * 1.3)
    }
  }
  x <- seq(xlow, xup, length.out = xlength)##fixed phi values

  ## design matrix for TRRT and RTTR design (mu, p2, p3, p4, S, phi)
  X <- lapply(seq_along(TRnum), function(i) {
    design_matrix(TRnum[[i]], i)
  })
  names(X) <- TRname

  s1xy <- lapply(seq(n1), function(i) {
    yi <- e$Y[e$subject==subject1[i]]
    miss.pos <- which(is.na(yi)) # missing position
    if(length(miss.pos) > 0){
      yi <- yi[-miss.pos]
      X.m <- X[[1]][-miss.pos,]
    } else{
      X.m <- X[[1]]
    }
    list(yi, miss.pos, X.m)
  })
  s2xy <- lapply(seq(n2), function(i) {
    yi <- e$Y[e$subject==subject2[i]]
    miss.pos <- which(is.na(yi)) # missing position
    if(length(miss.pos) > 0){
      yi <- yi[-miss.pos]
      X.m <- X[[2]][-miss.pos,]
    } else{
      X.m <- X[[2]]
    }
    list(yi, miss.pos, X.m)
  })

  # use TRname to generalize var-cov matrix
  varcov_blueprint <- lapply(TRname, varcov_matrix)
  names(varcov_blueprint) <- TRname

  ##negative log likelihood ###
  mnormNLL <- function(theta, val) {
    p <- sigma_vals(theta, method, nPeriods, val)
    beta <- p[[1]]
    s_vals <- p[[2]]

    ##var/cov matrix####
    vmat <- lapply(varcov_blueprint, function(i) {
      matrix(s_vals[i], nPeriods, nPeriods)
    })

    l <- 0 ##variable for the sum of negative log likelihood##

    for (i in seq(n1)){
      s_i <- s1xy[[i]]
      yi <- s_i[[1]]
      miss.pos <- s_i[[2]]
      Xi1.m <- s_i[[3]]
      if(length(miss.pos) > 0) {
        vi1.m <- vmat[[1]][-miss.pos, -miss.pos, drop = FALSE]
      } else {
        vi1.m <- vmat[[1]]
      }
      y.pred1 <- Xi1.m %*% beta
      l <- l - mvtnorm::dmvnorm(yi, mean=y.pred1, sigma=vi1.m, log=TRUE, checkSymmetry = FALSE)
    } ## sum of negative log likelihood for subjects in seq 1

    for(i in seq(n2)){
      s_i <- s2xy[[i]]
      yi <- s_i[[1]]
      miss.pos <- s_i[[2]]
      Xi2.m <- s_i[[3]]
      if(length(miss.pos) > 0) {
        vi2.m <- vmat[[2]][-miss.pos, -miss.pos, drop = FALSE]
      } else {
        vi2.m <- vmat[[2]]
      }
      y.pred2 <- Xi2.m %*% beta
      l <- l - mvtnorm::dmvnorm(yi, mean=y.pred2, sigma=vi2.m, log=TRUE, checkSymmetry = FALSE)
    } ## sum of log likelihood for subjects in seq 2
    return(l)
  }

  #################get the profilelikelihood################

  ###get the profile likelihood for fixed phi##
  maxLik <- rep(NA, xlength)

  for(i in seq(xlength)) {
    ## get minimized negative log likelihood using nlm###
    op <- stats::nlm(mnormNLL, theta, x[i])
    # would `optim` work here?

    ## convert to maximum likelihood##
    if (op$code <= 2) {
      maxLik[i]<- exp(-op$minimum)
    }
  }
  lik.norm <- maxLik / max(maxLik, na.rm = TRUE)
  xmax <- x[which.max(lik.norm)]
  xli4_5 <- range(x[lik.norm >=1/4.5], na.rm = TRUE)
  xli8 <- range(x[lik.norm >=1/8], na.rm = TRUE)
  xli32 <- range(x[lik.norm >=1/32], na.rm = TRUE)
  li <- rbind(xli4_5, xli8, xli32)
  rownames(li) <- c('1/4.5 LI', '1/8 LI', '1/32 LI')
  colnames(li) <- c('lower', 'upper')
  obj <- list(poi = x, maxLik = maxLik, MAX = xmax, LI = li, method = m)
  class(obj) <- 'proLikelihood'
  obj
}

#' @rdname proLikelihood
#' @export

averageBE <- function(dat, colSpec = list(), theta = NULL, xlow, xup, xlength) {
  proLikelihood(dat, colSpec, theta, xlow, xup, xlength, 'average')
}

#' @rdname proLikelihood
#' @export

totalVarianceBE <- function(dat, colSpec = list(), theta = NULL, xlow, xup, xlength) {
  proLikelihood(dat, colSpec, theta, xlow, xup, xlength, 'total')
}

#' @rdname proLikelihood
#' @export

withinVarianceBE <- function(dat, colSpec = list(), theta = NULL, xlow, xup, xlength) {
  proLikelihood(dat, colSpec, theta, xlow, xup, xlength, 'within')
}
