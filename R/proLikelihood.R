#' Run profile likelihood
#'
#' function summary
#'
#' function details
#'
#' @param dat data.frame
#' @param colSpec a named list that should specify columns in data; \sQuote{subject},
#' \sQuote{formula}, and \sQuote{y} are required. \sQuote{period} and \sQuote{seq}
#' may also be provided. The \sQuote{formula} column should identify treatment with R/T.
#' @param theta numeric vector DESCRIBE
#' @param xlow numeric value DESCRIBE
#' @param xup numeric value DESCRIBE
#' @param xlength numeric value Defaults to 100; DESCRIBE
#' @param method character value Should be one of \sQuote{average}, \sQuote{total},
#' or \sQuote{within}. DESCRIBE
#'
#' @return A \sQuote{proLikelihood} object, with elements \sQuote{poi},
#' \sQuote{maxLik}, \sQuote{MAX}, \sQuote{LI}, and \sQuote{method}.
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
#' p3a <- proLikelihood(dd3, colSpec = cols, xlength = 300, method = 'average')
#' p3t <- proLikelihood(dd3, colSpec = cols, xlength = 300, method = 'total')
#' p3w <- proLikelihood(dd3, colSpec = cols, xlength = 300, method = 'within')
#' # two period case
#' dd2 <- dat[dat$period < 3,]
#' p2a <- proLikelihood(dd2, colSpec = cols, xlength = 300, method = 'average')
#' p2t <- proLikelihood(dd2, colSpec = cols, xlength = 300, method = 'total')
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
  xmax <- x[lik.norm==1]
  xli4_5 <- range(x[lik.norm >=1/4.5])
  xli8 <- range(x[lik.norm >=1/8])
  xli32 <- range(x[lik.norm >=1/32])
  li <- rbind(xli4_5, xli8, xli32)
  rownames(li) <- c('1/4.5 LI', '1/8 LI', '1/32 LI')
  colnames(li) <- c('lower', 'upper')
  obj <- list(poi = x, maxLik = maxLik, MAX = xmax, LI = li, method = m)
  class(obj) <- 'proLikelihood'
  obj
}

#' Run profile likelihood
#'
#' function summary
#'
#' function details
#'
#' @param dat data.frame
#' @param colSpec a named list that should specify columns in data; \sQuote{subject},
#' \sQuote{formula}, and \sQuote{y} are required. \sQuote{period} and \sQuote{seq}
#' may also be provided. The \sQuote{formula} column should identify treatment with R/T.
#' @param theta numeric vector DESCRIBE
#' @param xlow numeric value DESCRIBE
#' @param xup numeric value DESCRIBE
#' @param xlength numeric value Defaults to 100; DESCRIBE
#'
#' @return A \sQuote{proLikelihood} object, with elements \sQuote{poi},
#' \sQuote{maxLik}, \sQuote{MAX}, \sQuote{LI}, and \sQuote{method}.
#'
#' @examples
#' \donttest{
#' data(dat)
#' cols <- list(subject = 'subject', formula = 'formula', y = 'AUC')
#' l <- averageBE(dat, colSpec = cols, xlength = 300)
#' }
#'
#' @export

averageBE <- function(dat, colSpec = list(), theta = NULL, xlow, xup, xlength) {
  proLikelihood(dat, colSpec, theta, xlow, xup, xlength, 'average')
}

#' Run profile likelihood
#'
#' function summary
#'
#' function details
#'
#' @param dat data.frame
#' @param colSpec a named list that should specify columns in data; \sQuote{subject},
#' \sQuote{formula}, and \sQuote{y} are required. \sQuote{period} and \sQuote{seq}
#' may also be provided. The \sQuote{formula} column should identify treatment with R/T.
#' @param theta numeric vector DESCRIBE
#' @param xlow numeric value DESCRIBE
#' @param xup numeric value DESCRIBE
#' @param xlength numeric value Defaults to 100; DESCRIBE
#'
#' @return A \sQuote{proLikelihood} object, with elements \sQuote{poi},
#' \sQuote{maxLik}, \sQuote{MAX}, \sQuote{LI}, and \sQuote{method}.
#'
#' @examples
#' \donttest{
#' data(dat)
#' cols <- list(subject = 'subject', formula = 'formula', y = 'AUC')
#' tv <- totalVarianceBE(dat, colSpec = cols, xlength = 300)
#' }
#'
#' @export

totalVarianceBE <- function(dat, colSpec = list(), theta = NULL, xlow, xup, xlength) {
  proLikelihood(dat, colSpec, theta, xlow, xup, xlength, 'total')
}

#' Run profile likelihood
#'
#' function summary
#'
#' function details
#'
#' @param dat data.frame
#' @param colSpec a named list that should specify columns in data; \sQuote{subject},
#' \sQuote{formula}, and \sQuote{y} are required. \sQuote{period} and \sQuote{seq}
#' may also be provided. The \sQuote{formula} column should identify treatment with R/T.
#' @param theta numeric vector DESCRIBE
#' @param xlow numeric value DESCRIBE
#' @param xup numeric value DESCRIBE
#' @param xlength numeric value Defaults to 100; DESCRIBE
#'
#' @return A \sQuote{proLikelihood} object, with elements \sQuote{poi},
#' \sQuote{maxLik}, \sQuote{MAX}, \sQuote{LI}, and \sQuote{method}.
#'
#' @examples
#' \donttest{
#' data(dat)
#' cols <- list(subject = 'subject', formula = 'formula', y = 'AUC')
#' wv <- withinVarianceBE(dat, colSpec = cols, xlength = 300)
#' }
#'
#' @export

withinVarianceBE <- function(dat, colSpec = list(), theta = NULL, xlow, xup, xlength) {
  proLikelihood(dat, colSpec, theta, xlow, xup, xlength, 'within')
}
