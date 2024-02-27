#' Internal functions for pk analysis
#' 
#' code{design_matrix}: design matrix for given sequence
#'
#' code{select_theta}: select good starting values for theta
#'
#' code{setup_env}: create environment with required variables
#'
#' code{sigma_vals}: calculate sigma values used in varcov matrix
#'
#' code{validateColumns}: ensure accurate column specification
#'
#' code{varcov_matrix}: design varcov matrix
#'
#' @name be-internal
#' @aliases design_matrix select_theta setup_env
#' sigma_vals validateColumns varcov_matrix
#' @keywords internal
NULL

# taken from EHR::validateColumns
validateColumns <- function(df, columnSpecs, defaultSpecs = list()) {
  # KEY = colname(s)
  # if default is NULL, not required
  # if default is NA, required
  n <- names(df)
  u <- names(columnSpecs)
  if(length(defaultSpecs) == 0) {
    defaultSpecs <- as.list(rep(NA, length(columnSpecs)))
    names(defaultSpecs) <- u
  }
  x <- names(defaultSpecs)
  cst <- function(t, v) {
    sprintf('%s"%s"', t, paste(v, collapse = '", "'))
  }
  errors <- character(4)
  # user should provide named list, or list of given length
  if(is.null(u)) {
    if(length(columnSpecs) == length(x)) {
      u <- x
      names(columnSpecs) <- x
    } else {
      errors[1] <- cst('column specification is incorrect; please identify all columns: ', x)
    }
  }
  # user should not provide unexpected columns
  bad_col <- setdiff(u, x)
  if(length(bad_col)) {
    errors[2] <- cst('column specification is incorrect; the following column(s) should not be present: ', bad_col)
  }
  # provide defaults, including NULL/NA
  add_col <- setdiff(x, u)
  columnSpecs[add_col] <- defaultSpecs[add_col]
  # safely remove NULL
  columnSpecs <- columnSpecs[lengths(columnSpecs, FALSE) > 0L]
  u <- names(columnSpecs)
  # require NA
  req_col <- u[is.na(sapply(columnSpecs, `[`, 1))]
  if(length(req_col)) {
    errors[3] <- cst('column specification is incorrect; please identify the following columns: ', req_col)
  }
  # check for missing columns
  mycols <- unlist(columnSpecs)
  mycols <- mycols[!is.na(mycols)]
  miss_col <- setdiff(mycols, c(n, seq_along(n)))
  if(length(miss_col)) {
    errors[4] <- cst('data set is missing expected columns; the following column(s) are missing: ', miss_col)
  }
  err <- paste(errors[errors != ''], collapse = '\n  ')
  if(err != '') stop(err)
  # convert any numeric columns into names
  for(i in seq_along(columnSpecs)) {
    csix <- match(columnSpecs[[i]], seq_along(n))
    columnSpecs[[i]][!is.na(csix)] <- n[csix[!is.na(csix)]]
  }
  columnSpecs
}

setup_env <- function(dat, colSpec = list()) {
    dat.req <- list(subject = NA, formula = NA, y = NA, period = NULL, seq = NULL)
    dat.col <- validateColumns(dat, colSpec, dat.req)

    if('period' %in% names(colSpec)) {
        # need to reorder by subject|period
        dat <- dat[order(dat[,colSpec$subject], dat[,colSpec$period]),]
        e <- list2env(list(subject = dat[,colSpec$subject], formula = dat[,colSpec$formula]))
        id_rows <- tapply(seq(nrow(dat)), e$subject, I)
        e$period <- dat[,colSpec$period]
    } else {
        e <- list2env(list(subject = dat[,colSpec$subject], formula = dat[,colSpec$formula]))
        id_rows <- tapply(seq(nrow(dat)), e$subject, I)
        e$period <- unsplit(lapply(id_rows, seq_along), e$subject)
    }
    e$Y <- log(dat[,colSpec$y]) ##take log of outcome

    ind <- match(e$formula, c('R','T'))
    e$Rind <- c(1,0)[ind]
    e$Tind <- ind-1

    TRseq <- vapply(id_rows, function(i) paste(e$Tind[i], collapse = ''), character(1))
    uTRseq <- unname(unique(TRseq))
    genseq <- unsplit(match(TRseq, uTRseq), e$subject)
    if('seq' %in% names(colSpec)) {
        ## this needs a solution
        provseq <- dat[,colSpec$seq]
        if(any(provseq != genseq)) {
            uTRseq <- uTRseq[c(2,1)]
            genseq <- unsplit(match(TRseq, uTRseq), e$subject)
        }
        if(any(provseq != genseq)) {
            stop('sequence provided does not match generated')
        }
        e$seq <- provseq
    } else {
        e$seq <- genseq
    }
    e$TRname <- chartr('01', 'RT', uTRseq)
    e$TRnum <- lapply(strsplit(uTRseq, ''), as.numeric)
    e
}

select_theta <- function(e, nPeriods) {
    if(!requireNamespace("nlme", quietly = TRUE)) {
        stop('please install "nlme" package -- install.packages(\'nlme\')')
    }

    ###########
    f <- stats::as.formula(Y ~ factor(seq)+factor(period)+formula)
    m <- nlme::lme(f, data = e, method="REML", random=list(~0+Rind+Tind|subject), weights= nlme::varIdent(form= ~1|formula), na.action = stats::na.exclude)
    ##########

    ##extract ratio of swt/swr
    devratio <- 1/unique(nlme::varWeights(m$modelStruct))[2]
    if(is.na(devratio)) devratio <- 1

    ##extract sbr and sbt and correlation coefficient etc.##
    vcorr <- nlme::VarCorr(m)
    sbr2 <- as.numeric(vcorr[1,1])
    sbt2 <- as.numeric(vcorr[2,1])
    rho <- as.numeric(vcorr[2,3])
    swr2 <- as.numeric(vcorr[3,1])
    swr <- as.numeric(vcorr[3,2])
    swt2 <- (swr*devratio)^2

    ##extract the fixed effects##
    sumTT <- summary(m)$tTable
    sumTT1 <- sumTT[seq(nPeriods + 1), 1]
    # theta1 (mu, p2, p3, p4, S, phi)
    theta1 <- unname(c(sumTT1[-2], sumTT1[2], sumTT[nPeriods+2]))
    c(theta1, log(sbt2), log(sbr2), log(swt2), log(swr2), rho)
}

sigma_vals <- function(theta, method = c('average','total','within'), nPeriods, val) {
    meth <- match.arg(method)
    # theta
    ## seq1: mu, p2, p3, p4, S, phi
    ## seq2: logsigmaBT2, logsigmaBR2, logsigmaWT2, logsigmaWR2
    ## seq3: rho
    seq1 <- seq(nPeriods+2)
    seq2 <- seq(nPeriods+3, length.out = 4)
    seq3 <- nPeriods + 7
    beta <- theta[seq1]
    rho <- theta[seq3]
    sigma2 <- exp(theta[seq2])
    bt <- sigma2[1]
    br <- sigma2[2]
    wt <- sigma2[3]
    wr <- sigma2[4]
    tt <- bt + wt
    tr <- br + wr
    s_vals <- c(
        BT = bt,
        BR = br,
        TT = tt,
        TR = tr,
        BRT = NA
    )
    if(method == 'average') {
        beta[nPeriods+2] <- val
    } else if(method == 'total') {
        s_vals['TT'] <- val^2 * tr
        s_vals['BT'] <- s_vals['TT'] - wt
        # if "px" is small, this could be negative (and a problem)
#         if(s_vals['BT'] < 0) warning('possible "xlow" issue')
    } else if(method == 'within') {
        s_vals['TT'] <- bt + val^2 * wr
    }
    s_vals['BRT'] <- rho * prod(sqrt(s_vals[c('BT','BR')]))
    list(beta, s_vals)
}

design_matrix <- function(x, seqno = 1) {
    n <- length(x)
    m <- cbind(diag(n), seqno-1, x)
    m[,1] <- 1
    colnames(m) <- c('mu', paste0('p', seq(2, n)), 'S', 'phi')
    m
}

varcov_matrix <- function(x) {
    RTname <- strsplit(x, '')[[1]]
    n <- length(RTname)
    mmm <- matrix('', n, n)
    diag(mmm) <- paste0('T', RTname)

    j1 <- matrix(RTname[which(upper.tri(mmm), arr.ind = TRUE)], ncol = 2)
    j2 <- matrix(c('BR', 'BRT', 'BRT', 'BT'), 2, 2, dimnames = list(c('R','T'), c('R','T')))
    j3 <- j2[j1]
    mmm[upper.tri(mmm)] <- j3
    mmm[lower.tri(mmm)] <- j3
    mmm
}
