% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/printplot.R
\name{print.proLikelihood}
\alias{print.proLikelihood}
\title{Print method for proLikelihood object}
\usage{
\method{print}{proLikelihood}(x, ...)
}
\arguments{
\item{x}{proLikelihood object}

\item{\dots}{unused}
}
\value{
Output the mean difference, total standard deviation ratio or within-subject standard deviation ratio values (depending on the \sQuote{method}) with 
the calculated corresponding profile likelihood values.
}
\description{
Print \sQuote{poi} (mean difference, total standard deviation ratio or within-subject standard deviation ratio) and \sQuote{maxLik} (corresponding profile likelihood) elements from a proLikelihood object.
}
\examples{
\donttest{
data(dat)
cols <- list(subject = 'subject', formula = 'formula', y = 'AUC')
l <- averageBE(dat, colSpec = cols, xlength = 300)
l
}

}
