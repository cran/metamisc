\name{metamisc-package}
\alias{metamisc-package}
\docType{package}
\title{Diagnostic and prognostic meta-analysis}
\description{This package provides functions for diagnostic and prognostic meta-analyses. It estimates univariate, bivariate and multivariate models, and allows the aggregation of previously published prediction models with new data.}
\details{
\tabular{ll}{
Package: \tab metamisc\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2013-05-15\cr
License: \tab GPL-2 \cr
}
The package provides tools for the meta-analysis of individual participant (IPD) and/or aggregate data (AD). At this stage, it is possible to pool univariate (with \code{\link{uvmeta}}) and bivariate (with \code{\link{riley}}) summary data using frequentist and Bayesian approaches.}
\author{Thomas Debray <thomas.debray@gmail.com>}

\references{
DerSimonian R, Laird N. Meta-analysis in clinical trials. \emph{Controlled Clinical Trials} 1986; \bold{7}: 177--188.

Reitsma J, Glas A, Rutjes A, Scholten R, Bossuyt P, Zwinderman A. Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews. \emph{Journal of Clinical Epidemiology} 2005; \bold{58}: 982--990.

Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
}

\keyword{ package }
\keyword{regression}
\keyword{meta-analysis}

\seealso{ \code{\link{Daniels}}, \code{\link{Kertai}}, \code{\link{riley}}, \code{\link{Scheidler}}, \code{\link{uvmeta}}}
