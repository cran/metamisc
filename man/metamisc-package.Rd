\name{riley-package}
\alias{riley-package}
\docType{package}
\title{Diagnostic and prognostic meta-analysis}
\description{This package provides functions for diagnostic and prognostic meta-analyses. It estimates univariate, bivariate and multivariate models, and allows the aggregation of previously published prediction models with new data.}
\details{
\tabular{ll}{
Package: \tab metamisc\cr
Type: \tab Package\cr
Version: \tab 0.0.2\cr
Date: \tab 2012-08-30\cr
License: \tab GPL-2 \cr
}
The package provides tools for the meta-analysis of individual participant (IPD) and/or aggregate data (AD). Currently, it is possible to pool univariate summary (with \code{\link{uvmeta}}) and diagnostic accuracy (with \code{\link{riley}}) data. Whereas the former applies a univariate meta-analysis using DerSimonian and Laird's method (method-of-moment estimator), the latter implements a bivariate meta-analysis (Restricted Maximum Likelihood) using the alternative model for bivariate random-effects meta-analysis by Riley et al. (2008). For this the number true positives (TP), false negatives (FN), true negatives (TN) and false postives (FP) for each study must be known.}
\author{Thomas Debray <thomas.debray@gmail.com>}
\references{
DerSimonian, R., & Laird, N. (1986). 
\dQuote{Meta-analysis in clinical trials.} \emph{Controlled Clinical Trials}, \bold{7}, 177--188.

Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.

Riley, RD., Thompson, JR., & Abrams, KR. (2008).
\dQuote{An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown.} \emph{Biostatistics}, \bold{9}, 172--186.
}
\keyword{ package }
\keyword{regression}
\keyword{meta-analysis}
\seealso{
\code{\link{riley}}
\code{\link{uvmeta}}
}
