\name{metamisc-package}
\alias{metamisc-package}
\docType{package}
\title{Diagnostic and Prognostic Meta-Analysis}
\description{Meta-analysis of diagnostic and prognostic modeling studies. Summarize estimates of diagnostic test accuracy and prediction model performance. Validate, update and combine published prediction models.}

\details{
%\tabular{ll}{
%Package: \tab metamisc\cr
%Type: \tab Package\cr
%Version: \tab 0.1.7\cr
%Date: \tab 2017-11-19\cr
%License: \tab GPL-2 \cr
%}
The following functionality is currently implemented: univariate meta-analysis of summary data (\code{\link{uvmeta}}), bivariate meta-analysis of correlated outcomes (\code{\link{riley}}), meta-analysis of prediction model performance (\code{\link{valmeta}}), evaluation of funnel plot asymmetry (\code{\link{fat}}). \cr
\cr
The \code{metamisc} package also provides a comprehensive framework for developing prediction models when patient-level data from multiple studies or settings are available (\code{\link{metapred}}).
}

\author{Thomas Debray <thomas.debray@gmail.com>, Valentijn de Jong}

\references{
Debray TPA, Moons KGM, Ahmed I, Koffijberg H, Riley RD. A framework for developing, implementing, and evaluating clinical prediction models in an individual participant data meta-analysis. \emph{Stat Med}. 2013;32(18):3158--80. 

Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review and meta-analysis of prediction model performance. \emph{BMJ}. 2017; \bold{356}:i6460. 

Debray TPA, Moons KGM, Riley RD. Detecting small-study effects and funnel plot asymmetry in meta-analysis of survival data: a comparison of new and existing tests. \emph{Res Syn Meth}. 2017.

Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
}

\seealso{
\code{\link{fat}}, \code{\link{metapred}}, \code{\link{riley}}, \code{\link{uvmeta}}, \code{\link{valmeta}}  
}

\keyword{package}
\keyword{regression}
\keyword{meta-analysis}
\keyword{validation}

