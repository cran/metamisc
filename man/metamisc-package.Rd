\name{metamisc-package}
\alias{metamisc-package}
\docType{package}
\title{Meta-Analysis of Diagnosis and Prognosis Research Studies}
\description{Facilitate meta-analysis of diagnosis and prognosis research studies. It includes functions to  summarize multiple estimates of prediction model discrimination and calibration performance, as described by Debray et al. (2019) <doi:10.1177/0962280218785504>. It also includes functions to evaluate funnel plot asymmetry, as described by Debray et al. (2018) <doi:10.1002/jrsm.1266>. Finally, the package provides functions for developing multivariable prediction models from datasets with clustering.}

\details{
%\tabular{ll}{
%Package: \tab metamisc\cr
%Type: \tab Package\cr
%Version: \tab 0.2.3\cr
%Date: \tab 2020-09-23\cr
%License: \tab GPL-2 \cr
%}
The following functionality is currently implemented: univariate meta-analysis of summary data (\code{\link{uvmeta}}), bivariate meta-analysis of correlated outcomes (\code{\link{riley}}), meta-analysis of prediction model performance (\code{\link{valmeta}}), evaluation of funnel plot asymmetry (\code{\link{fat}}). \cr
\cr
The \code{metamisc} package also provides a comprehensive framework for developing prediction models when patient-level data from multiple studies or settings are available (\code{\link{metapred}}).
}

\author{Thomas Debray <thomas.debray@gmail.com>, Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>}

\references{
Debray TPA, Moons KGM, Ahmed I, Koffijberg H, Riley RD. A framework for developing, implementing, and evaluating clinical prediction models in an individual participant data meta-analysis. Stat Med. 2013;32(18):3158--80. 

Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review and meta-analysis of prediction model performance. BMJ. 2017;356:i6460. 

Debray TPA, Moons KGM, Riley RD. Detecting small-study effects and funnel plot asymmetry in meta-analysis of survival data: a comparison of new and existing tests. Res Syn Meth. 2018;9(1):41--50. 

Debray TPA, Damen JAAG, Riley R, Snell KIE, Reitsma JB, Hooft L, et al. A framework for meta-analysis of prediction model studies with binary and time-to-event outcomes. Stat Methods Med Res. 2019 Sep;28(9):2768--86. 

Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. Biostatistics 2008; 9: 172--186.

Riley RD, Moons K, Snell KIE, Ensor J, Hooft L, Altman D, et al. A guide to systematic review and meta-analysis of prognostic factor studies. BMJ. 2019;364:k4597. 

Steyerberg EW, Nieboer D, Debray TPA, Van Houwelingen JC. Assessment of heterogeneity in an individual participant data meta-analysis of prediction models: An overview and illustration. Stat Med. 2019;38(22):4290--309. 
}

\seealso{
\code{\link{fat}}, \code{\link{metapred}}, \code{\link{riley}}, \code{\link{uvmeta}}, \code{\link{valmeta}}  
}

\keyword{package}
\keyword{regression}
\keyword{meta-analysis}
\keyword{validation}

