\name{EuroSCORE}
\alias{EuroSCORE}
\docType{data}
\title{
Predictive performance of EuroSCORE II
}
\description{
This data set contains estimates on the predictive performance of the European system for cardiac operative risk evaluation (EuroSCORE II) in patients undergoing cardiac surgery. Results are based on the original development study and 22 validations identified by Guida \emph{et al}.
}
\usage{data("EuroSCORE")}
\format{
  A data frame with 23 observations on the following 13 variables.
  \describe{
    \item{\code{Study}}{The first author of the study}
    \item{\code{n}}{The total number of patients on which performance estimates are based}
    \item{\code{n.events}}{The total number of observed events}
    \item{\code{c.index}}{The estimated concordance statistic}
    \item{\code{se.c.index}}{The standard error of the concordance statistic}
    \item{\code{c.index.95CIl}}{The lower bound of the 95\% confidence interval of the estimated concordance statistic}
    \item{\code{c.index.95CIu}}{The upper bound of the 95\% confidence interval of the estimated concordance statistic}
    \item{\code{Po}}{The overall observed event probability (on a scale from 0 to 100)}
    \item{\code{Pe}}{The overall expected event probability (on a scale from 0 to 100)}
    \item{\code{SD.Pe}}{The standard error of \code{Pe}}
    \item{\code{e.events}}{The total number of expected events}
    \item{\code{multicentre}}{An indicator of whether the study was a multicentre study}
    \item{\code{mean.age}}{The mean age of the patients}
    \item{\code{sd.age}}{The spread of the age of the patients}
    \item{\code{pts.before.2010}}{An indicator of whether the study included patients before 2010 (i.e., before EuroSCORE II was developed)}
  }
}
\details{
Published in 2012, EuroSCORE II was developed using logistic regression in a dataset comprising 16,828 adult patients undergoing major cardiac surgery from 154 hospitals in 43 countries over a 12-weewk period (May-July) in 2010.  EuroSCORE II was developed to predict in-hospital mortality for patients undergoing any type of cardiac surgery. In 2014, a systematic review of published evidence on the performance value of the euroSCORE II was undertaken by Guida \emph{et al}. Twenty-two validations, including more 145,592 patients from 21 external validation articles (one study included two validations) and a split-sample validation contained within original development article were included in the review; 23 validation studies in total.
}
\source{
Guida P, Mastro F, Scrascia G, Whitlock R, Paparella D. Performance of the European System for Cardiac Operative Risk Evaluation II: a meta-analysis of 22 studies involving 145,592 cardiac surgery procedures. \emph{J Thorac Cardiovasc Surg}. 2014; \bold{148}(6):3049--3057.e1. 

Nashef SAM, Roques F, Sharples LD, Nilsson J, Smith C, Goldstone AR, et al. EuroSCORE II. \emph{Eur J Cardiothorac Surg}. 2012; \bold{41}(4):734-744; discussion 744-745. 

}

\examples{
data(EuroSCORE)
}
\keyword{datasets}
\keyword{EuroSCORE}