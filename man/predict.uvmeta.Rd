\name{predict.uvmeta}
\alias{predict.uvmeta}

\title{ Prediction Interval }

\description{ Calculates a prediction interval for the effect size. This interval predicts in what range future effect sizes will fall given what has already been observed. }

\usage{ \method{predict}{uvmeta}(object, level = 0.95, \ldots) }

\arguments{
  \item{object}{ an \code{uvmeta} object.}
  \item{level}{ numeric, the level for calculations of confidence intervals }
  \item{\dots}{ arguments to be passed on to other functions }
}

\details{ Prediction intervals are based on Student's t-distribution with (numstudies - 2) degrees of freedom. }

\value{ Vector containing the prediction interval of the effect size. }

\author{ Thomas Debray <thomas.debray@gmail.com> }

\references{
Riley, RD., Higgins, JPT., & Deeks, JJ. (2011). 
\dQuote{Interpretation of random effects meta-analyses.} \emph{British Medical Journal}, \bold{342}, d549.
}

\seealso{
   \code{\link{uvmeta}}
}


