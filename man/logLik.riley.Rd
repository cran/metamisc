\name{logLik.riley}
\alias{logLik.riley}

\title{Print the log-likelihood}

\description{This function provides the (restricted) log-likelihood of a fitted model.}

\usage{\method{logLik}{riley}(object, \ldots)}

\arguments{
  \item{object}{a \code{riley} object, representing a fitted alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown.}
  \item{\dots}{arguments to be passed on to other functions, currently ignored}
}

\value{ Returns an object of class \code{logLik}. This is the (restricted) log-likelihood of the model represented by \code{object} evaluated at the estimated coefficients. It contains at least one attribute, "\code{df}" (degrees of freedom), giving the number of (estimated) parameters in the model. }

\references{
Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
}

\author{Thomas Debray <thomas.debray@gmail.com>}

\examples{
data(Scheidler)
ds <- Scheidler[which(Scheidler$modality==1),]
fit <- riley(ds, type="test.accuracy")
logLik(fit)
}

\seealso{
\code{\link{plot.riley}}, \code{\link{predict.riley}}, \code{\link{summary.riley}}, \code{\link{riley}}, \code{\link{rileyDA}}, \code{\link{rileyES}}
}

\keyword{ likelihood }
\keyword{ riley }
\keyword{ bivariate }
\keyword{ meta-analysis }
