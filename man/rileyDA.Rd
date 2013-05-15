\name{rileyDA}
\alias{rileyDA}

\concept{random effects model}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}

\title{ Fit the alternative model for bivariate random-effects meta-analysis (Riley) }

\description{ This function fits the alternative model for bivariate random-effects meta-analysis on diagnostic test accuracy data when the within-study correlations are unknown assumed to be different from zero. A transformation is applied to the sensitivities ans false positive rates of each study, in order to meet the normality assumptions of the model.}

\usage{
rileyDA(X = NULL, TP, FN, FP, TN, correction = 0.5, 
           correction.control = "all", optimization = "Nelder-Mead", 
           control = list(), ...)
}

\arguments{
  \item{X}{any object that can be converted to a data frame with integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN}.}
  \item{TP}{vector of integers representing the number of true positives, ignored if \code{X} is not \code{NULL}}
  \item{FN}{vector of integers representing the number of false negatives, ignored if \code{X} is not \code{NULL}}
  \item{FP}{vector of integers representing the number of false positives, ignored if \code{X} is not \code{NULL}}
  \item{TN}{vector of integers representing the number of true negatives, ignored if \code{X} is not \code{NULL}}
  \item{correction}{numeric, continuity correction applied if zero cells}
  \item{correction.control}{character, if set to \code{"all"} (the default) the continuity correction is added to the whole data if only one cell in one study is zero. If set to \code{"single"} the correction is only applied to rows of the data which have a zero.}
  \item{optimization}{The optimization method that should be used for minimizing the negative (restricted) log-likelihood function. The default method is an implementation of that of Nelder and Mead (1965), that uses only function values and is robust but relatively slow. Other methods are described in \code{\link{optim}}.}
  \item{control}{ A list of control parameters to pass to \code{\link{optim}}.}
  \item{\dots}{arguments to be passed on to other functions, currently ignored}
}

\details{ The following parameters are estimated using \code{\link{rileyES}}: logit of sensitivity (\code{beta1}), logit of false positive rate (\code{beta2}), additional variation of \code{beta1} beyond sampling error (\code{psi1}), additional variation of \code{beta2} beyond sampling error (\code{psi2}) and a transformation of the correlation between \code{psi1} and \code{psi2} (\code{rhoT}). The original correlation is given as \code{inv.logit(rhoT)*2-1}. The results from a univariate random-effects meta-analysis with a method-of-moments estimator are used as starting values for \code{beta1}, \code{beta2}, \code{psi1} and \code{psi2} in the \code{optim} command. The starting value for \code{rhoT} is 0. Standard errors for all parameters are obtained from the inverse Hessian matrix. }

\value{ An object of the class \code{riley} for which many standard methods are available. A warning message is casted when the Hessian matrix contains negative eigenvalues, which implies that the identified solution is a saddle point and thus not optimal.}

\author{Thomas Debray <thomas.debray@gmail.com>}


\keyword{regression}
\keyword{multivariate}
\keyword{bivariate}
