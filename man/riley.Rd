\name{riley}
\alias{riley}
\alias{riley.default}
\alias{print.riley}

\concept{random effects model}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}

\title{ Fit the alternative model for bivariate random-effects meta-analysis (Riley) }

\description{ This function fits the alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. This bivariate model was proposed by Riley et al. (2008) and is similar to the bivariate random-effects model from Reitsma et al. (2005), but includes an overall correlation parameter rather than separating the (usually unknown) within- and between-study correlation. As a consequence, the alternative model is not fully hierarchical, and estimates of additional variation beyond sampling error (\code{psi}) are not directly equivalent to the between-study variation (\code{tau}) from the general model. Furthermore, it has been argumented that assuming zero within-study correlations (i.e. applying Reitsma's approach) is reasonable when summarizing the sensitivities and false positive rates of a diagnostic test (Reitsma et al. 2005, Daniels and Hughes 1997, Korn et al. 2005, Thompson et al. 2005, Van Houwelingen et al. 2002). The alternative model for bivariate random-effects meta-analysis may, however, be useful when there is large within-study variability, few primary studies are available or the general model estimates the between-study correlation as 1 or -1.}

\usage{
riley(X, ...)
\method{riley}{default}(X = NULL, TP, FN, FP, TN, correction = 0.5, 
    correction.control = "all", optimization = "Nelder-Mead", control = list(), ...)
}

\arguments{
  \item{X}{any object that can be converted to a data frame with integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN}, alternatively a matrix with column names including \code{TP}, \code{FN}, \code{FP} and \code{TN}.}
  \item{TP}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{FN}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{FP}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{TN}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{correction}{numeric, continuity correction applied if zero cells}
  \item{correction.control}{character, if set to \code{"all"} (the default) the continuity correction is added to the whole data if only one cell in one study is zero. If set to \code{"single"} the correction is only applied to rows of the data which have a zero.}
  \item{optimization}{The optimization method that should be used for minimizing the negative (restricted) log-likelihood function. The default method is an implementation of that of Nelder and Mead (1965), that uses only function values and is robust but relatively slow. It will work reasonably well for non-differentiable functions. Other methods are described in \code{\link{optim}}.}
  \item{control}{ A list of control parameters to pass to \code{optim}.}
  \item{\dots}{arguments to be passed on to other functions, currently ignored}
}

\details{ The following parameters are estimated by iteratively maximizing the restriced log-likelihood using the Newton-Raphson procedure: logit of sensitivity (\code{mu1}), logit of false positive rate (\code{mu2}), additional variation of \code{mu1} beyond sampling error (\code{psi1}), additional variation of \code{mu2} beyond sampling error (\code{psi2}) and a transformation of the correlation between \code{psi1} and \code{psi2} (\code{rhoT}). The original correlation is given as \code{inv.logit(rhoT)*2-1}. The results from a univariate random-effects meta-analysis with a method-of-moments estimator are used as starting values for \code{mu1}, \code{mu2}, \code{psi1} and \code{psi2} in the \code{optim} command. The starting value for \code{rhoT} is 0. Standard errors for all parameters are obtained from the inverse Hessian matrix. }

\value{ An object of the class \code{riley} for which many standard methods are available. A warning message is casted when the Hessian matrix contains negative eigenvalues, which implies that the identified solution is not optimal.}

\references{
Nelder, JA., & Mead, R. (1965) 
\dQuote{A simplex algorithm for function minimization.} 
\emph{Computer Journal}, \bold{7}, 308--313. 

Daniels, MJ., & Hughes, MD. (1997). 
\dQuote{Meta-analysis for the evaluation of potential surrogate markers.} 
\emph{Statistics in Medicine}, \bold{16}, 1965--1982.

van Houwelingen, HC., Arends, LR., & Stijnen, T. (2002). 
\dQuote{Advanced methods in meta-analysis: multivariate approach and meta-regression.} 
\emph{Statistics in Medicine}, \bold{21}, 589--624.

Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.

Korn, EL., Albert, PS., & McShane, LM. (2005). 
\dQuote{Assessing surrogates as trial endpoints using mixed models.} 
\emph{Statistics in Medicine}, \bold{24}, 163--182.

Thompson, JR., Minelli, C., Abrams, KR., Tobin, MD., & Riley, RD. (2005). 
\dQuote{Meta-analysis of genetic studies using mendelian randomization--a multivariate approach.} 
\emph{Statistics in Medicine}, \bold{24}, 2241-2254.

Riley, RD., Thompson, JR., & Abrams, KR. (2008).
\dQuote{An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown.} \emph{Biostatistics}, \bold{9}, 172--186.

}
\author{
Thomas Debray <thomas.debray@gmail.com>
}

\seealso{
\code{\link{logLik.riley}} \code{\link{plot.riley}} \code{\link{predict.riley}} \code{\link{summary.riley}} \code{\link{vcov.riley}} 
}

\examples{
data(Scheidler)

#Select computed tomography data
ds <- Scheidler[which(Scheidler$modality==1),]

#Perform the analysis
fit <- riley(ds,correction.control="single")
summary(fit)
plot(fit)
}

\keyword{regression}
\keyword{multivariate}
\keyword{bivariate}
