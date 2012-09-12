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
riley(X, type="effect.size", optimization = "Nelder-Mead", control = list(), ...)
\method{riley}{default}(X, type="effect.size", optimization = "Nelder-Mead", control = list(), ...)
}

\arguments{
  \item{X}{data frame containing integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN} (for diagnostic test accuracy data, cfr. \code{\link{rileyDA}}) or numeric variables \code{Y1}, \code{vars1}, \code{Y2} and \code{vars2} (for effect size data, cfr. \code{\link{rileyES}}).}
  \item{type}{a character string defining the type of data that is being summarized. Defaults to \code{effect.size} for summarizing effect sizes for which the normality assumption holds (for more details see \code{\link{rileyES}}). Diagnostic test accuracy data (i.e. sensitivities and specificities) can be pooled by choosing \code{test.accuracy} (for more details see \code{\link{rileyDA}}).}
  \item{optimization}{The optimization method that should be used for minimizing the negative (restricted) log-likelihood function. The default method is an implementation of that of Nelder and Mead (1965), that uses only function values and is robust but relatively slow. It will work reasonably well for non-differentiable functions. Other methods are described in \code{\link{optim}}.}
  \item{control}{ A list of control parameters to pass to \code{optim}.}
  \item{\dots}{arguments to be passed on to other functions.}
}

\details{ Parameters are estimated by iteratively maximizing the restriced log-likelihood using the Newton-Raphson procedure. Algorithms for dealing with missing data are currently not implemented, but Bayesian approaches will become available in later versions. }

\value{ An object of the class \code{riley} for which many standard methods are available.}

\references{
Nelder JA, Mead R. A simplex algorithm for function minimization. \emph{Computer Journal} (1965); \bold{7}: 308--313. 

Daniels MJ, Hughes MD. Meta-analysis for the evaluation of potential surrogate markers. \emph{Statistics in Medicine} 1997; \bold{16}: 1965--1982.

van Houwelingen HC, Arends LR, Stijnen T. Advanced methods in meta-analysis: multivariate approach and meta-regression. \emph{Statistics in Medicine} 2002; \bold{21}: 589--624.

Reitsma J, Glas A, Rutjes A, Scholten R, Bossuyt P, Zwinderman A. Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews. \emph{Journal of Clinical Epidemiology} 2005; \bold{58}: 982--990.

Korn EL, Albert PS, McShane LM. Assessing surrogates as trial endpoints using mixed models. \emph{Statistics in Medicine} 2005; \bold{24}: 163--182.

Thompson JR, Minelli C, Abrams KR, Tobin MD, Riley RD. Meta-analysis of genetic studies using mendelian randomization--a multivariate approach. \emph{Statistics in Medicine} 2005; \bold{24}: 2241--2254.

Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
}

\author{Thomas Debray <thomas.debray@gmail.com>}

\seealso{
\code{\link{logLik.riley}}, \code{\link{plot.riley}}, \code{\link{predict.riley}}, \code{\link{rileyDA}}, \code{\link{rileyES}}, \code{\link{summary.riley}}, \code{\link{vcov.riley}}
}

\examples{
data(Scheidler)
data(Daniels)
data(Kertai)

#Meta-analysis of potential surrogate markers data
fit1 <- riley(Daniels) #Maxit reached, try again with more iterations
fit1 <- riley(Daniels,control=list(maxit=10000))
summary(fit1)

#Meta-analysis of prognostic test studies
fit2 <- riley(Kertai,type="test.accuracy")
summary(fit2)

#Meta-analysis of computed tomography data 
ds <- Scheidler[which(Scheidler$modality==1),]
fit3 <- riley(ds,type="test.accuracy")
summary(fit3)
}

\keyword{regression}
\keyword{multivariate}
\keyword{bivariate}
