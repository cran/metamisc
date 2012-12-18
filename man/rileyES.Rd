\name{rileyES}
\alias{rileyES}

\concept{random effects model}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}

\title{ Fit the alternative model for bivariate random-effects meta-analysis (Riley) }

\description{ This function fits the alternative model for bivariate random-effects meta-analysis on effect size data when the within-study correlations are unknown. This bivariate model was proposed by Riley et al. (2008) and is similar to the bivariate random-effects model from Reitsma et al. (2005), but includes an overall correlation parameter rather than separating the (usually unknown) within- and between-study correlation. As a consequence, the alternative model is not fully hierarchical, and estimates of additional variation beyond sampling error (\code{psi}) are not directly equivalent to the between-study variation (\code{tau}) from the general model. Furthermore, it has been argumented that assuming zero within-study correlations (i.e. applying Reitsma's approach) is reasonable when summarizing the sensitivities and false positive rates of a diagnostic test (Reitsma et al. 2005, Daniels and Hughes 1997, Korn et al. 2005, Thompson et al. 2005, Van Houwelingen et al. 2002). The alternative model for bivariate random-effects meta-analysis may, however, be useful when there is large within-study variability, few primary studies are available or the general model estimates the between-study correlation as 1 or -1.}

\usage{
rileyES(X = NULL, Y1, Y2, vars1, vars2, optimization = "Nelder-Mead", control = list(),...)
}

\arguments{
  \item{X}{any object that can be converted to a data frame with integer variables \code{Y1}, \code{vars1}, \code{Y2} and \code{vars2}.}
  \item{Y1}{vector of numerics representing the effect sizes of outcome 1, ignored if \code{X} is not \code{NULL}}
  \item{vars1}{vector of numerics representing the error variances of \code{Y1}, ignored if \code{X} is not \code{NULL}}
  \item{Y2}{vector of numerics representing the effect sizes of outcome 2, ignored if \code{X} is not \code{NULL}}
  \item{vars2}{vector of numerics representing the error variances of \code{Y2}, ignored if \code{X} is not \code{NULL}}
  \item{optimization}{The optimization method that should be used for minimizing the negative (restricted) log-likelihood function. The default method is an implementation of that of Nelder and Mead (1965), that uses only function values and is robust but relatively slow. Other methods are described in \code{\link{optim}}.}
  \item{control}{ A list of control parameters to pass to \code{\link{optim}}.}
  \item{\dots}{arguments to be passed on to other functions, currently ignored}
}

\details{ The following parameters are estimated by iteratively maximizing the restriced log-likelihood using the Newton-Raphson procedure: logit of sensitivity (\code{beta1}), logit of false positive rate (\code{beta2}), additional variation of \code{beta1} beyond sampling error (\code{psi1}), additional variation of \code{beta2} beyond sampling error (\code{psi2}) and a transformation of the correlation between \code{psi1} and \code{psi2} (\code{rhoT}). The original correlation is given as \code{inv.logit(rhoT)*2-1}. The results from a univariate random-effects meta-analysis with a method-of-moments estimator are used as starting values for \code{beta1}, \code{beta2}, \code{psi1} and \code{psi2} in the \code{optim} command. The starting value for \code{rhoT} is 0. Standard errors for all parameters are obtained from the inverse Hessian matrix. }

\value{ An object of the class \code{riley} for which many standard methods are available. A warning message is casted when the Hessian matrix contains negative eigenvalues, which implies that the identified solution is a saddle point and thus not optimal.}

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


\keyword{regression}
\keyword{multivariate}
\keyword{bivariate}
