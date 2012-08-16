\name{mada-data}
\docType{data}
\alias{mada-data}
\alias{Scheidler}

\title{Diagnostic accuracy data}

\description{One data frame with diagnostic accuracy data from binary test outcomes.}

\usage{
data("Scheidler")
}


\format{One data frames with 6 variables.
  \describe{
  \item{author}{string . author of article}
  \item{modality}{integer . type of test (1=CT, 2=LAG, 3=MRI)}
  \item{TP}{integer. number of true positives}  
  \item{FN}{integer. number of false negatives}  
  \item{FP}{integer. number of false positives}  
  \item{TN}{integer. number of true negatives}  
  }
}


\details{
The \code{Scheidler} data comprises the results from a meta-analysis where three imaging techniques for the diagnosis of lymph node metastasis in women with cervical cancer are compared. Forty-four studies in total were included: 17 studies evaluated lymphangiography, another 17 studies examined computed tomography and the remaining 10 studies focused on magnetic resonance imaging. Diagnosis of metastatic disease by lymphangiography (LAG) is based on the presence of nodal-filling defects, whereas computed tomography (CT) and magnetic resonance imaging (MRI) rely on nodal enlargement.
}

\source{
Scheidler, J., Hricak, H., Yu, K.K., Subak, L., & Segal M.R. (1997).
\dQuote{Radiological evaluation of lymph node metastases in patients with cervical cancer. A meta-analysis.} \emph{Journal of the American Medical Association}, \bold{278}, 1096--1101.

Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.
}

\keyword{datasets}
