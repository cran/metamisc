\name{plot.valmeta}
\alias{plot.valmeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forest Plots
}
\description{
Function to create forest plots for objects of class \code{"valmeta"}.
}
\usage{
\method{plot}{valmeta}(x, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{an object of class \code{"valmeta"}}
    \item{\ldots}{Additional arguments which are passed to \code{forest} from the package \code{metafor}.}
}
\details{
Plots are generated using functionalities provided by the \code{metafor} pacakge. The forest plot shows the performance estimates of each validation with corresponding confidence intervals. A polygon is added to the bottom of the forest plot, showing the summary estimate based on the model (with the outer edges of the polygon indicating the confidence interval limits). A 95\% prediction interval is added by default,  the dotted line indicates its (approximate) bounds.
}

\references{
Debray TPA, Damen JAAG, Snell KIE, Ensor J, Hooft L, Reitsma JB, et al. A guide to systematic review and meta-analysis of prediction model performance. \emph{BMJ}. 2017;356:i6460. 

Lewis S, Clarke M. Forest plots: trying to see the wood and the trees. \emph{BMJ}. 2001; 322(7300):1479--80. 

Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. \emph{BMJ}. 2011 342:d549--d549. 

Viechtbauer W. Conducting Meta-Analyses in R with the metafor Package. \emph{Journal of Statistical Software}. 2010; 36(3). Available from: \url{http://www.jstatsoft.org/v36/i03/}

}
\author{Thomas Debray <thomas.debray@gmail.com>}
\note{
As indicated by \code{metafor}, the labels, annotations, and symbols may become quite small and impossible to read when the number of studies is quite large. Stretching the plot window vertically may then provide a more readable figure (one should call the function again after adjusting the window size, so that the label/symbol sizes can be properly adjusted). Also, the \code{cex}, \code{cex.lab}, and \code{cex.axis} arguments are then useful to adjust the symbol and text sizes.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(EuroSCORE)
attach(EuroSCORE)

fit <- valmeta(cstat=c.index, cstat.se=se.c.index, 
               cstat.95CI=cbind(c.index.95CIl,c.index.95CIu),
               N=n, O=n.events)
plot(fit)

detach(EuroSCORE)
}
\keyword{ meta-analysis }
\keyword{ discrimination }
\keyword{ calibration }
