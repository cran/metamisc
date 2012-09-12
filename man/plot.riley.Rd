\name{plot.riley}
\alias{plot.riley}

\title{
Plot the summary of the bivariate model from Riley et al. (2008).
}
\description{
This function plots the summary sensitivity and false positive rate with their corresponding confidence regions.
}
\usage{\method{plot}{riley}(x, plotsumm = TRUE, plotnumerics = TRUE, level = 0.95, main="",
                       ylim = c(0,1), xlim = c(0,1), pch = 1, lty = 1, lwd = 1, cex.numerics=0.45,
                       add=FALSE, \ldots)}

\arguments{
  \item{x}{a \code{riley} object.}
  \item{plotsumm}{logical, should the plot draw the summary pair of sensitivity and false positive rate?}
  \item{plotnumerics}{logical, should the plot contain a summary table of sensitivity and false positive rate?}
  \item{level}{numeric, the level for calculations of confidence intervals}
  \item{main}{string, title of the plot}
  \item{ylim}{numeric of length 2, which section of the sensitivities to plot?}
  \item{xlim}{numeric of length 2, which section of the false positive rates to plot?}
  \item{pch}{integer, symbol for the pair of mean sensitivity and false positive rate}
  \item{lty}{integer, line type of confidence curve}
  \item{lwd}{integer, line width of the confidence curve}
  \item{cex.numerics}{numeric, text size}
  \item{add}{logical, should the confidence region be added to the current plot?}
  \item{\dots}{arguments to be passed on to other functions}
}

\references{
Riley RD, Thompson JR, Abrams KR. An alternative model for bivariate random-effects meta-analysis when the within-study correlations are unknown. \emph{Biostatistics} 2008; \bold{9}: 172--186.
}

\author{Thomas Debray <thomas.debray@gmail.com>}

\seealso{
\code{\link{riley}}
}

\examples{
data(Scheidler)

ds1 <- Scheidler[which(Scheidler$modality==1),]
ds2 <- Scheidler[which(Scheidler$modality==2),]
ds3 <- Scheidler[which(Scheidler$modality==3),]

#Perform the analyses
fit1 <- riley(ds1, type="test.accuracy")
fit2 <- riley(ds2, type="test.accuracy")
fit3 <- riley(ds3, type="test.accuracy")

plot(fit1,plotnumerics=FALSE,pch=0) #CT
plot(fit2,plotnumerics=FALSE,add=TRUE,pch=1) #LAG
plot(fit3,plotnumerics=FALSE,add=TRUE,pch=2) #MRI
}

\keyword{riley}
\keyword{bivariate}
