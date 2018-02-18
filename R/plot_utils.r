#' Forest plot
#' 
#' Generate a forest plot by specifying the various effect sizes, confidence intervals and summary estimate.
#' @param theta Numeric vector with effect size for each study
#' @param theta.ci Two-dimensional array specifying the lower bound (first column) and upper bound (second column) of the 
#' confidence interval of the effect sizes
#' @param theta.slab Character vector specifying the study labels
#' @param theta.summary Meta-analysis summary estimate of the effect sizes
#' @param theta.summary.ci Numeric vector specifying the lower bound (first item) and upper bound (second item) of the 
#' confidence interval of the summary estimate
#' @param theta.summary.pi Numeric vector specifying the lower bound (first item) and upper bound (second item) of the 
#' prediction interval of the summary estimate. 
#' @param title Title of the forest plot
#' @param sort By default, studies are sorted by ascending effect size (\code{sort="asc"}). Set to \code{"desc"} for 
#' sorting in reverse order, or any other value to ignore sorting.
#' @param theme Theme to generate the forest plot. By default, the classic dark-on-light ggplot2 theme is used. 
#' See \link[ggplot2]{theme_bw} for more information.
#' @param predint.linetype The linetype of the prediction interval
#' @param xlim The \code{x} limits \code{(x1, x2)} of the forest plot
#' @param xlab Optional character string specifying the X label
#' @param refline Optional numeric specifying a reference line
#' @param label.summary Optional character string specifying the label for the summary estimate
#' @param \dots Additional arguments, which are currently ignored.
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
forest <- function (theta, theta.ci, theta.slab, theta.summary, 
                    theta.summary.ci, 
                    theta.summary.pi = c(NA, NA),
                    title,
                    sort = "asc",
                    theme = theme_bw(),
                    predint.linetype=1,
                    xlim,
                    xlab="", 
                    refline=0,
                    label.summary="Summary Estimate", ...) {

  if (missing(theta)) stop("Study effect sizes are missing!")
  if (missing(theta.ci)) stop("Confidence intervals of effect sizes missing!")
  if (missing(theta.slab)) stop("Study labels are missing!")
  
  num.studies <- unique(c(length(theta), dim(theta.ci)[1], length(theta.slab)))
  if (length(num.studies)>1) stop(paste("Mismatch in data dimensions!"))
  
  label.predint <- "Prediction Interval"
  
  #Extract data
  yi <- theta
  slab <- theta.slab
  
  #Sort data
  if (sort=="asc") {
    i.index <- order(yi)
  } else if (sort=="desc") {
    i.index <- order(yi, decreasing=TRUE)
  } else {
    i.index <- 1:length(yi)
  }
  
  
  # Add meta-analysis results
  if (NA %in% theta.summary.pi) {
    scat  <- c(rep(1, num.studies), 0) #indicator variable for data points
    slab  <- c(slab[i.index], label.summary)
    yi    <- c(yi[i.index], theta.summary)
    ci.lb <- c(theta.ci[i.index,1], theta.summary.ci[1])
    ci.ub <- c(theta.ci[i.index,2], theta.summary.ci[2])
  } else {
    scat  <- c(rep(1,length(i.index)), 0, 0)
    slab  <- c(slab[i.index], label.summary, label.predint)
    yi    <- c(yi[i.index], theta.summary, theta.summary)
    ci.lb <- c(theta.ci[i.index,1], theta.summary.ci[1], theta.summary.pi[1])
    ci.ub <- c(theta.ci[i.index,2], theta.summary.ci[2], theta.summary.pi[2])
  }
  
  
  ALL <- data.frame(study=slab, mean=yi, m.lower=ci.lb, m.upper=ci.ub, order=length(yi):1, scat=scat)
  
  # reorder factor levels based on another variable (HPD.mean)
  ALL$study.ES_order <- reorder(ALL$study, ALL$order, mean)    
  
  p <- with(ALL, ggplot(ALL[!is.na(ALL$mean), ], 
                        aes(x = study.ES_order, y = mean, ymin = m.lower, ymax = m.upper)) +
              geom_pointrange(data = subset(ALL, scat == 1)) + 
              scale_x_discrete(limits=rev(slab)) + #change order of studies
              coord_flip() + 
              theme +
              ylab(xlab) + 
              xlab(""))
  
  if (!missing(xlim)) {
    p <- p + ylim(xlim)
  }
  
  # Add title
  if (!missing(title)) {
    p <- p + ggtitle(title)
  }
  
  # Add refline
  if (!missing(refline)) {
    if (is.numeric(refline)) {
      p <- p + geom_hline(yintercept = refline,  linetype = "dotted") 
    }
  }
  
  # Add meta-analysis summary
  g2 <- with(ALL, subset(ALL, study == label.summary))
  g2$pi.upper <- theta.summary.pi[2]
  g2$pi.lower <- theta.summary.pi[1]
  g2$ci.upper <- theta.summary.ci[2]
  g2$ci.lower <- theta.summary.ci[1]
  
  # Add meta-analysis summary
  g3 <- with(ALL, subset(ALL, study == label.predint))
  g3$pi.upper <- theta.summary.pi[2]
  g3$pi.lower <- theta.summary.pi[1]
  
  # Prediction interval
  if (!is.na(g3$pi.lower) & !is.na(g3$pi.upper)) {
    p <- p +with(g3, geom_errorbar(data=g3, aes(ymin = pi.lower, ymax = pi.upper, x=label.predint), 
                                   width = 0.5, size=1.0, linetype=predint.linetype))
  }

  # Confidence interval
  p <- p + with(g2, geom_errorbar(data=g2, aes(ymin = ci.lower, ymax = ci.upper, x=label.summary), width = 0.5, size=1.0))

  # Summary estimate
  p <- p + with(g2, geom_point(data=g2, shape=23, size=3, fill="white"))

  p
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  #library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}