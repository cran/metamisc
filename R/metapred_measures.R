# TODO:
# add transforms to fema

##############################                Performance / error functions                  ###############################
# ### By convention, all performance measures:
# # Arguments:
# p     numeric vector of predicted probabilities
# y     numeric/integer vector of observed outcome, of same length as p.
# ...   for future compatibility
# 
# # Return, either:
# numeric of length 1. (old version)
# Or:
# object of class "mp.perf", where the first element, [[1]], is the estimate of performance
# Note that performance functions that do not produce normally distributed statistics,
# such as the auc, should not have "mp.perf" as its first class! The first class is used as a
# check of necessity of conversion to another scale.

##############################                

# Error function: Mean Squared Error # Replaced! See below!
# mse <- brier <- function(p, y, ...) mean((p - y)^2)

# rmse <- function(p, y, ...)
  # sqrt(mse(p = p, y = y, ...))

# Error function: Variance of prediction error
# var.e <- function(p, y, ...) 
#   var(p - y, ...)

# library(moments)
var.e <- var.e.with.se <- function(p, y, ...) 
  var.with.se(p - y)

# Necessary for var of var estimation.
# See https://math.stackexchange.com/questions/72975/variance-of-sample-variance
# and https://en.wikipedia.org/wiki/Variance , section Distribution of the sample variance
# and moments:kurtosis
sigma4 <- function(x, ...) 
  mean((x - mean(x, ...))^2, ...)^2

# Asymptotically unbiased estimate of variance of sample variance.
# https://math.stackexchange.com/questions/72975/variance-of-sample-variance
# x vector
# Returns variance, and se and variance of variance
var.with.se <- function(x, ...) {
  est <- var(x)
  v <-  2 * sigma4(x) / (length(x) - 1)
  out <- data.frame(estimate = est, se = sqrt(v), variances = v, n = length(x))
  class(out) <- c("mp.perf", class(out))
  out
}
# Measure 1: Coefficient of variation of prediction error.
  # abs logical absolute value
coef.var.pred <- function(p, y, abs = TRUE, ...)
  coef.var(x = p - y, abs = abs) 

#' @importFrom pROC auc
auc <- AUC <- AUROC <-  function(p, y, ...) {
  if (is.matrix(p))
    p <- p[, 1]
  if (is.matrix(y))
    y <- y[, 1]
  pROC::auc(response = y, predictor = p)
}

calibration.intercept <- cal.int <- function(p, y, estFUN, family, ...)
  pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "intercept")

bin.cal.int <- function(p, y, ...)
  pred.recal(p = p, y = y, estFUN = "glm", family = binomial, which = "intercept")

# Slope.only is a trick to make this functin work for metapred.
# Slope.only should otherwise always be false! Also: this messes up the variances,
# making meta-analysis impossible!
# multiplicative slope!
calibration.slope <- cal.slope <- function(p, y, estFUN, family, slope.only = TRUE, ...) {
  # refit <- pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "slope")
  # if (slope.only) {
  #   refit[[1]] <- refit[[1]][[2]]
  # }
  # refit
  
  refit <- pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "slope")
  if (slope.only) {
    refit$estimate <- refit[[1]] <- refit[[1]][2]
    refit$variances <- variances(refit)[2]
    class(refit) <- "mp.perf" # This should make it call the right confint method.
  }
  refit
}

# additive slope!
calibration.add.slope <- cal.add.slope <- function(p, y, estFUN, family, slope.only = TRUE, ...)  {
  
  refit <- pred.recal(p = p, y = y, estFUN = estFUN, family = family, which = "add.slope")
  if (slope.only) {
    refit$estimate <- refit[[1]] <- refit[[1]][2]
    refit$variances <- variances(refit)[2]
    class(refit) <- "mp.perf" # This should make it call the right confint method.
  }
  refit
}


# se.method character name of method for se calculation
# bs.n integer number of bootstrap samples
# ... Compatiblility only
# For asymptotic, see https://journals.ametsoc.org/doi/full/10.1175/2007WAF2007049.1
mse <- brier <- mse.with.se <- function(p, y, se.method = "asymptotic", bs.n = 10000, ...) {
  er <- p - y
  est <- mean(er^2)
  out <- data.frame(estimate = est, se = NA, variances = NA)
  out$n <- n <- length(p)
  
  if (se.method == "bootstrap") {
    ses <- rep(NA, bs.n)
    for (i in seq_len(bs.n))
      ses[i] <- mean(sample(er, length(er), replace = T)^2)
    out$se <- sd(ses)
    out$variances <- out$se^2
  } else if (se.method == "asymptotic") {
    out$variances <- var(er^2)/n
    out$se <- sqrt(out$variances)
  }
  
  class(out) <- c("mp.perf", class(out))
  out
}

#' @export
coef.mp.perf <- function(object, ...) 
  object$estimate

# Object mp.perf object
# use.fallback ignored
# compatibility only
#' @export
nobs.mp.perf <- function(object, use.fallback = FALSE, ...) 
  object$n

# Object mp.perf object
# parm "estimate" is the only viable obtion
# level confidence level
#' @export
confint.mp.perf <- function(object, parm = "estimate", level = .95, ...) { 
  ses <- se(object, ...)
  est <- object[[parm]]
  if(level < 0 || level > 1)
    stop("Impossible confidence level. Possible levels: 0 < level < 1")
  z <- qt(1 - (1 - level)/2, df = object$n - 1) # z = t distributed
  data.frame("ci.lb" = est - z * ses,"ci.ub" = est + z * ses)
}

# Object pROC::auc object
# parm "estimate" is the only viable obtion
# level ignored
# ... passed on to pROC::ci
#' @export
confint.auc <- function(object, parm = "estimate", level = .95, method = "delong", ...) {
  if(level < 0 || level > 1)
    stop("Impossible confidence level. Possible levels: 0 < level < 1")
  bounds <- pROC::ci(object, method = method)
  data.frame("ci.lb" = bounds[1],"ci.ub" = bounds[3]) # 2 is the point estimate.
}

# A generic function that calls the generic confint method
# Note that it does not have the parm parameter, because this never changes in metapred.
# Returns the confidence intervals with the specified names: ci.lb and ci.ub
get_confint <- function(object, level = 0.95, ...) {
  ci.bounds <- confint(object = object, level = level, ...)
  out <- data.frame(ci.lb = NA, ci.ub = NA)
  
  if (any(names(ci.bounds) == "ci.lb")) # as in metamisc and metafor (?)
    out$ci.lb <- ci.bounds[["ci.lb"]]
  else if (any(names(ci.bounds) == "2.5 %")) # as in glm
    out$ci.lb <- ci.bounds[["2.5 %"]]
  else if (any(names(ci.bounds) == "lower 0.95")) # as in logistf
    out$ci.lb <- ci.bounds[["lower 0.95"]]
  else if (any(names(as.data.frame(ci.bounds)) == "Lower 95%"))  # as in confint(logistf(...))
    out$ci.lb <- as.data.frame(ci.bounds)[["Lower 95%"]]
  
  if (any(names(ci.bounds) == "ci.ub")) # as in metamisc and metafor (?)
    out$ci.ub <- ci.bounds[["ci.ub"]]
  else if (any(names(ci.bounds) == "97.5 %")) # as in glm
    out$ci.ub <- ci.bounds[["97.5 %"]]
  else if (any(names(ci.bounds) == "upper 0.95")) # as in logistf
    out$ci.ub <- ci.bounds[["upper 0.95"]]
  else if (any(names(as.data.frame(ci.bounds)) == "Upper 95%"))  # as in confint(logistf(...))
    out$ci.ub <- as.data.frame(ci.bounds)[["Upper 95%"]]
  out
}

# m <- metamisc:::auc(runif(50, 0, 1), rbinom(50, 1, .5))
# aucci <- pROC:::ci(m)
# aucci <- pROC:::ci(m, method = "bootstrap")
# metamisc:::confint.auc(m)

# Binomial Log likelihood
# p Numeric vector, 0 <= p <= 1, predicted probability under the model
# y observed outcome, 0 or FALSE is no outcome, >= 1 or TRUE is outcome
# na.rm logical. should missing values be removed?
bin.ll <- function(p, y, na.rm = TRUE)
  sum(log(p[y]), na.rm = na.rm) + sum(log(1 - p[-y]), na.rm = na.rm)


############################## Heterogeneity, generalizability, pooled performance functions ###############################
# ### By convention, all generalizability measures:
# # Current required arguments:
# object  data.frame containing at least a column with 'estimate', and preferably more statistics: 
#             se, var, ci.lb     ci.ub measure  n   class
# # Required arguments: (OLD)
# x       list of class "listofperf", list of performance in different strata. 
#           Note that it has its own unlist method. Practically all (except plot) should call unlist first!
# ...     for compatibility.
#
# # Possible arguments, that are always passed through successfully:
# coef    data.frame containing coefficients of stratified models. Rows are strata, columns are coefs.
# coef.se data.frame containing se of coefficients of stratified models. Rows are strata, columns are coefs.
#
# 
# # Return:
# numeric of length 1.
##############################   

# Measure 0: mean
abs.mean <- function(object, ...)
  abs(mean(object$estimate)) 

## also possible the other way around (e.g. for cal slopes and intercepts)
mean.abs <- function(object, ...) 
  mean(abs(object$estimate))

# Measure 1: Coefficient of variation (=scaled sd)
# In general sense, abs needs not be TRUE, but for metapred it should,
# such that higher values are worse performance.
coef.var <- function(x, abs = TRUE, ...) {
  x <- unlist(x) 
  cv <- sd(x)/mean(x)
  if (isTRUE(abs)) abs(cv) else cv
}

# var.x.mean.with.se <- function(x, abs = TRUE, ...) {
#   x <- unlist(x) 
#   v <- var.with.se(x)
#   s2 <- v$estimate
#   vs2 <- v$variances
#   # s <- sqrt(v$estimate)
#   
#   m <- mean(x)
#   vm <- v$variances / length(x)
#   
#   est <- s2 * mean(x)
#   vest <- vs2 * vm + vs2 * m^2 + vm * s2^2
#   
#   out <- data.frame(estimate = est, variances = vest, se = sqrt(vest), n = length(x))
#   if (isTRUE(abs)) abs(out) else out
# }

coef.var.with.se <- function(object, abs = TRUE, ...) 
  cbind(data.frame(estimate = coef.var(object[["estimate"]])), bootstrap.se(object, coef.var))


bootstrap.se <- function(object, fun, k = 2000, ...) {
  x <- unlist(object[["estimate"]])
  fun <- match.fun(fun)
  
  est <- rep(NA, k)
  for (i in seq_len(k)) {
    est[i] <- fun(sample(x, size = length(x), replace = TRUE))
  }
  v <- var(est)
  out <- data.frame(variances = v, se = sqrt(v), n = length(x), k = k)
  class(out) <- c("mp.perf", class(out))
  out
}

coef.var.mean <- function(object, abs = TRUE, ...)  {
  x <- unlist(object[["estimate"]])
  coef.var(x, abs = abs) + if (abs) abs(mean(x)) else mean(x)
}
  

# Measure 2 (?): GINI coefficient
# #' @importFrom Hmisc GiniMd
GiniMd <- function(object, ...) 
  GiniMd(object[["estimate"]], na.rm = T)

# Also from Hmisc:
gmd <- function(object, ...) {
  x <- object[["estimate"]]
  n <- length(x)
  sum(outer(x, x, function(a, b) abs(a - b))) / n / (n - 1)
}

weighted.abs.mean <- function(object, ...) 
  abs(mean((object[["estimate"]] * object$n))) / sum(object$n)

# Fixed-Effects Meta-Analysis, Inverse Variance Method.
# NOTE: AUC is on wrong scale!!
fema <- function(object, ...) {
  # if (object$class[[1]] == "auc") {
  #   x <- logit(object$perf)
  #   v <- some other transform (object$var)
  # }
  # return(inv.logit(sum(unlist(x) / unlist(v)) / sum(1/unlist(v)) ))
  # else
  x <- object[["estimate"]]
  v <- object$var
  sum(unlist(x) / unlist(v)) / sum(1/unlist(v))
}

rema.beta <- rema.mean <- function(object, method = "REML", ...) 
  rema.perf(object, method = method, ...)$est

# valmeta does not produce tau!
# so rema.tau cannot be used on auc!
rema.tau <- function(object, method = "REML", ...)
  rema.perf(object, method = method, ...)$tau # Note: Intentionally selects tau2 if only that one is available.

# pooled.var <- function(x, n, ...) {
#   x <- unlist(x)
#   ## TODO: Extract sample size for each cluster and apply corresponding to the right performance measures
#   ## TODO: use rubins rules.
# }

pooled.var <- function(object, ...) {
  x <- unlist(object[["estimate"]])
  mean(x) + var(x) * (1 + 1/length(nrow(object)))
}

# squared.diff #a penalty equal to the mean squared differences 
squared.diff <- function(object, ...) {
  x <- unlist(object[["estimate"]])
  mse(x, mean(x))
}

# Mean of largest half of values
mean.of.large <- function(object, ...) {
  x <- unlist(object[["estimate"]])
  mean(x[x >= median(x)])
}

#  Forest plot of list of performance measures: AUC, cal intercept or slope, or mse/brier.
#' @importFrom metafor rma.uni
#' @export
plot.listofperf <- function(x, pfn, ...) { # xlab tbi from perfFUN
  # print("get perf name")
  # print("pfn:")
  # print(pfn)
  # print("...:")
  # print(list(...))
  # if (!is.null(pfn <- list(...)$pfn) && is.character(pfn)) {
    xlab <- paste(pfn, "in validation strata")
  # } else {
    # xlab <- "Performance in validation strata."  
  # }
  # print("get strata names")
  if (is.null(names(x))) # The # is to show users that the numbers are not their own. (no longer necessary)
    names(x) <- paste("#", seq_along(x), sep = "") 
  # print("compute ci now:")
  z <- ci.listofperf(object = x, ...)
  # z <<- z
  # print("meta-analyze performance:")
  
  # Thomas: I changed the implementation to uvmeta to ensure our prediction intervals are bsaed on Student T distribution
  # and to ensure we are using REML everywhere.
  if (inherits(x[[1]], "auc")) { # To be replaced by child function.
    # print("by valmeta")
    ma <- valmeta(measure = "cstat", cstat = z$theta, cstat.cilb = z[,"theta.ci.lb"], 
                  cstat.ciub=z[,"theta.ci.ub"], cstat.cilv=0.95, method = "REML")
    est <- ma$est
    pi.lb <- ma$pi.lb
    pi.ub <- ma$pi.ub
  } else if (inherits(x[[1]], c("lm"))) { # 
    # print("by rma.uni")
    ma <- uvmeta(r = sapply(x, coef), r.vi = sapply(x, variances), method = "REML") 
    est <- ma$est
    pi.lb <- ma$pi.lb
    pi.ub <- ma$pi.ub
  } else if (inherits(x[[1]], "mse")) {
  #   # print("by rma.uni")
    ma <- uvmeta(r = sapply(x, `[[`, "estimate"), r.vi = sapply(x, variances), method = "REML")
    est <- ma$est
    pi.lb <- ma$pi.lb
    pi.ub <- ma$pi.ub
  }
  # This is the same for both methods:
  ci.lb <- ma$ci.lb
  ci.ub <- ma$ci.ub

  # print("Make forest plot.")
  fp <- metamisc::forest(theta       = z$theta,
                         theta.ci.lb = z$theta.ci.lb,
                         theta.ci.ub = z$theta.ci.ub,
                         theta.slab  = names(x),
                         theta.summary       = est,
                         theta.summary.ci.lb = ci.lb,
                         theta.summary.ci.ub = ci.ub,
                         theta.summary.pi.lb = pi.lb,
                         theta.summary.pi.ub = pi.ub,
                         xlab  = xlab,
                         ...)
  plot(fp)
  invisible(NaN)
}

# x mp.cv.val object
# y ignored, compatibility only
#' @export
plot.mp.cv.val <- function(x, y, ...)
  plot.listofperf(x$perf.full, x$perf.name, ...)

rema.perf <- function(object, method = "REML", ...) {
  if (object$class[[1]] == "mp.perf" || object$class[[1]] == "recal") {
    ma <- uvmeta(r = object[["estimate"]], r.vi = object$var, method = method) # uvmeta uses a Student T distribution, in contrast to metafor
    return(list(est = ma$est,     
                pi.lb = ma$pi.lb,
                pi.ub = ma$pi.ub,
                ci.lb = ma$ci.lb,
                ci.ub = ma$ci.ub,
                tau2  = ma$tau2,
                tau   = sqrt(ma$tau2)))
  } else if (object$class[[1]] == "auc") {
    ma <- valmeta(measure = "cstat", cstat = object[["estimate"]], 
                  cstat.cilb = object[,"ci.lb"], cstat.ciub = object[,"ci.ub"],
                  cstat.cilv = 0.95, method = method)
    return(list(est = ma$est,
                pi.lb = ma$pi.lb,
                pi.ub = ma$pi.ub,
                ci.lb = ma$ci.lb,
                ci.ub = ma$ci.ub)) # valmeta does not produce tau! (but can be obtained from ma$fit if "ret.fit=T")
  }
  stop("class not recognized")
}

rema.mp.cv.val <- function(object, method = "REML", ...)
  rema.perf(object[["perf"]], method = method)



#' Forest plot of a metapred fit
#' 
#' Draw a forest plot of the performance of an internally-externally cross-validated model. By default the final model is shown.
#' 
#' @author Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
#' 
#' @param object A \code{metapred} fit object.
#' @param step Which step should be plotted? Defaults to the best step. numeric is converted to name of the step: 0 for 
#' an unchanged model, 1 for the first change...
#' @param model Which model change should be plotted? NULL (default, best change) or character name of variable or (integer) 
#' index of model change.
#' @param perfFUN Numeric or character. Which performance statistic should be plotted? Defaults to the first.
#' @param method character string specifying whether a fixed- or a random-effects model should be used to summarize the
#' prediction model performance. A fixed-effects model is fitted when using method="FE". Random-effects models are fitted 
#' by setting method equal to one of the following: "DL", "HE", "SJ", "ML", "REML", "EB", "HS", or "GENQ". Default is "REML".
#' @param ...  Other arguments passed to plotting internals. E.g. \code{title}. See \link{forest.default} for details.
#' 
#' @examples 
#' data(DVTipd)
#' 
#' # Internal-external cross-validation of a pre-specified model 'f'
#' f <- dvt ~ histdvt + ddimdich + sex + notraum
#' fit <- metapred(DVTipd, strata = "study", formula = f, scope = f, family = binomial)
#' 
#' # Display the model's external performance (expressed as mean squared error by default) 
#' # for each study
#' forest(fit)
#' 
#' @export
forest.metapred <- function(object, perfFUN = 1, step = NULL, method = "REML", model = NULL, ...)
  forest.mp.cv.val(subset(object, step = step, model = model), perfFUN = perfFUN, method = method, ...)


#' Forest plot of a validation object.
#' 
#' Draw a forest plot of the performance of an internally-externally cross-validated model.
#' 
#' @author Valentijn de Jong <Valentijn.M.T.de.Jong@gmail.com>
#'  
#' @param object An \code{mp.cv.val} or \code{perf} object.
#' @param perfFUN Numeric or character. Which performance statistic should be plotted? Defaults to the first.
#' @param method character string specifying whether a fixed- or a random-effects model should be used to summarize the
#' prediction model performance. A fixed-effects model is fitted when using method="FE". Random-effects models are fitted 
#' by setting method equal to one of the following: "DL", "HE", "SJ", "ML", "REML", "EB", "HS", or "GENQ". Default is "REML".
#' @param ... Other arguments passed to plotting internals. E.g. \code{title}. See \link{forest.default} for details.
#' 
#' @aliases forest.mp.cv.val forest.perf
#' 
#' @method forest mp.cv.val
#' 
#' @export
# forest.mp.cv.val <- function(object, statistic = 1, ...)
#   forest.perf(object[["perf.all"]][[statistic]],
#               xlab = if (is.character(statistic)) statistic else
#                 object[["perf.names"]][[statistic]], ...)

forest.mp.cv.val <- function(object, perfFUN = 1, method = "REML", ...) 
  forest.perf(perf(object, perfFUN = perfFUN, ...),
              xlab = if (is.character(perfFUN)) perfFUN else
                object$perf.names[[perfFUN]], method = method, ...)

forest.perf <- function(object, method = "REML", ...) {
  if (is.null(theta.slab <- list(...)$theta.slab))
    theta.slab <- as.character(object$val.strata)
  ma <- rema.perf(object, method = method)
  fp <- metamisc::forest(theta       = object[["estimate"]],
                         theta.ci.lb = object$ci.lb,
                         theta.ci.ub = object$ci.ub,
                         theta.slab  = theta.slab,
                         theta.summary       = ma$est,
                         theta.summary.ci.lb = ma$ci.lb,
                         theta.summary.ci.ub = ma$ci.ub,
                         theta.summary.pi.lb = ma$pi.lb,
                         theta.summary.pi.ub = ma$pi.ub,
                         ...)
  plot(fp)
  invisible(NaN) # To be replaced with fp, when metapred() can handle it.
}
# 
# sampleBinary <- function(n = 50, J = 1, b = rep(log(2), J), alpha = NULL, col.names = NULL ) {
#   J <- length(b)
#   if (is.null(alpha)) alpha <- -log(sqrt(prod(exp(b))))
#   if (is.null(col.names)) col.names <- c("Y", paste("X", 1:J, sep = ""))
#   coefficientss <- c(alpha, b)
#   x  <- cbind(1, matrix(rbinom(n * J, size = 1, prob = .5), nrow = n, ncol = J))
#   lp <- coefficientss %*% t(x)
#   p  <- metamisc:::inv.logit(lp)
#   y  <- stats::rbinom(length(lp), size = 1, prob = p)
# 
#   out <- data.frame(cbind(y,x[ , -1]))
#   colnames(out) <- col.names
#   out
# }
# 
# da <- sampleBinary(n = 2000, J = 3)
# da$X3[sample(seq_len(nrow(da)), size = nrow(da)/2, replace = F)] <- 3
# 
# mp <- metapred(da, strata = "X3", family = binomial,
#                perfFUN = list(mse = "mse", auc = "auc", cal.slope = "cal.slope"))
# cv <- subset(mp)
# # perf0 <- cv$perf.all[[2]]
# 
# forest(cv, statistic = 3)
# forest(mp, statistic = "mse", title = "Hallo")
# # metamisc:::forest.perf(perf0)

fat.perf <- function(object, ...)
  fat(b = object[["estimate"]], b.se = object[["se"]], n.total = sum(object[["n"]]), ...)

fat.mp.cv.val <- function(object, ...) 
  fat(object[["perf"]], ...)

fat.metapred <- function(object, ...)
  fat.mp.cv.val(subset(object, ...))

funnel.perf <- function(object, ...)
  plot(fat.perf(object, ...), ...)

funnel.mp.cv.val <- function(object, ...)
  plot(fat.mp.cv.val(object))

funnel.metapred <- function(object, ...)
  plot(fat.metapred(object, ...))