### To add / change:
# is.metapred()
# variances for intercept recalibration.
# One-stage MA: predStep: is  type = response correct?
# perf method
# penalties
# performance measurement, current is Brier for binomial
# add more options to penalty step
# , cl.name in modelStep
# Remove Reduce() from perfStep, and make perfStep() compatible with multiple data sets.
# This means that cvFUN = fixed does not work correctly either.

# Changed:
# Added more robust tests for stepwise.

#' Generalized Stepwise Regression for prediction models
#'
#' Generalized stepwise regression for obtaining a prediction model with adequate performance across data sets. Requires
#' data from individuals in multiple studies.
#' 
#' @author Valentijn de Jong
#' 
#' @references Debray TPA, Moons KGM, Ahmed I, Koffijberg H, Riley RD. A framework for developing, implementing, 
#' and evaluating clinical prediction models in an individual participant data meta-analysis. 
#' \emph{Stat Med}. 2013;32(18):3158-80. 
#'
#' @param data data.frame containing the datasets.
#' @param strata Name of the strata (e.g. studies or clusters) variable, as character. Used for two-stage MA only.
#' @param formula Formula of the full model to be evaluated, and possibly reduced. If not supplied,
#' it is assumed the first column in the data set is the outcome, and all remaining columns
#' (except \code{strata}) are predictors. See \link[stats]{formula} for details.
#' @param estFUN Function for estimating the model in the first stage. Currently "lm" and "glm" are supported.
#' @param stepwise Logical. Should stepwise selection be performed?
#' @param center.out Logical. Should the outcome be centered within studies?
#' @param center.cov Logical. Should covariates be centered within studies?
#' @param recal.int Logical. Should the intercept be recalibrated?
#' @param cvFUN Cross-validation method, on the study (i.e. cluster or stratum) level. "
#' l1o" for leave-one-out cross-validation (default). "bootstrap" for bootstrap. Or "fixed", for one or more data sets
#' which are only used for validation. A user written function may be supplied as well.
#' @param cv.k Parameter for cvFUN. For \code{cvFUN="bootstrap"}, this is the number of bootstraps. For \code{cvFUN="fixed"}, 
#' this is a vector of the indices of the (sorted) data sets.
#' @param metaFUN Function for computing the meta-analytic coefficient estimates in two-stage MA. Default: \link[metafor]{rma}
#' from the metafor package is used. Default settings are univariate random effects, estimated with "REML". Method can be
#' passed trough the \code{meta.method} argument.
#' @param meta.method Name of method for meta-analysis. Default is "REML". For more options see \link[metafor]{rma}.
#' @param predFUN Function for predicting new values. Defaults to the appropriate link functions for two-stage MA where
#' \code{glm()} or \code{lm()} is used in the first stage. For one-stage models \code{predict()} is used.
#' @param perfFUN Function for computing the performance of the prediction models. Default: mean squared error.
#' @param genFUN Function computing generalizability measure using the performance measures. Default: (absolute) mean.
#' \code{squareddiff} for a penalty equal to the mean squared differences between coefficients.
#' @param selFUN Function for selecting the best method. Default: lowest value for \code{genFUN}. Should be set to
#' "which.max" if high values for \code{genFUN} indicate a good model.
#' @param ... To pass arguments to estFUN (e.g. family = "binomial"), or other methods.
#'
#' @return \code{metapred} A list of class \code{metapred}, containing the final coefficients in \code{coefficients}, and the stepwise
#' tree of estimates of the coefficients \code{(coef)}, performance measures \code{(perf)}, generalizability measures
#' \code{(gen)} in \code{stepwise}, and more.
#'
#' @examples 
#' data(DVTipd)
#' DVTipd$cluster <- 1:4 # Add a fictional clustering to the data set.
#' metamisc:::metapred(DVTipd, strata = "cluster", f = dvt ~ sex + vein + malign, family = binomial)
#' 
#'\dontrun{
#' # Some additional examples:
#' metamisc:::metapred(DVTipd, strata = "cluster", f = dvt ~ sex + vein + malign
#' , family = binomial, stepwise = FALSE)
#' metamisc:::metapred(DVTipd, strata = "cluster", f = dvt ~ sex + altdiagn + histdvt
#' , family = binomial, recal.int = TRUE)
#' metamisc:::metapred(DVTipd, strata = "cluster", f = dvt ~ sex + altdiagn + histdvt
#' , family = binomial, meta.method = "DL")
#'}
#' # By default, metapred assumes the first column is the outcome.
#' DVTipd.reordered <- DVTipd[c("dvt", "ddimdich", "histdvt", "cluster")]
#' mp <- metamisc:::metapred(DVTipd.reordered, strata = "cluster", family = binomial)
#' fitted <- predict(mp, newdata = DVTipd.reordered)
#' 
#' 
#' @import stats
#'
#' @importFrom stats formula
#'
#' @export

metapred <- function(data, strata, formula = NULL, estFUN = "glm", stepwise = TRUE, center.out = FALSE,
                     center.cov = FALSE, recal.int = FALSE, cvFUN = NULL, cv.k = NULL,
                     metaFUN = NULL, meta.method = "REML", predFUN = NULL, perfFUN = NULL, genFUN = NULL, selFUN = "which.min",
                     ...) {
  call   <- match.call()
  data   <- as.data.frame(data)
  estFUN <- match.fun(estFUN)
  two.stage <- TRUE

  if (is.null(formula)) formula <- stats::formula(data[ , -which(colnames(data) == strata)])
  strata.i  <- data[, which(colnames(data) == strata)]
  data      <- stats::model.frame(formula, data = data)
  data      <- centerData(data, center.in = strata.i, center.1st = center.out, center.rest = center.cov) # change for 1 stage?
  data.list <- asDataList(data, strata.i)

  if (is.null(cvFUN))   cvFUN   <- "l1o"
  if (is.null(metaFUN)) metaFUN <- "urma"
  if (is.null(perfFUN)) perfFUN <- "mse"
  if (is.null(genFUN))  genFUN  <- "absmean"
  # Change to "-" when perfFUN <- mcfadden or some other measure for which greater = better.

  # Do not add metaFUN to this list!
  # probably not predFUN either.
  cvFUN.name <- if (!missing(cvFUN) && is.character(cvFUN)) cvFUN else "cvFUN"
  estFUN.name <- estFUN
  cvFUN   <- get(cvFUN)
  perfFUN <- get(perfFUN)
  genFUN  <- get(genFUN)
  selFUN  <- get(selFUN)

  J <- ncol(data) - 1 # because 1 is outcome. For surv will have to be 2
  ccs <- seq_len(J) + 1
  b <- v <- b.recal <- perf <- gen <- list()
  best.gen <- NULL
  folds <- cvFUN(1:length(data.list), k = cv.k)
  if (!isTRUE(length(folds$dev) > 0) || !isTRUE(length(folds$dev[[1]]) > 0))
    stop("At least 1 cluster must be used for development.")
  family <- NULL

  for (j in (seq_len(ncol(data)) - 1)) {

    fold.b <- fold.v <- fold.perf <- fold.b.recal <- list()

    for (fold.i in seq_len(length(folds$dev)))
    {
      # Generate new model(s)
      step <- modelStep(data.list = data.list, ccs = ccs, estFUN = estFUN, perfFUN = perfFUN,
                        metaFUN = metaFUN, meta.method = meta.method, drop = j, cl = folds$dev.i[[fold.i]],
                        cl.name = folds$dev[[fold.i]], ...)

      fold.b[getFoldName(ds = folds$dev[[fold.i]], f = fold.i, type = cvFUN.name)] <- list(step$meta.b)
      fold.v[getFoldName(ds = folds$dev[[fold.i]], f = fold.i, type = cvFUN.name)] <- list(step$meta.v)

      if (is.null(family) && (!is.null(step$dummy.model$family)))
        family <- step$dummy.model$family

      # Compute performance measures
      step.perf <- perfStep(newdata = Reduce(rbind, data.list[folds$val[[fold.i]]]), b = step$meta.b,
                            fit = step$dummy.model, two.stage = two.stage, ccs = step$step.ccs, f = step$f,
                            recal.int = recal.int, estFUN = estFUN, predFUN = predFUN, perfFUN = perfFUN, ...)

      fold.perf[   getFoldName(ds = folds$val[[fold.i]], f = fold.i, type = cvFUN.name)] <- list(step.perf$perf)
      fold.b.recal[getFoldName(ds = folds$val[[fold.i]], f = fold.i, type = cvFUN.name)] <- list(step.perf$b)
    }

    b[      getStepName(j)] <- list(fold.b)
    v[      getStepName(j)] <- list(fold.v)
    b.recal[getStepName(j)] <- list(fold.b.recal)

    fold.perf.df <- t(data.frame(fold.perf))


    # row.names(fold.perf.df) <- sapply(folds$val, getFoldName)
    step.gen <- apply(fold.perf.df, MARGIN = 2, FUN = genFUN)
    gen[getStepName(j)]  <- list(step.gen)
    perf[getStepName(j)] <- list(fold.perf.df)

    # Select a model
    selected.model <- selFUN(unlist(c(best.gen, step.gen)))

    if (!stepwise) break
    if (j) {
      if (!(selected.model - 1)) break
      ccs      <- ccs[-(selected.model - 1)]
      best.gen <- step.gen[[selected.model - 1]]
    } else     best.gen <- step.gen[[1]]
  }
  # End loop covariate selection. = up to one "step" for each covariate

  coefficients.recal <- if (recal.int) b.recal else NULL

  # Generate final model
  final.predictors    <- getCovariateNames(data.list, ccs)
  final.formula.names <- getFormula(data, ccs)                       # To select columns from a model.frame
  final.formula       <- removeFormulaBackticks(final.formula.names) # To create a model.frame
  final.model         <- modelStep(data.list, ccs = ccs, estFUN = estFUN, perfFUN = perfFUN,
                                   metaFUN = metaFUN, meta.method = meta.method,
                                   cl = NULL, drop = FALSE, ...)
  final.b             <- final.model$meta.b[[1]]
  final.v             <- final.model$meta.v[[1]]

  out <- list(stepwise = list(gen = gen, perf = perf, coefficients = b, n.steps = j, v = v,
                              coefficients.recal = coefficients.recal),
              # final = list(
              gen = best.gen, predictors = final.predictors, formula = final.formula,
              formula.names = final.formula.names, coefficients = final.b, variance = final.v, #),
              call = call,
              family = family,
              FUN = list(cvFUN = cvFUN, perfFUN = perfFUN, metaFUN = metaFUN, genFUN = genFUN,
                         selFUN = selFUN, predFUN = step.perf$predictMethod, estFUN = estFUN.name),
              options = list(cv.k = cv.k, meta.method = meta.method, recal.int = recal.int,
                             center.cov = center.cov, center.out = center.out)
  )

  class(out) <- c("metapred")
  return(out)
}

#' @export
print.metapred <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call: ");                       print(x$call); cat("\n")
  cat("Steps attempted:",              x$stepwise$n.steps, "\n")
  cat("Final predictors:",             x$predictors, "\n")
  cat("Final generalisability value:", x$gen, "\n"); cat("\n")

  print(round(x$coefficients, digits = digits))
}

perfStep <- function(newdata, b, fit, two.stage, ccs = rep(list(1:ncol(newdata), ncol(b))), f = formula(newdata),
                     recal.int = FALSE, estFUN = NULL, predFUN = NULL, perfFUN = NULL, ...) {
  b.original <- b <- as.list(b)

  if (length(b) > 0) {
    if (isTRUE(recal.int))
      for (i in 1:length(b))
        b[[i]][1] <- b[[i]][1] + computeRecal(object = fit, newdata = newdata, b = b[[i]], estFUN = estFUN)[1]
  } else
    stop("length(b) must be > 0; to estimate performance.")


  predictMethod <- getPredictMethod(fit = fit, two.stage = two.stage, predFUN = predFUN, ...)
  p <- list()
  for (i in 1:length(b)) {
    p[[i]] <- predictMethod(object = fit, newdata = newdata, type = "response", b = b[[i]],
                            f = f[[i]], two.stage = two.stage, ...)
  }

  perf <- sapply(p, FUN = perfFUN, y = newdata[ , 1], data = newdata, fit = fit)
  names(perf) <- rownames(b)

  out <- list(perf = perf, predictMethod = predictMethod, b = b, b.original = b.original, recal.int = recal.int)
  out
}

mse <- function(p, y, data = NULL, ...) mean((p - y)^2)

absmean <- function(perf.measures, ...) {
  pm <- unlist(perf.measures)
  abs(mean(pm))
  }

squareddiff <- function(perf.measures, ...) {
  pm <- unlist(perf.measures)
  abs(mean(pm)) + mean((mean(pm) - pm)^2)
  }

# Gets the predict method.
# fit Model fit object.
# two.stage logical. Is the model a two-stage model?
# predFUN Optional function, which is immediately returned
# ... For compatibility only.
getPredictMethod <- function(fit, two.stage, predFUN = NULL, ...) {
  # A user written function may be supplied:
  if (!is.null(predFUN))
    if (is.function(predFUN))
      return(get(predFUN))
  else stop("predFUN should be a function.")

  # If two-stage, the fit is used only to extract the link function.
  # If one-stage, fit's prediction method may be used.
  if (two.stage) {
    if (inherits(fit, c("glm", "lm")))
      return(predictGLM)
    else stop("No prediction method has been implemented for this model type yet for two-stage
              meta-analysis. You may supply one with the predFUN argument.")
  } else return(predict)
  }

# Prediction function for two-stage GLM objects
# object glm model fit object
# newdata newdata to predict for, of class "data.frame"
# b vector of coefficients. Overrides coefficients of object
# f formula used for selecting relevant variables from newdata.
# ... For compatibility only.
predictGLM <- function(object, newdata, b = NULL, f = NULL,  ...) {
  if (is.null(b)) b <- coef(object)
  if (is.null(f)) f <- formula(object)
  f <- stats::as.formula(f)
  X <- as.matrix(stats::model.frame(f, data = newdata))
  X[ , 1] <- 1
  lp <- as.matrix(X) %*% b
  if (is.null(object$family)) lp
  else object$family$linkinv(lp)
}

# For prediction of newdata. Not used internally.
# object Model fit object
# newdata New data set of class "data.frame"
# type Currently unused name of prediction type.
# ... For compatibility only.
#' @export
predict.metapred <- function(object, newdata = NULL, type = "response", recal.int = FALSE, ...)
{
  if (isTRUE(is.null(newdata)))
    stop("A newdata argument should be supplied. The colnames should match variable names of the metapred object.")

  if (isTRUE(recal.int))
    object <- recalibrate(object = object, newdata = newdata)

  object$FUN$predFUN(object = object, newdata = newdata, type = type, ...)
}

modelStep <- function(data.list, ccs, estFUN, perfFUN, metaFUN, meta.method, genFUN, cl = NULL, cl.name = NULL, drop = TRUE, ...) {
  if (isTRUE(is.null(cl))) cl <- seq_len(length(data.list))
  if (isTRUE(is.null(cl.name))) cl.name <- as.character(cl)
  metaFUN <- match.fun(metaFUN)
  step.b <- step.v <- step.covar <- step.ccs <- step.f <- step.meta.b <- step.meta.v <- list()

  if (drop && (is.null(ccs) || length(ccs) == 0))
    stop("model already empty. Nothing left to  remove.")

  ccs.test <- if (drop) ccs else Inf # Inf leads to the full model
  for (ccs.ex in ccs.test) {
    ccs.t <- ccs[which(!ccs.ex == ccs)]
    study.b <- study.v <- study.covar <- list()
    f <- getFormula(data.list[[1]], ccs.t)

    for (s in cl) {
      study.model <- estFUN(formula = f, data = data.list[[s]], ...)
      study.b[[getFoldName(cl.name)]]     <- getCoefs(study.model)
      study.v[[getFoldName(cl.name)]]     <- getVars(study.model)
      study.covar[[getFoldName(cl.name)]] <- getCoVars(study.model)
    }
    step.b[getModelName(data.list, ccs.t)]     <- list(study.b)
    step.v[getModelName(data.list, ccs.t)]     <- list(study.v)
    step.covar[getModelName(data.list, ccs.t)] <- list(study.covar)
    step.ccs[getModelName(data.list, ccs.t)]   <- list(ccs.t)
    step.f[getModelName(data.list, ccs.t)]     <- list(f)


    b <- as.data.frame(t(as.data.frame(study.b)))
    v <- as.data.frame(t(as.data.frame(study.v)))

    step.meta <- metaFUN(b = b, v = v, method = meta.method) # , covar = step.covar)
    step.meta.b[getModelName(data.list, ccs.t)]  <- list(step.meta$b)
    step.meta.v[getModelName(data.list, ccs.t)]  <- list(step.meta$v)
  }

  return(list(meta.b = step.meta.b, meta.v = step.meta.v,
              first.stage.b = step.b, first.stage.v = step.v, first.stage.covar = step.covar,
              f = step.f, step.ccs = step.ccs,
              dummy.model = study.model))
}

# Univariate Random Effects Meta-Analysis
# b data.frame or matrix, containing coefficients
# v data.frame or matrix, containing variances
# method Method for meta-analysis.
# ... Optional arguments for rma().
#' @importFrom metafor rma
urma <- function(b, v, method = "REML", ...)
{
  if (!(is.data.frame(b) || is.matrix(b)) || !(is.data.frame(v) || is.matrix(v)) )
    stop("b and v must both be a data.frame or matrix.")
  if (!identical(dim(b), dim(v)))
    stop("b and v must have the same dimensions.")

  meta.b <- meta.se <- rep(NA, ncol(b))
  for (col in 1:ncol(b)) {
    r <- metafor::rma(b[ , col] , v[ , col], method = method, ...)
    meta.b[col]  <- r$beta
    meta.se[col] <- r$se
  }

  meta.v <- meta.se^2

  names(meta.b) <- names(meta.v) <- names(meta.se) <- colnames(b)
  list(b = meta.b, v = meta.v, se = meta.se)
}

#' @export
coef.metapred <- function(object, ...)
  object$coefficients

#' @export
family.metapred <- function(object, ...)
  object$family

#' @export
formula.metapred <- function(x, ...)
  x$formula
