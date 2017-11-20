# Internal function for centering
# x numeric vector to be centered in data sets.
# center.in indices of data sets.
center <- function(x, center.in) {
  if (length(center.in) != length(x))
    stop("length(center.in) should match length(x).")
  for (trial in sort(unique(center.in)))
  {
    selection.id <- center.in == trial
    selection <- x[selection.id]
    x[selection.id] <- selection - mean(selection, na.rm = T)
  }
  x
}

# Centers data within studies / clusters
# data data.frame. data set.
# center.in numeric vector corresponding to cluster indices.
# center.1st logical. Should the 1st variable (response) be centered
# center.rest logical. Should the other variables (predictors) be centered?
centerData <- function(data, center.in, center.1st = FALSE, center.rest = FALSE) {
  if (!is.data.frame(data) && !is.matrix(data))
    stop("data should be a data.frame or matrix.")
  if (length(center.in) != nrow(data))
    stop("length(center.in) should match nrow(data).")
  if (isTRUE(center.1st))
    data[ , 1] <- center(data[ , 1], center.in)
  if (isTRUE(center.rest) && ncol(data) > 2)
    for (col in 2:ncol(data))
      data[ , col] <- center(data[ , col], center.in)
    data
}

# coerces data set to data.list
# data data.frame. data set.
# strata.i numeric. stratum indicators.
# Returns data.list
asDataList <- function(data, strata.i) {
  data.list <- list()
  strata <- sort(unique(strata.i))
  for (i in 1:length(strata))
    data.list[[i]] <- data[strata.i == strata[i], ]
  names(data.list) <- strata
  data.list
}

# gets only the relevant data from a data.list
# data.list list of data sets
# ccs numeric. covariate column selection
# cl numeric. indices of clusters to be selected.
# returns data.list.
getDataList <- function(data.list, ccs, cl) 
  lapply(data.list[cl], getData, predictors = ccs)




# These functions are used for making various names. Any change should be made here, such that
# these functions can also be used to retrieve objects from a list.
# ds vector. indices of datasets.
# f numeric. fold index.
# type character. type of cv (just a chosen name)
# data.list list of data sets.
# data data.frame.
# covariate.columns indices of covariate columns.
# All return character.
getFoldName <- function(ds, f = NULL, type = NULL)
  paste(getclName(ds = ds), getcvName(f = f, type = type), sep = if (length(f) > 0 || length(type) > 0) ". " else "")

getclName <- function(ds)
  paste("cl", toString(ds), sep = " ")

getcvName <- function(f, type = NULL)
  paste(type, f, sep = " ")

getCovariateNames <- function(data.list, covariate.columns) {
  if (!length(covariate.columns))  return(NULL)
  if (!is.null(colnames(data.list[[1]]))) return(colnames(data.list[[1]])[covariate.columns])
  warning("Covariate names could not be found.")
  return(paste("X", covariate.columns - 1, sep = ""))
}

getModelName <- function(data.list, covariate.columns)
  if (is.null(name <- getCovariateNames(data.list, covariate.columns) ))
    return("intercept only") else return(toString(name))

getStepName <- function(x)
  paste("step", x, sep = "")

getPredictorNames <- function(f, data)
  colnames(stats::model.frame.default(formula = f, data = data))[-1]

# This one is mostly important for getting the intercept-only formula:
# data data.frame
# predictors indices of predictors
# returns formula
getFormula <- function(data, predictors = NULL) {
  if (length(predictors) == 0 || identical(predictors, 0)) {
    f <- stats::formula(data)
    f[3] <- 1 # [3] is the right hand side.
  } else f <- stats::formula(data[ , c(1, predictors)])
  stats::as.formula(f)
}

# Gets only the requested data, as well as the response variable.
# data data.frame or matrix
# pred.indices indexes of preditor columns
# Returns data.frame.
getData <- function(data, pred.indices) {
  if (!is.data.frame(data) && !is.matrix(data))
    stop("data must be a data.frame or matrix.")
  sel <- c(1, pred.indices)
  d <- data.frame(data[ , sel])
  names(d) <- names(data)[sel]
  d
}


# gets formula of for data set. Replacement of formula.data.frame, which returns bogus for
# 1 column data.frames.
# data data.frame
# Returns formula.
getFullFormula <- function(data) {
  data <- as.data.frame(data)
  # if (!is.data.frame(data) && !is.matrix(data))
  #   stop("data must be a data.frame or matrix.")
  
  if (identical(ncol(data), 1L))
  {
    xnames <- names(data)[-1]
    if (identical(length(xnames), 0L))
      xnames <- "1"
    left <- paste(names(data)[1], "~")
    return(paste(left, xnames))
    
  } else return(stats::formula(data))
}

getCoefs  <- function(fit, ...) {
  if (inherits(fit, "multinom"))
    return(coefMultinom(fit, ...))
  else return(coef(fit))
}
# Needs some work. Should also return some coefficient names.
coefMultinom <- function(fit, ...)
  as.vector(t(coef(fit)))

# Perhaps unnecessary:
getVars   <- function(fit, ...) diag(vcov(fit))
getCoVars <- function(fit, ...) vcov(fit)
getSE     <- function(fit, ...) sqrt(getVars(fit))


# Coerces l to one string
# l list or vector of strings.
# returns 1 string, character.
oneStr <- function(l, sep = "") {
  if (length(l) > 0) out <- l[[1]] else return("")
  
  if (length(l) > 1) {
    for (i in 2:length(l)) out <- paste(out, l[[i]], sep = sep)
  }
  return(out)
}


# Necessary for making a formula from a model.frame.
# f formula
# returns: formula without backticks.
removeFormulaBackticks <- function(f)
{
  g <- gsub("`", "", f)
  h <- oneStr(list(g[[2]], " ~ ", g[[3]]), sep = "")
  stats::update.formula(f, h)
}

# f formula
# p.name character. name of predictor
# Returns: formula
removePredictor <- function(f, p.name)
  stats::update.formula(f, paste(". ~ . -", p.name) )

# f formula
# p.name list or vector of character names of predictors
# Returns: formula
removePredictors <- function(f, p.names)
  removePredictor(f, paste(oneStr(unlist(p.names), sep = " - "), sep = " - "))


### The following functions are for generating the folds for the cross-validation in metapred
# ds Numeric or character vector. Names of the data sets.
# k Numeric. Differs per function
#   bootstrap: number of bootstraps
#   fixed: indices of data sets for validation.
#   successive (still a hidden function): Number of validation sets.
#   leaveOneOut: ignored.
leaveOneOut <- l1o <- function(ds, ...) {
  ds <- sort(ds)
  if (length(ds) < 2)
    stop("iecv not possible for less than 2 data sets.")
  indexes <- seq_len(length(ds))
  out <- list(dev = list(), dev.i = list(), val = as.list(ds[indexes]), val.i = as.list(indexes))
  
  for (i in indexes) {
    out$dev[[i]]   <- ds[-indexes[i]]
    out$dev.i[[i]] <- indexes[-i]
  }
  
  out
}

bs <- bootstrap <- function(ds, k = NULL, ...) {
  ds <- sort(ds)
  if (is.null(k))
    k <- 200
  if (length(ds) < 2)
    stop("Bootstrapping data sets is impossible for < 2 data sets.")
  dev <- dev.i <- val <- val.i <- list()
  
  i <- 1
  while (length(dev) < k) {
    indexes <- sample(1:length(ds), replace = T)
    if (length(unique(indexes)) >= length(ds))
      next
    dev[[i]] <- ds[indexes]
    val[[i]] <- ds[-indexes]
    dev.i[[i]] <- indexes
    val.i[[i]] <- seq_len(length(ds))[-indexes]
    i <- i + 1
  }
  list(dev = dev, dev.i = dev.i , val = val, val.i = val.i)
}

fixed <- function(ds, k = NULL, ...)
{
  ds <- sort(ds)
  if (length(ds) < 2)
    stop("Selecting a validation set is impossible for < 2 data sets.")
  if (is.null(k))
    k <- length(ds)
  indexes <- seq_len(length(ds))
  list(dev = list(ds[-k]), dev.i = list(indexes[-k]), val = list(ds[k]), val.i = list(indexes[k]))
}

successive <- function(ds, k = NULL, ...) {
  ds <- sort(ds)
  if (is.null(k)) k <- 1
  k <- as.integer(k)
  if (k < 1) stop("k must be >= 1")
  if (length(ds) < (k + 1)) stop("Cross-validation requires k + 1 data sets, where k is the number of test sets.")
  out <- list(dev = list(), dev.i = list(), val = list(), val.i = list())
  
  for (i in seq_len(length(ds) - k) ) {
    sel <- 1:i
    out$dev.i[[i]] <- sel
    out$dev[[i]] <- ds[sel]
  }
  for (i in seq_len(length(out$dev)) ) {
    sel <- (i + 1):(i + k)
    out$val.i[[i]] <- sel
    out$val[[i]] <- ds[sel]
  }
  out
}
