#' Calculate the total O:E ratio
#'
#' This function calculates (transformed versions of) the ratio of total number of observed versus expected events with the 
#' corresponding sampling variance. 
#' 
#' @param OE vector with the estimated ratio of total observed versus total expected events
#' @param OE.se vector with the standard errors of the estimated O:E ratios
#' @param OE.cilb vector to specify the lower limits of the confidence interval for \code{OE}.
#' @param OE.ciub vector to specify the upper limits of the confidence interval for \code{OE}.
#' @param OE.cilv vector to specify the levels of aformentioned confidence interval limits. 
#' (default: 0.95, which corresponds to the 95\% confidence interval).
#' @param citl vector with the estimated calibration-in-the-large statistics
#' @param citl.se vector with the standard error of the calibration-in-the-large statistics
#' @param N vector to specify the sample/group sizes.
#' @param O vector to specify the total number of observed events.
#' @param E vector to specify the total number of expected events
#' @param Po vector to specify the (cumulative) observed event probabilities.
#' @param Po.se vector with the standard errors of \code{Po}.
#' @param Pe vector to specify the (cumulative) expected event probabilites
#' (if specified, during time \code{t.val})
#' @param data optional data frame containing the variables given to the arguments above.
#' @param slab optional vector with labels for the studies.
#' @param add a non-negative number indicating the amount to add to zero counts. See `Details'
#' @param g a quoted string that is the function to transform estimates of the total O:E ratio; see the details below.
#' @param level level for confidence interval, default \code{0.95}.
#' @param \ldots Additional arguments.
#' 
#' @examples 
#' ######### Validation of prediction models with a binary outcome #########
#' data(EuroSCORE)
#' 
#' # Calculate the total O:E ratio and its standard error
#' oecalc(O=n.events, E=e.events, N=n, data=EuroSCORE, slab=Study)
#' 
#' # Calculate the log of the total O:E ratio and its standard error
#' oecalc(O=n.events, E=e.events, N=n, data=EuroSCORE, slab=Study, g="log(OE)")
#' 
#' @keywords meta-analysis calibration performance
#' 
#' @author Thomas Debray <thomas.debray@gmail.com>
#' 
#' @export
oecalc <- function(OE, OE.se, OE.cilb, OE.ciub, OE.cilv, citl, citl.se, N, O, E, Po, Po.se, Pe, 
                   data, slab, add=1/2, g=NULL, level=0.95, ...) {
  
  ### check if data argument has been specified
  if (missing(data))
    data <- NULL
  
  ### need this at the end to check if append=TRUE can actually be done
  no.data <- is.null(data)
  
  ### check if data argument has been specified
  if (is.null(data)) {
    data <- sys.frame(sys.parent())
  } else {
    if (!is.data.frame(data))
      data <- data.frame(data)
  }
  
  #######################################################################################
  # Retrieve all data
  #######################################################################################
  mf <- match.call()
  
  mf.slab       <- mf[[match("slab",   names(mf))]]
  slab          <- eval(mf.slab,   data, enclos=sys.frame(sys.parent()))
  mf.OE         <- mf[[match("OE", names(mf))]]
  OE            <- eval(mf.OE, data, enclos=sys.frame(sys.parent()))
  mf.OE.se      <- mf[[match("OE.se", names(mf))]]
  OE.se         <- eval(mf.OE.se, data, enclos=sys.frame(sys.parent()))
  mf.OE.cilb    <- mf[[match("OE.cilb", names(mf))]]
  OE.cilb       <- eval(mf.OE.cilb, data, enclos=sys.frame(sys.parent()))
  mf.OE.ciub    <- mf[[match("OE.ciub", names(mf))]]
  OE.ciub       <- eval(mf.OE.ciub, data, enclos=sys.frame(sys.parent()))
  mf.OE.cilv    <- mf[[match("OE.cilv", names(mf))]]
  OE.cilv       <- eval(mf.OE.cilv, data, enclos=sys.frame(sys.parent()))
  mf.citl       <- mf[[match("citl", names(mf))]]
  citl          <- eval(mf.citl, data, enclos=sys.frame(sys.parent()))
  mf.citl.se    <- mf[[match("citl.se", names(mf))]]
  citl.se       <- eval(mf.citl.se, data, enclos=sys.frame(sys.parent()))
  mf.N          <- mf[[match("N", names(mf))]]
  N             <- eval(mf.N, data, enclos=sys.frame(sys.parent()))
  mf.O          <- mf[[match("O", names(mf))]]
  O             <- eval(mf.O, data, enclos=sys.frame(sys.parent()))
  mf.E          <- mf[[match("E", names(mf))]]
  E             <- eval(mf.E, data, enclos=sys.frame(sys.parent()))
  mf.Po         <- mf[[match("Po", names(mf))]]
  Po            <- eval(mf.Po, data, enclos=sys.frame(sys.parent()))
  mf.Pe         <- mf[[match("Pe", names(mf))]]
  Pe            <- eval(mf.Pe, data, enclos=sys.frame(sys.parent()))
  mf.Po.se      <- mf[[match("Po.se", names(mf))]]
  Po.se         <- eval(mf.Po.se, data, enclos=sys.frame(sys.parent()))
  
  #######################################################################################
  # Count number of studies
  #######################################################################################
  k <- 0
  
  if (!no.data) {
    k <- dim(data)[1]
  } else if (!is.null(OE)) {
    k <- length(OE)
  } else if (!is.null(OE.se)) {
    k <- length(OE.se)
  } else if (!is.null(OE.cilb)) {
    k <- length(OE.cilb)
  } else if (!is.null(OE.ciub)) {
    k <- length(OE.ciub)
  } else if (!is.null(OE.cilv)) {
    k <- length(OE.cilv)
  } else if (!is.null(citl)) {
    k <- length(citl)
  } else if (!is.null(citl.se)) {
    k <- length(citl.se)
  } else if (!is.null(N)) {
    k <- length(N)
  }  else if (!is.null(O)) {
    k <- length(O)
  } else if (!is.null(E)) {
    k <- length(E)
  } else if (!is.null(Po)) {
    k <- length(Po)
  } else if (!is.null(Po.se)) {
    k <- length(Po.se)
  } else if (!is.null(Pe)) {
    k <- length(Pe)
  }

  if (k<1) stop("No data provided!")
  
  if(is.null(OE))  OE <- rep(NA, times=k)
  if(is.null(OE.se)) OE.se <- rep(NA, times=k)
  if(is.null(OE.cilb)) OE.cilb <- rep(NA, times=k)
  if(is.null(OE.ciub)) OE.ciub <- rep(NA, times=k)
  if(is.null(OE.cilv)) OE.cilv <- rep(0.95, times=k) # Assume 95% CI by default 
  if(is.null(O)) O <- rep(NA, times=k)
  if(is.null(E)) E <- rep(NA, times=k)
  if(is.null(N)) N <- rep(NA, times=k)
  if(is.null(Po)) Po <- rep(NA, times=k)
  if(is.null(Pe)) Pe <- rep(NA, times=k)
  
  
  #######################################################################################
  # Assign study labels
  # taken from escalc
  #######################################################################################
  if (!is.null(slab)) {
    
    if (anyNA(slab))
      stop("NAs in study labels.")
    
    if (class(slab)=="factor") {
      slab <- as.character(slab)
    }
    
    ### check if study labels are unique; if not, make them unique
    if (anyDuplicated(slab))
      slab <- make.unique(slab)
    
    if (length(slab) != k)
      stop("Study labels not of same length as data.")
    
    ### add slab attribute to the cstat vector
    attr(OE, "slab") <- slab
  }
  

  #######################################################################################
  # Derive the OE ratio and its error variance
  # The order defines the preference for the final result
  #######################################################################################
  results <- list()
  results[[1]] <- data.frame(est=resoe.OE.se(OE=OE, OE.se=OE.se, g=g), method="OE and SE(OE)")
  results[[2]] <- data.frame(est=resoe.OE.ci(OE=OE, OE.cilb=OE.cilb, OE.ciub=OE.ciub, OE.cilv=OE.cilv, g=g), method="OE and CI(OE)")                                   #                                model=pars.default$model.oe)
  results[[3]] <- data.frame(est=resoe.O.E.N(O=O, E=E, N=N, correction = add, g=g), method="O, E and N")
  results[[4]] <- data.frame(est=resoe.O.Pe.N(O=O, Pe=Pe, N=N, correction = add, g=g), method="O, Pe and N")
  results[[5]] <- data.frame(est=resoe.E.Po.N(E=E, Po=Po, N=N, correction = add, g=g), method="E, Po and N")
  results[[6]] <- data.frame(est=resoe.Po.Pe.N(Po=Po, Pe=Pe, N=N, correction = add, g=g), method="Po, Pe and N")
  results[[7]] <- data.frame(est=resoe.O.Po.E(O=O, Po=Po, E=E, correction = add, g=g), method="O, E and Po") 
  results[[8]] <- data.frame(est=resoe.O.Pe.E(O=O, Pe=Pe, E=E, correction = add, g=g), method="O, E and Pe") 
  results[[9]] <- data.frame(est=resoe.O.E(O=O, E=E, correction = add, g=g), method="O and E")
  results[[10]] <- data.frame(est=resoe.Po.Pe(Po=Po, Pe=Pe, g=g), method = "Po and Pe")


  #t.citl   <- restore.oe.citl(citl=citl, citl.se=citl.se, O=O, Po=Po, N=N, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, 
  #                            model=pars.default$model.oe) 
  
  # Select appropriate estimate for 'theta' and record its source
  dat.est <- dat.se <- dat.method <-NULL
  for (i in 1:length(results)) {
    dat.est <- cbind(dat.est, (results[[i]])[,1]) 
    dat.se <- cbind(dat.se, sqrt(results[[i]][,2])) #take square root of error variance
    dat.method <- cbind(dat.method, as.character(results[[i]][,3]))
  }

  myfun = function(dat) { which.min(is.na(dat)) }
  sel.theta1 <- apply(dat.est, 1, myfun)
  sel.theta2 <- apply(dat.se, 1, myfun)
  sel.theta <- apply(cbind(sel.theta1, sel.theta2), 1, max) # only take the estimate for which we have an SE available
  
  theta <- dat.est[cbind(seq_along(sel.theta), sel.theta)] # Preferred estimate for theta 
  theta.se <- dat.se[cbind(seq_along(sel.theta), sel.theta)]# Preferred estimate for SE(theta)
  theta.source <-  dat.method[cbind(seq_along(sel.theta), sel.theta)] # Method used for estimating theta and its SE
  
  #######################################################################################
  # Derive confindence intervals
  #######################################################################################
  theta.cil <- theta.ciu <- rep(NA, k)
  
  # Directly transform the provided confidence limits for studies where the reported level is equal to the requested level
  if (is.null(g)) {
    theta.cil[OE.cilv==level] <- OE.cilb[OE.cilv==level]
    theta.ciu[OE.cilv==level] <- OE.ciub[OE.cilv==level]
  } else {
    theta.cil[OE.cilv==level] <- eval(parse(text=g), list(OE = OE.cilb[OE.cilv==level]))
    theta.ciu[OE.cilv==level] <- eval(parse(text=g), list(OE = OE.ciub[OE.cilv==level]))
  } 
  
  # Calculate the desired confidence intervals
  theta.cil[is.na(theta.cil)] <- (theta+qnorm((1-level)/2)*theta.se)[is.na(theta.cil)]
  theta.ciu[is.na(theta.ciu)] <- (theta+qnorm((1+level)/2)*theta.se)[is.na(theta.ciu)]
  
  #######################################################################################
  # Attempt to restore O, E and N
  #######################################################################################
  O[is.na(O)] <- (Po*N)[is.na(O)]
  E[is.na(E)] <- (Pe*N)[is.na(E)]
  N[is.na(N)] <- (O/Po)[is.na(N)]
  N[is.na(N)] <- (E/Pe)[is.na(N)]
  
  
  #######################################################################################
  # Sore results
  #######################################################################################
  ds <- data.frame(theta=theta, theta.se=theta.se, 
                   theta.cilb=theta.cil, theta.ciub=theta.ciu, 
                   theta.source=theta.source, O=O, E=E, N=N)
  
  if(is.null(slab) & !no.data) {
    slab <- rownames(data)
    rownames(ds) <- slab
  } else if (!is.null(slab)) {
    slab <- make.unique(as.character(slab))
    rownames(ds) <- slab
  }
  
  
  return(ds)
}

# Calculate OE and its error variance from O, E and N
resoe.O.E.N <- function(O, E, N, correction, g=NULL) {
  
  k <- length(O)
  out <- array(NA, dim=c(k,2))
  cc <- which(E==0)
  E[cc] <- E[cc]+correction
  O[cc] <- O[cc]+correction
  N[cc] <- N[cc]+correction
  out[,1] <- O/E 
  out[,2] <- ((O*(1-O/N))/(E**2)) # Error variance
  
  if(is.null(g)) {
    return (out)
  }

  toe <- toe.var <- rep(NA, k) # Transformed OE and its error variance

  for (i in 1:k) {
    oei <- out[i,1]
    toe[i] <- eval(parse(text=g), list(OE = oei))
    vi  <- out[i,2]
    names(oei) <- names(vi) <- "OE"
    toe.deriv <- deltaMethod(object=oei, g=g, vcov.=vi)
    toe.var[i] <- as.numeric(toe.deriv["SE"])**2
  }
  
  out <- cbind(toe, toe.var)
  return (out)
}

resoe.O.Pe.N <- function(O, Pe, N, correction, g=NULL) {
  return(resoe.O.E.N(O=O, E=Pe*N, N=N, correction=correction, g=g))
}

resoe.E.Po.N <- function(E, Po, N, correction, g=NULL) {
  return(resoe.O.E.N(O=Po*N, E=E, N=N, correction=correction, g=g))
}

resoe.Po.Pe.N <- function(Po, Pe, N, correction, g=NULL) {
  return(resoe.O.E.N(O=Po*N, E=Pe*N, N=N, correction=correction, g=g))
}

resoe.O.Po.E <- function(O, Po, E, correction, g=NULL) {
  return(resoe.O.E.N(O=O, E=E, N=O/Po, correction=correction, g=g))
}

resoe.O.Pe.E <- function(O, Pe, E, correction, g=NULL) {
  return(resoe.O.E.N(O=O, E=E, N=E/Pe, correction=correction, g=g))
}
  
resoe.OE.se <- function(OE, OE.se, g=NULL) {
  k <- length(OE)
  out <- array(NA, dim=c(k,2))
  out[,1] <- OE
  out[,2] <- OE.se**2
  
  if(is.null(g)) {
    return (out)
  }
  
  toe <- toe.var <- rep(NA, k) #Transformed OE and its error variance
  
  for (i in 1:k) {
    oei <- out[i,1]
    toe[i] <- eval(parse(text=g), list(OE = oei))
    vi  <- out[i,2]
    names(oei) <- names(vi) <- "OE"
    toe.var[i] <- as.numeric((deltaMethod(object=oei, g=g, vcov.=vi))["SE"])**2
  }
  
  out <- cbind(toe, toe.var)
  return (out)
}

resoe.OE.ci <- function(OE, OE.cilb, OE.ciub, OE.cilv, g=NULL) {
  k <- length(OE)
  out <- array(NA, dim=c(k,2))
  out[,1] <- OE
  out[,2] <- ((OE.ciub - OE.cilb)/(2*qnorm(0.5+OE.cilv/2)))**2 #Derive from 95% CI
  
  if(is.null(g)) {
    return (out)
  }
  
  toe <- toe.var <- rep(NA, k) # Transformed OE and its error variance
  
  for (i in 1:k) {
    toe[i] <- eval(parse(text=g), list(OE = out[i,1]))
    toe.cilb <- eval(parse(text=g), list(OE = OE.cilb[i]))
    toe.ciub <- eval(parse(text=g), list(OE = OE.ciub[i]))
    toe.var[i] <- ((toe.ciub - toe.cilb)/(2*qnorm(0.5+OE.cilv[i]/2)))**2 #Derive from 95% CI
  }
  
  out <- cbind(toe, toe.var)
  return (out)
}

resoe.O.E <- function(O, E, correction, g=NULL) {
  k <- length(O)
  out <- array(NA, dim=c(k,2))
  
  cc <- which(E==0)
  E[cc] <- E[cc]+correction
  O[cc] <- O[cc]+correction
  out[,1] <- O/E
  out[,2] <- (O/(E**2))
  
  if(is.null(g)) {
    return (out)
  }
  
  toe <- toe.var <- rep(NA, k) #Transformed OE and its error variance
  
  for (i in 1:k) {
    oei <- out[i,1]
    vi  <- out[i,2]
    names(oei) <- names(vi) <- "OE"
    toe[i] <- eval(parse(text=g), list(OE = oei))
    toe.var[i] <- as.numeric((deltaMethod(object=oei, g=g, vcov.=vi))["SE"])**2
  }
  
  out <- cbind(toe, toe.var)
  return (out)
}

resoe.citl <- function(citl, citl.se, Po, O, N, correction, g=NULL) {
  k <- length(citl)
  out <- array(NA, dim=c(k,2))
  
  cc <- which(O==0)
  O[cc] <- O[cc]+correction
  N[cc] <- N[cc]+correction
  Po[is.na(Po)] <- (O/N)[is.na(Po)]
  
  out[,1] <- -(exp(citl)*Po-exp(citl)-Po)
  
  for (i in 1:k) {
    citli <- citl[i]
    vi  <- citl.se[i]**2
    names(citli) <- names(vi) <- "citl"
    expr = paste("-(exp(citl)*",Po[i], "-exp(citl)-", Po[i],")")
    out[i,2] <- as.numeric((deltaMethod(object=citli, g=expr, vcov.=vi)["SE"]))**2
  }
  
  warning("Implementation not finalized!")
  
  if(is.null(g)) {
    return (out)
  }
  
  warning("Implementation still needed")
}

resoe.Po.Pe <- function (Po, Pe, g=NULL) {
  k <- length(Po)
  out <- array(NA, dim=c(k,2))
  
  out[,1] <- Po/Pe
  out[,2] <- NA # SE cannot be estimated from Po and Pe alone
  
  if(is.null(g)) {
    return (out)
  }
  
  toe <- toe.var <- rep(NA, k) #Transformed OE and its error variance
  
  for (i in 1:k) {
    oei <- out[i,1]
    vi  <- out[i,2]
    names(oei) <- names(vi) <- "OE"
    toe[i] <- eval(parse(text=g), list(OE = oei))
    toe.var[i] <- NA  # SE cannot be estimated from Po and Pe alone
  }
  
  out <- cbind(toe, toe.var)
  return (out)
}