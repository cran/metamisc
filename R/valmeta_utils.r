# Adapted from R package 'car' to avoid package loading issues in R-forge
deltaMethod <- function (object, g, vcov., func = g, constants, level=0.95, ...) {
  if (!is.character(g)) 
    stop("The argument 'g' must be a character string")
  
  para <- object         
  para.names <- names(para)
  g <- parse(text = g)
  q <- length(para)
  for (i in 1:q) {
    assign(names(para)[i], para[i])
  }
  if(!missing(constants)){
    for (i in seq_along(constants)) assign(names(constants[i]), constants[[i]])}
  est <- eval(g)
  names(est) <- NULL
  gd <- rep(0, q)
  for (i in 1:q) {
    gd[i] <- eval(D(g, names(para)[i]))
  }
  se.est <- as.vector(sqrt(t(gd) %*% vcov. %*% gd))
  result <- data.frame(Estimate = est, SE = se.est, row.names = c(func))
  result
}



restore.c.var.hanley <- function(cstat, N.subjects, N.events, restore.method=4, g=NULL) {
  
  fHanley <- function(cstat, nstar, mstar, m, n) {
    ((cstat*(1-cstat)*(1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n))
  }
  
  n <- N.events #Number of events
  m <- N.subjects-N.events #Number of non-events
  
  if (restore.method==2) {
    mstar <- m-1
    nstar <- n-1
  } else if (restore.method==4) {
    mstar <- nstar <- N.subjects/2-1
  } else {
    stop ("Method not implemented yet!")
  }
  
  # Crude estimate
  cstat.var <- fHanley(cstat=cstat, m=m, n=n, mstar=mstar, nstar=nstar)
  if (is.null(g)) {
    return(cstat.var)
  }
  
  # Apply delta method if a transformation for c is defined
  ti <- rep(NA, length(cstat))
  for (i in 1:length(cstat)) {
    ci <- cstat[i]
    vi <- cstat.var[i]
    names(ci) <- names(vi) <- "cstat"
    ti[i] <- as.numeric((deltaMethod(object=ci, g=g, vcov.=vi))["SE"])**2
  }
  return(ti)
}


restore.c.var.se <- function(cstat, c.se, g=NULL) {
  if(is.null(g)) {
    return (c.se**2)
  }
  
  ti <- rep(NA, length(cstat))
  
  for (i in 1:length(cstat)) {
    ci <- cstat[i]
    vi <- c.se[i]**2
    names(ci) <- names(vi) <- "cstat"
    ti[i] <- as.numeric((deltaMethod(object=ci, g=g, vcov.=vi))["SE"])**2
  }
  return(ti)
}

restore.c.var.hanley2<- function(sd.LP, N.subjects, N.events, restore.method=4, g=NULL) {
  cstat <- calculate.cstat.sdPI(sd.LP, g=g)
  return(restore.c.var.hanley(cstat=cstat, N.subjects=N.subjects, N.events=N.events, restore.method=restore.method, g=g))
}

restore.c.var.se <- function(cstat, c.se, g=NULL) {
  if(is.null(g)) {
    return (c.se**2)
  }
  
  ti <- rep(NA, length(cstat))
  
  for (i in 1:length(cstat)) {
    ci <- cstat[i]
    vi <- c.se[i]**2
    names(ci) <- names(vi) <- "cstat"
    ti[i] <- as.numeric((deltaMethod(object=ci, g=g, vcov.=vi))["SE"])**2
  }
  return(ti)
}

restore.c.var.ci <- function(cil, ciu, level, g=NULL) {
  if(!is.null(g)) {
    lower <- eval(parse(text=g), list(cstat = cil))
    upper <- eval(parse(text=g), list(cstat = ciu))
  } else {
    lower <- cil
    upper <- ciu
  }
  if(missing(level)) level <- rep(0.95, length(cil))
  
  return(((upper - lower)/(2*qnorm((1-level)/2)))**2)
}

calculate.cstat.theta <- function(cstat, g=NULL) {
  if(is.null(g)) {
    return(cstat)
  }
  return(eval(parse(text=g), list(cstat = cstat)))
}

calculate.cstat.sdPI <- function (sdPI, g=NULL) {
  myfun <- function(x, sd.lp) {
    inv.logit(sqrt(2)*sd.lp*x)*dnorm(x, mean=0, sd=1)
  }
  cstat <- rep(NA, length(sdPI))
  for (i in 1:length(sdPI)) {
    cstat[i] <- ifelse(is.na(sdPI[i]), NA, 2*(integrate(myfun, lower=0, upper=+Inf, sd.lp=sdPI[i]))$value)
  }
  if(is.null(g)) {
    return(cstat)
  }
  return(eval(parse(text=g), list(cstat = cstat)))
}

# Calculate OE and its SE from O, E and N
restore.oe.OE <- function(OE, OE.se, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  out <- array(NA, dim=c(length(OE),2))
  
  if (model == "normal/identity") {
    out[,1] <- OE
    out[,2] <- OE.se
  } else if (model %in% c("normal/log", "poisson/log")) {
    out[,1] <- log(OE)
    out[,2] <- OE.se/OE
  } 
  
  # Extrapolation not possible
  if (!is.na(t.ma) & class(t.val)=="numeric") {
    out[which(t.val!=t.ma),] <- NA
  }
  
  return (out)
}

restore.oe.OE.95CI <- function(OE, OE.95CI, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  out <- array(NA, dim=c(length(OE),2))
  
  if (model == "normal/identity") {
    OE.se <- abs((OE.95CI[,2] - OE.95CI[,1])/(2*qnorm(0.975))) #Derive from 95% CI
    return (restore.oe.OE(OE=OE, OE.se=OE.se, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
  } else if (model %in% c("normal/log", "poisson/log")) {
    out[,1] <- log(OE)
    out[,2] <- abs((log(OE.95CI[,2]) - log(OE.95CI[,1]))/(2*qnorm(0.975)))
  }
  
  # Extrapolation not possible
  if (!is.na(t.ma) & class(t.val)=="numeric") {
    out[which(t.val!=t.ma),] <- NA
  }
  return (out)
}

# Calculate OE and its SE from O, E and N
restore.oe.O.E.N <- function(O, E, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  out <- array(NA, dim=c(length(O),2))
  
  if (model == "normal/identity") {
    cc <- which(E==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    N[cc] <- N[cc]+correction
    out[,1] <- O/E 
    out[,2] <- sqrt((0*(1-O/N))/(E**2))
  } else if (model %in% c("normal/log", "poisson/log")) {
    cc <- which(E==0 | O==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    N[cc] <- N[cc]+correction
    out[,1] <- (log(O)-log(E))
    out[,2] <- sqrt((1-(O/N))/O)
  } 
  
  # Apply extrapolation or omit study where t.val != t.ma
  if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
    out[which(t.val!=t.ma),] <- NA
  } else if (t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
    Po.new <- 1-exp(t.ma*log(1-(O/N))/t.val)
    Pe.new <- 1-exp(t.ma*log(1-(E/N))/t.val)
    
    if (model == "normal/identity") {
      theta.new <- (Po.new/Pe.new)
      theta.se.new <- sqrt(((t.ma**2)*exp(2*t.ma*log(1-(O/N))/t.val)*(O/N))/((t.val**2)*N*(1-(O/N))*(1-exp(t.ma*log(1-(E/N))/t.val))**2))
    } else if (model %in% c("normal/log", "poisson/log")) {
      theta.new <- log(Po.new/Pe.new)
      theta.se.new <- sqrt(((t.ma**2)*(O/N)*exp(2*t.ma*log(1-(O/N))/t.val))/((t.val**2)*N*(1-(O/N))*((1-exp(t.ma*log(1-(O/N))/t.val))**2)))
    } else {
      stop("Model not supported")
    }
    out[which(t.val!=t.ma),1] <- theta.new[which(t.val!=t.ma)]
    out[which(t.val!=t.ma),2] <- theta.se.new[which(t.val!=t.ma)]
  }
  
  #Don't provide O/E for studies where N is missing
  out[is.na(N),] <- NA
  
  return (out)
}

# Calculate OE and its SE from O and E
restore.oe.O.E <- function(O, E, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  out <- array(NA, dim=c(length(O),2))
  
  if (model == "normal/identity") {
    cc <- which(E==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    out[,1] <- O/E
    out[,2] <- sqrt(0/(E**2))
  } else if (model %in% c("normal/log", "poisson/log")) {
    cc <- which(E==0 | O==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    out[,1] <- (log(O)-log(E))
    out[,2] <- sqrt(1/O)
  } 
  
  # Extrapolation not possible
  if (!is.na(t.ma) & class(t.val)=="numeric") {
    out[which(t.val!=t.ma),] <- NA
  } 
  
  return (out)
}

restore.oe.PoPe <- function (Po, Pe, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  out <- array(NA, dim=c(length(Po),2))
  
  if (missing(t.val)) {
    t.val <- rep(NA, length(Po))
  }
  
  if (model == "normal/identity") {
    
    out[,1] <- Po/Pe
    # Extrapolation not possible
    if (!is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma),] <- NA
    }
  } else if (model %in% c("normal/log", "poisson/log")) {
    out[,1] <- log(Po)-log(Pe)
    
    if (!is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma),] <- NA
    }
  } 
  
  return (out)
}

restore.oe.OPoE <- function(O, Po, E, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.O.E.N(O=O, E=E, N=O/Po, correction=correction, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}

restore.oe.OPeE <- function(O, Pe, E, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.O.E.N(O=O, E=E, N=E/Pe, correction=correction, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}

restore.oe.OPeN <- function(O, Pe, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.O.E.N(O=O, E=Pe*N, N=N, correction=correction, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}

restore.oe.EPoN <- function(E, Po, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.O.E.N(O=Po*N, E=E, N=N, correction=correction, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}

restore.oe.PoPeN <- function (Po, Pe, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.O.E.N(O=Po*N, E=Pe*N, N=N, correction=correction, t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}

# Restore OE ratio from calibration-in-the-large
restore.oe.citl <- function(citl, citl.se, O, Po, N, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  out <- array(NA, dim=c(length(citl),2))
  
  if (missing(t.val)) {
    t.val <- rep(NA, length(citl))
  }
  
  Po[is.na(Po)] <- (O/N)[is.na(Po)]
  
  # Apply extrapolation using Poisson distribution
  if (t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
    Po[which(t.val!=t.ma)] <- (1-exp(t.ma*log(1-Po)/t.val))[which(t.val!=t.ma)]
  } else if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
    Po[which(t.val!=t.ma)] <- NA # Omit studies where follow-up duration is improprer
  }
  
  if (model == "normal/identity") {
    out[,1] <- -(exp(citl)*(Po)-exp(citl)-(Po))

  }
  if (model %in% c("normal/log", "poisson/log")) {
    out[,1] <- (log(-(exp(citl)*(Po)-exp(citl)-(Po))))
    out[,2] <- sqrt((((Po-1)**2)*((Po**2)+1)*(exp(Po+citl))**2)*(citl.se**2)/((-Po*exp(citl))+Po+exp(citl))**2)
  }
  return (out)
}

restore.oe.var.seOE1 <- function(se, OE, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  if (model == "normal/identity") {
    out <- se**2
    
    if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  } else if (model %in% c("normal/log", "poisson/log")) {
    out <- (se**2)/(OE**2) # Equation 16 in appendix BMJ paper
    
    if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  }
  return(out)
}

restore.oe.var.seOE2 <- function(se, O, E, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  return(restore.oe.var.seOE1(se=se, OE=(O/E), t.extrapolate=t.extrapolate, t.ma=t.ma, t.val=t.val, model=model))
}
  
#In most situations, O and E are reported separately without any estimate of uncertainty.
#In the following derivations, we regard E as a fixed constant. We treat O as a binomially distributed
#variable since O is given as the number of successes (events) from N subjects
restore.oe.var.OEN <- function(O, E, N, correction = 0.5, t.extrapolate=F, t.ma=NA, t.val, model="normal/log") {
  if (model == "normal/identity") {
    cc <- which(E==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    N[cc] <- N[cc]+correction
    Po <- O/N
    
    out <- O*(1-Po)/(E**2)
    
    if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  } else if (model %in% c("normal/log", "poisson/log")) {
    cc <- which(E==0 | O==0)
    E[cc] <- E[cc]+correction
    O[cc] <- O[cc]+correction
    N[cc] <- N[cc]+correction
    Po <- O/N
    
    out <- (1-Po)/(O)
    
    if (!t.extrapolate & !is.na(t.ma) & class(t.val)=="numeric") {
      out[which(t.val!=t.ma)] <- NA
    }
  }
  return(out)
}



# See params of geom_smooth for more details
# Smoothed calibration plot: use  formula = obsy ~ splines::bs(predy, 3)
plotCalibration <- function(predy, obsy, modelname="Model", 
                            formula = obsy ~ predy, 
                            method="glm", se = se, 
                            level=0.95, 
                            fam=binomial, ...) {
  #require(ggplot2)
  #require(ggExtra) ## Do the extra outside the package
  
  # Use formula instead
  
  # Add density plot under gg plot
  #predy<-rnorm(300)
  #obsy<-rt(300,df=10)
  #family <- gaussian
  
  xy <- data.frame(predy, obsy)
  
  scatter <- ggplot(xy, aes(x=predy, y=obsy)) +
    labs(x = "Predicted", y="Observed") + 
    scale_x_continuous(limits=c(min(predy),max(obsy))) + 
    scale_y_continuous(limits=c(min(predy),max(obsy))) #+ 
  
  
  # Add reference line for perfect calibration
  scatter <- scatter + geom_abline(aes(slope=1, intercept=0, linetype="Perfect calibration"), size=1)
  scatter <- scatter + geom_smooth(aes(linetype=modelname), method = method, 
                                   se = se, level=level, method.args = list(family = fam), ...)

  
  #scatter <- scatter + ggMarginal(data = xy, x = "predy", y = "obsy", margins = "x", type="histogram", size=4)
  #scatter <- scatter + labs(x = "Predicted", y="Observed") 
  scatter
  
  
}
                                                                                                                               