# Meta-analysis is a statistical technique by which information from
# independent studies is assimilated. This function allows to perform
# fixed-effects and random-effects meta-analysis.
# r    : Vector of the effect sizes
# vars : Vector of the effect variances
###############################################################################
# Example
# example.r = c(0.10,0.30,0.35,0.65,0.45,0.15)
# example.var = c(0.03,0.03,0.05,0.01,0.05,0.02)
# uvmeta(example.r,example.var)
###############################################################################

#TODO: allow data transformations
uvmeta <- function(r, vars, model="random", method="MOM", labels, na.action,
                   pars, verbose=FALSE, ...) 
  UseMethod("uvmeta")

uvmeta.default <- function(r, vars, model="random", method="MOM", labels, na.action, 
                           pars, verbose=FALSE, ...)
{
  calcProfile <- function (mleObj, pars) {
    levels = pars$quantiles
    levels[which(pars$quantiles<0.5)]  = 1-(pars$quantiles[which(pars$quantiles<0.5)]*2)
    levels[which(pars$quantiles>=0.5)] = 1-(1-pars$quantiles[which(pars$quantiles>=0.5)])*2
    levels=unique(levels)
    pci = array(NA,dim=c(length(coef(mleObj)),length(pars$quantiles)))
    colnames(pci) = paste(pars$quantiles*100,"%",sep=" ")
    for (i in 1:length(levels)) {
      pcint <- confint(mle,level=levels[i], quietly=T)
      
      if (length(coef(mleObj))>1) {
        cols.select <- which(colnames(pcint) %in% colnames(pci))
        pci[,colnames(pcint)[cols.select]] <- pcint[,cols.select]
      } else {
        cols.select <- which(names(pcint) %in% colnames(pci))
        pci[,names(pcint)[cols.select]] <- pcint[cols.select]
      }
            
    }
    return(pci)
  }
  
  
  pars.default <- list(quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), 
                       n.chains=4, #JAGS (# chains)
                       n.adapt=5000, #JAGS
                       n.init=1000,  #JAGS
                       n.iter=10000) #JAGS
  
  if (length(r)!=length(vars)) {
    stop("The vectors 'r' and 'vars' have different lengths!")
  }
  
  ds <- as.data.frame(cbind(as.vector(r),as.vector(vars)))
  colnames(ds) <- c("theta","v")
  
  if (!missing(labels)) {
    if (length(labels) != length(r))
      stop("The vectors 'labels' and 'r' have different lengths!")
    rownames(ds) = labels
  } 
  
  if (missing(na.action)) 
    na.action <- "na.fail"
  if (length(na.action)) 
    ds <- do.call(na.action, list(ds))
  
  if (!missing(pars)) {
    for (i in 1:length(pars)) {
      element <- ls(pars)[i]
      pars.default[[element]] <- pars[[i]]
    }
  }

  est <- NA 
  
  
  #############################################################################
  # Start analyses
  #############################################################################
  numstudies = dim(ds)[1]
  dfr = numstudies-1
  
  if(numstudies < 3) {
    warning("There are very few primary studies!")
  }
  
  if (method == "MOM") { 
    results = as.data.frame(array(NA,dim=c(4, length(pars.default$quantiles)+2)))
    colnames(results) = c("Estimate","Var",paste(pars.default$quantiles*100,"%",sep=""))
    rownames(results) = c("mu","tausq","Q","Isq")
    
    # FIXED EFFECTS MODEL
    w = 1/ds$v
    
    #Combined effect
    weighted_Tbar = sum(ds$theta*w)/sum(w)
    
    # Variance of the combined effect
    var_T = 1/sum(w)
    
    # Standard error combined effect
    se_T = sqrt(var_T)
    
    # The Z-value
    z_T = weighted_Tbar/se_T
    
    # RANDOM EFFECTS MODEL
    Q = sum(w*(ds$theta-weighted_Tbar)**2)
    results["Q",] = c(Q,NA,rep(NA,length(pars.default$quantiles)))
    
    
    # Between-study variance
    if (model=="random" & Q > dfr) {
      re_C =  sum(w) - sum(w**2)/sum(w)
      between_study_var = (Q - dfr)/re_C
    } else {
      between_study_var = 0
    }

    # Within-study plus between-study variance
    re_v = vars + between_study_var
    
    # Updated weights
    re_w = 1/re_v
    
    # Combined effect
    re_weighted_Tbar =  sum(ds$theta*re_w)/sum(re_w)
    
    # Variance of the combined effect
    re_var_T  = 1/sum(re_w)
    
    # Standard error of combined effect
    re_se_T = sqrt(re_var_T)
    
    # The Z-value
    re_z_T = re_weighted_Tbar/re_se_T
    
    
    
    if (model=="random") {
      results["mu",] = c(re_weighted_Tbar,re_var_T,re_weighted_Tbar+qnorm(pars.default$quantiles)*sqrt(re_var_T))
      results["tausq",] = c(between_study_var,NA,rep(NA,length(pars.default$quantiles)))
      
      # Calculate I2 and its confidence limits
      Isq <- (results["Q",]-dfr)/results["Q",]
      Isq[which(Isq>1)] <- 1
      Isq[which(Isq<0)] <- 0
      results["Isq",] = Isq
    } else if (model=="fixed") {
      results["mu",] = c(weighted_Tbar,var_T,weighted_Tbar+qnorm(pars.default$quantiles)*sqrt(var_T))
      results["tausq",] = c(0,0,rep(NA,length(pars.default$quantiles)))
    }
    pred.int <- results["mu","Estimate"] + qt(pars.default$quantiles,df=(numstudies-2))*sqrt(results["tausq","Estimate"]+results["mu","Var"])
    names(pred.int) <- paste(pars.default$quantiles*100,"%",sep="")
    
    est <- list(results=results,model=model,df=dfr,numstudies=numstudies, pred.int=pred.int)
    
  } else if (method== "mle") {
    results = as.data.frame(array(NA,dim=c(4, length(pars.default$quantiles)+2)))
    colnames(results) = c("Estimate","Var",paste(pars.default$quantiles*100,"%",sep=""))
    rownames(results) = c("mu","tausq","Q","Isq")
    
    #mle.loglik <- function( theta, tausq, ds) {
    #  loglik <- -0.5*sum(log(2*pi*(ds$v+tausq)))-0.5*sum(((ds$theta-theta)**2)/(ds$v+tausq))
    #  return (-loglik) #return negative log-likelihood
    #}
    mle.loglik.random <- function(theta, tausq, ds) { #random effects
          loglik <- sum(dnorm(x=ds$theta,mean=theta, sd=sqrt(tausq+ds$v),log=T))
          return (-loglik)
    }
    mle.loglik.fixed <- function(theta, ds) { #fixed effects
      loglik <- sum(dnorm(x=ds$theta,mean=theta, sd=sqrt(ds$v),log=T))
      return (-loglik)
    }
    
    if (model=="random") {
      mle <- mle2(minuslogl=mle.loglik.random,start=list(theta=0, tausq=0),data=list(ds=ds),method="L-BFGS-B",lower=list(theta=-Inf,tausq=0))
      
      profile = calcProfile(mle,pars.default)  #Use profile likelihood confidence intervals
      results["mu",] = c(coef(mle)["theta"],diag(vcov(mle))["theta"],profile[1,])
      results["tausq",] = c(coef(mle)["tausq"],diag(vcov(mle))["tausq"],profile[2,])
      
      w <- 1/(ds$v+coef(mle)["tausq"])
      Q <- sum(w*(ds$theta-coef(mle)["theta"])**2)
      results["Q",] = c(Q,NA,rep(NA,length(pars.default$quantiles)))
      
      # Calculate I2 
      Isq <- (results["Q",]-dfr)/results["Q",]
      Isq[which(Isq>1)] <- 1
      Isq[which(Isq<0)] <- 0
      results["Isq",] = Isq
    } else {
      mle <- mle2(minuslogl=mle.loglik.fixed,start=list(theta=0),data=list(ds=ds))
      profile = calcProfile(mle,pars.default)  #Use profile likelihood confidence intervals
      results["mu",] = c(coef(mle)["theta"],diag(vcov(mle))["theta"],profile)
      results["tausq",] = c(0,0,rep(0,length(pars.default$quantiles)))
    }
    
    pred.int <- results["mu","Estimate"] + qt(pars.default$quantiles,df=(numstudies-2))*sqrt(results["tausq","Estimate"]+results["mu","Var"])
    names(pred.int) <- paste(pars.default$quantiles*100,"%",sep="")
    
    est <- list(results=results,model=model,df=dfr,numstudies=numstudies, pred.int=pred.int, loglik=-attr(mle,"min"))
  } else if (method == "bayes") { 
    quiet = !verbose
    
    modelfile <-  if (model=="random") system.file(package="metamisc", "model", "uvmeta_ranef.bug") else system.file(package="metamisc", "model", "uvmeta_fixef.bug")
    jags <- jags.model(modelfile,
                       data = list('r' = ds$theta,
                                   'vars' = ds$v,
                                   'k' = numstudies), #prior precision matrix
                       n.chains = pars.default$n.chains,
                       n.adapt = pars.default$n.adapt,
                       quiet = quiet)
    update(jags, pars.default$n.init) #initialize
    samples <- coda.samples(jags, c('mu','tausq','Q','Isq','theta.new'),n.iter=pars.default$n.iter)
    
    results <- summary(samples,quantiles=pars.default$quantiles) 
    pred.int=(results[[2]])["theta.new",]
    
    results.overview = as.data.frame(array(NA,dim=c(dim(results[[1]])[1], length(pars.default$quantiles)+2)))
    colnames(results.overview) = c("Estimate","Var",paste(pars.default$quantiles*100,"%",sep=""))
    rownames(results.overview) = rownames(results[[2]])
    results.overview[,1] = (results[[1]])[,"Mean"]
    results.overview[,2] = (results[[1]])[,"SD"]**2
    for (i in 1:length(pars.default$quantiles)) {
      results.overview[,(i+2)] = (results[[2]])[,i]
    }
    results.overview = results.overview[c("mu","tausq","Q","Isq"),]
    
    est <- list(results=results.overview,model=model,df=dfr,numstudies=numstudies,pred.int=pred.int)
  } else {
    stop("Invalid meta-analysis method!")
  }
  # } else if (method=="pl") {
  #   results = as.data.frame(array(NA,dim=c(4, length(pars$quantiles)+2)))
  #  colnames(results) = c("Estimate","Var",paste(pars$quantiles*100,"%",sep=""))
  #  rownames(results) = c("mu","tausq","Q","Isq")
  #  
  #  mle.loglik <- function(theta, tausq, ds) {
  #    loglik <- sum(dnorm(x=ds$theta,mean=theta, sd=sqrt(tausq+ds$v),log=T))
  #    return (-loglik)
  #  }
  #      
  #  #### 5.99 ===> qchisq (0.95,df=2)
  #  if (model == "random") 
  #  {
  #    mle <- mle2(minuslogl=mle.loglik, start=list(theta=0,tausq=0), data=list(ds=ds),method="L-BFGS-B",lower=list(theta=-Inf,tausq=0))
  #    mle.cov <- vcov(mle)
  #    p0 <- profile(mle)
  #    
  #    levels = pars$quantiles
  #    levels[which(pars$quantiles<0.5)]  = 1-(pars$quantiles[which(pars$quantiles<0.5)]*2)
  #    levels[which(pars$quantiles>=0.5)] = 1-(1-pars$quantiles[which(pars$quantiles>=0.5)])*2
  #    pci = array(NA,dim=c(2,length(levels)))
  #    colnames(pci) = paste(pars$quantiles*100,"%",sep=" ")
  #    
  #    for (i in 1:length(levels)) {
  #      pcint <- confint(p0,level=levels[i])
  #      cols.select <- which(colnames(pcint) %in% colnames(pci))
  #      pci[,colnames(pcint)[cols.select]] <- pcint[,cols.select]
  #    }
  #    
  #    wt <- 1/(coef(mle)["tausq"]+ds$v)
  #    Q <- sum(wt*(ds$theta-coef(mle)["theta"])**2)
  #    results["Q",] = c(Q,NA,qchisq(pars$quantiles,df=dfr))
  #    
  #    
  #    # Use profile log-likelihood to calculate confidence intervals
  #    results["mu",] = c(coef(mle)["theta"],mle.cov[1,1],pci[1,])
  #    results["tausq",] = c(mle.tausq,mle.cov[2,2],pci[2,])
  #    
  #    # Calculate I2 and its confidence limits
  #    Isq <- (results["Q",]-dfr)/results["Q",]
  #    Isq[which(Isq>1)] <- 1
  #    Isq[which(Isq<0)] <- 0
  #    results["Isq",] = Isq
  # }
  #  pred.int <- results["mu","Estimate"] + qt(pars$quantiles,df=(numstudies-2))*sqrt(results["tausq","Estimate"]+results["mu","Var"])
  #  names(pred.int) <- paste(pars$quantiles*100,"%",sep="")

  #  est <- list(results=results,model=model,df=dfr,numstudies=numstudies, pred.int=pred.int)
 
  attr(est$results,"quantiles") = pars.default$quantiles
  est$data <- ds
  est$na.action <- na.action
  est$method <- method
  est$call <- match.call()
  class(est) <- "uvmeta"
  return(est)
}

plot.uvmeta <- function(x, ...) {
  
  quantiles <- attr(x$results,"quantiles")
  ci <- x$data[,"theta"]+t(qnorm(quantiles)*matrix(rep(sqrt(x$data[,"v"]),length(quantiles)),nrow=(length(quantiles)), ncol=dim(x$data)[1],byrow=T))
  
  xlim <- c(min(ci),max(ci))
  ylim <- c(2,(x$numstudies+5))

  par(mfrow=c(1,1), mar=( c(5, 12, 4, 4) + 0.1))
  
  loc = c((x$numstudies+4):3)
    
  lcol = "#EBEBEB"
  plot(-500,-500,xlim=xlim,ylim=ylim,xlab="", ylab="",yaxt="n",...)
  axis(2,at=c((x$numstudies+4):5,3),labels=c(rownames(x$data),"Pooled estimate"),las=1)
  
  abline(v=0.00,lty=2,col=lcol)
  for (i in 1:x$numstudies) {
    yloc = loc[i]
    points(x$data[i,"theta"],yloc,pch=15)
    lines(c(min(ci[i,]),max(ci[i,])),c(yloc,yloc))
    for (j in 1:dim(ci)[2])
      lines(c(ci[i,j],ci[i,j]),c((yloc-0.1),(yloc+0.1)),pch=3)
  }
  
  ci.bounds <- x$results["mu",-match(c("Estimate","Var"),colnames(x$results))]
  for (i in 1:length(ci.bounds))
   lines(c(ci.bounds[i],ci.bounds[i]),c(2.9,3.1),pch=3)
  
  lines(c(min(ci.bounds),max(ci.bounds)),c(3,3))
  points(x$results["mu","Estimate"],3,pch=23,bg="white")
  
  box()
}



print.uvmeta <- function(x, ...)
{
  out <- (x$results)
  text.model <- if (x$model=="fixed") "Fixed" else "Random"
  text.method <- if(x$method=="bayes") "credibility" else "confidence"
  cat(paste(text.model,"effects estimates with corresponding", text.method, "intervals:\n\n"))
	print(out)
  if (x$model=="random") {
    cat(paste("\n\nPrediction interval for mu:\n\n"))
    print(x$pred.int)
  }
  if(x$method=="mle") { #display MLE
    cat(paste("\nLog-likelihood: ", round(x$loglik,2),"\n"))
  }
  
	out
}


summary.uvmeta <- function(object, ...)
{
    cat("Call:\n")
    print(object$call)
    if (object$model=="fixed")  cat(paste("\nFixed effects summary:\t",round(object$results["mu","Estimate"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
    if (object$model=="random") {
        cat(paste("\nRandom effects summary:\t",round(object$results["mu","Estimate"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
        cat(paste("\n\nTau squared: \t\t",round(object$results["tausq","Estimate"],5),sep=""))
    }
    Q_p = 1-pchisq(object$results["Q","Estimate"],df=object$df)
    cat(paste("\nCochran's Q statistic: \t",round(object$results["Q","Estimate"],5)," (p-value: ",round(Q_p,5),")",sep=""))
    cat(paste("\nI-square index: \t", round(object$results["Isq","Estimate"]*100,3)," %\n",sep=""))
}



