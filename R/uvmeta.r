# Meta-analysis is a statistical technique by which information from
# independent studies is assimilated. This function allows to perform
# fixed-effects and random-effects meta-analysis.
# r    : Vector of the effect sizes
# vars : Vector of the effect variances
################################################################################
# Author  : Thomas Debray
# Version : 10 May 2011
################################################################################
# Example
# example.r = c(0.10,0.30,0.35,0.65,0.45,0.15)
# example.var = c(0.03,0.03,0.05,0.01,0.05,0.02)
# macc(example.r,example.var)
################################################################################
uvmeta <- function(r, vars, model="random", method="MOM", pars=list(quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), n.chains=4, n.adapt=5000, n.init=1000, n.iter=10000), ...) UseMethod("uvmeta")

uvmetaMOM <- function(r,vars, model="random", pars=list(quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975))) {
    # Degrees of freedom
    numstudies = length(r)
    dfr = numstudies-1
    
    if(numstudies < 3){warning("There are very few primary studies!")}

    ############################################################################
    # FIXED EFFECTS MODEL
    ############################################################################
    w = 1/vars

    #Combined effect
    weighted_Tbar = sum(r*w)/sum(w)

    # Variance of the combined effect
    var_T = 1/sum(w)

    # Standard error combined effect
    se_T = sqrt(var_T)

    # The Z-value
    z_T = weighted_Tbar/se_T

    ############################################################################
    # RANDOM EFFECTS MODEL
    ############################################################################
    # Q-statistic
    Q = sum(w*(r-weighted_Tbar)**2)
    I_sq = 0
    
    # Higgins and Thompson (2002) have also developed a confidence interval for
    # I2. The interval is formulated by calculating another of their proposed
    # measures of heterogeneity, the H2 index obtained by
    # (Higgins & Thompson, 2002, p. 1545, eq. 6) also known as Birges
    # ratio (Birge, 1932)
    H_sq = 1    #H2
    se_lnH = 0  #SE(ln(H2))

    # Between-study variance
    if (Q > dfr) {
        re_C =  sum(w) - sum(w**2)/sum(w)
        between_study_var = (Q - dfr)/re_C
        I_sq = (Q-dfr)/Q
        H_sq = Q/dfr
        se_lnH = (log(Q)-log(dfr))/(2*(sqrt(2*Q)-sqrt((2*length(r))-3)))
    } else {
        between_study_var = 0
        se_lnH = sqrt((1/(2*(length(r)-2)))*(1-(1/(3*((length(r)-2)**2)))))
    }
    #varQ = 2*dfr + 4*(sum(w)-sum(w**2)/sum(w))*between_study_var + 2*(sum(w**2)-2*((sum(w**3)/sum(w))+(sum(w**2)**2)/sum(w)**2))*between_study_var**2
    #varTauSq = varQ/(sum(w)-sum(w**2)/sum(w))**2

    # Within-study plus between-study variance
    re_v = vars + between_study_var

    # Updated weights
    re_w = 1/re_v

    # Combined effect
    re_weighted_Tbar =  sum(r*re_w)/sum(re_w)

    # Variance of the combined effect
    re_var_T  = 1/sum(re_w)

    # Standard error of combined effect
    re_se_T = sqrt(re_var_T)

    # The Z-value
    re_z_T = re_weighted_Tbar/re_se_T
    
    fixef.results = list(mean=weighted_Tbar,var=var_T)
    ranef.results = list(mean=re_weighted_Tbar,var=re_var_T,tauSq=between_study_var)
    H2.results    = list(H2=H_sq,se_lnH=se_lnH)
    I2.results    = list(I2=I_sq)
    
    ############################################################################
    # Q statistics
    ############################################################################
    # The Q statistic has a chi-square distribution with k - 1 degrees of
    # freedom, k being the number of studies. Thus, Q values higher than the
    # critical point for a given significance level alfa enable us to reject the
    # null hypothesis and conclude that there is statistically significant
    # between-study variation
    #Q.critical = qchisq(0.95,df=(length(r)-1))
    Q_p = 1-pchisq(Q,df=(length(r)-1))
    Q.results = list(Q=Q,p.value=Q_p)

    results = array(NA,dim=c(4, length(pars$quantiles)+2))
    colnames(results) = c("Mean","Var",paste(pars$quantiles*100,"%",sep=""))
    rownames(results) = c("mu","tausq","Q","Isq")
    results[3,] = c(Q,NA,qchisq(pars$quantiles,df=dfr))
    
    if (model=="random") {
        results[1,] = c(re_weighted_Tbar,re_var_T,re_weighted_Tbar+qnorm(pars$quantiles)*sqrt(re_var_T))
	results[2,] = c(between_study_var,NA,rep(NA,length(pars$quantiles)))
	results[4,] = c(I_sq,NA,rep(NA,length(pars$quantiles)))
    } else if (model=="fixed") {
        results[1,] = c(weighted_Tbar,var_T,weighted_Tbar+qnorm(pars$quantiles)*sqrt(var_T))
	results[2,] = c(0,0,rep(NA,length(pars$quantiles)))
    }
 
    ############################################################################
    # Output
    ############################################################################
    out <- list(results=results,model=model,df=dfr,numstudies=numstudies)
    class(out) <- "uvmeta"
    return (out)
}

uvmetaBayes <- function(r,vars, model="random",pars=list(quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), n.chains=4, n.adapt=5000, n.init=1000, n.iter=10000)) {
	numstudies = length(r)	
	dfr = numstudies-1

	modelfile <-  if (model=="random") system.file(package="metamisc", "model", "uvmeta_ranef.bug") else system.file(package="metamisc", "model", "uvmeta_fixef.bug")
	jags <- jags.model(modelfile,
                     data = list('r' = r,
                                 'vars' = vars,
                                 'k' = numstudies), #prior precision matrix
                     n.chains = pars$n.chains,
                     n.adapt = pars$n.adapt)
	update(jags, pars$n.init) #initialize
	samples <- coda.samples(jags, c('mu','tausq','Q','Isq'),n.iter=pars$n.iter)

	results <- summary(samples,quantiles=pars$quantiles)

	results.overview = array(NA,dim=c(dim(results[[1]])[1], length(pars$quantiles)+2))
	colnames(results.overview) = c("Mean","Var",paste(pars$quantiles*100,"%",sep=""))
	rownames(results.overview) = rownames(results[[2]])
	results.overview[,1] = (results[[1]])[,"Mean"]
	results.overview[,2] = (results[[1]])[,"SD"]**2
	for (i in 1:length(pars$quantiles)) {
		results.overview[,(i+2)] = (results[[2]])[,i]
	}
	
	out <- list(results=results.overview,model=model,df=dfr,numstudies=numstudies)
    	class(out) <- "uvmeta"
    	return (out)
}

uvmeta.default <- function(r,vars, model="random", method="MOM", ...)
{
    x <- as.vector(r)
    y <- as.vector(vars)
    est <- NA    
    if (length(x)!=length(y)) {warning("The vectors 'r' and 'vars' have a different size!")}
    if (method == "MOM") { est <- uvmetaMOM(x, y, model) }
    else if (method == "bayes") { est <- uvmetaBayes(x,y, model) }

    est$call <- match.call()
    class(est) <- "uvmeta"
    est
}

print.uvmeta <- function(x, ...)
{
	out <- (x$results)
	print(out)
	out
}


predict.uvmeta <- function(object, level = 0.95, ...)
{
  alpha = (1-level)/2

  #The correct number of degrees of freedom for this t distribution is complex, and we use a value of k–2 largely for pragmatic reasons. (Riley 2011)
  df = 2 
  
  pred.mean  <- object$results["mu","Mean"]
  pred.lower <- object$results["mu","Mean"] + qt(alpha,df=(object$numstudies-df))*sqrt(object$results["tausq","Mean"]+object$results["mu","Var"])
  pred.upper <- object$results["mu","Mean"] + qt((1-alpha),df=(object$numstudies-df))*sqrt(object$results["tausq","Mean"]+object$results["mu","Var"])
  predint <- c(pred.mean,pred.lower,pred.upper)
  names(predint) <- c("Estimate", paste((alpha*100),"%"),paste(((1-alpha)*100),"%"))
  predint
}

summary.uvmeta <- function(object, ...)
{
    cat("Call:\n")
    print(object$call)
    if (object$model=="fixed")  cat(paste("\nFixed effects summary:\t",round(object$results["mu","Mean"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
    if (object$model=="random") {
        cat(paste("\nRandom effects summary:\t",round(object$results["mu","Mean"],5))," (SE: ",round(sqrt(object$results["mu","Var"]),5), ")",sep="")
        cat(paste("\n\nTau squared: \t\t",round(object$results["tausq","Mean"],5),sep=""))
    }
    Q_p = 1-pchisq(object$results["Q","Mean"],df=object$df)
    cat(paste("\nCochran's Q statistic: \t",round(object$results["Q","Mean"],5)," (p-value: ",round(Q_p,5),")",sep=""))
    cat(paste("\nI-square index: \t", round(object$results["Isq","Mean"]*100,3)," %\n",sep=""))
}



