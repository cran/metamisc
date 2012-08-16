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
uvmeta <- function(r, vars, method="MOM", ...) UseMethod("uvmeta")

uvmetaMOM <- function(r,vars) {
    # Degrees of freedom
    numstudies = length(r)
    dfr = numstudies-1

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
    varQ = 2*dfr + 4*(sum(w)-sum(w**2)/sum(w))*between_study_var + 2*(sum(w**2)-2*((sum(w**3)/sum(w))+(sum(w**2)**2)/sum(w)**2))*between_study_var**2
    varTauSq = varQ/(sum(w)-sum(w**2)/sum(w))**2

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
    ranef.results = list(mean=re_weighted_Tbar,var=re_var_T,tauSq=between_study_var,varTauSq=varTauSq)
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
    Q.results = list(Q=Q,var=varQ,p.value=Q_p)
    
    ############################################################################
    # Output
    ############################################################################
    out <- list(fixef=fixef.results,ranef=ranef.results,Q=Q.results,H2=H2.results,I2=I2.results,df=dfr,numstudies=numstudies)
    return (out)
}

uvmeta.default <- function(r,vars, method="MOM", ...)
{
    x <- as.vector(r)
    y <- as.vector(vars)
    est <- NA    
    if (length(x)!=length(y)) {warning("The vectors 'r' and 'vars' have a different size!")}
    if (method == "MOM") { est <- uvmetaMOM(x, y) }

    est$call <- match.call()
    class(est) <- "uvmeta"
    est
}

print.uvmeta <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n\n")
    fe = array(NA,dim=c(2,2))
    fe[1,] = cbind(round(x$fixef$mean,5),round(sqrt(x$fixef$var),5))
    fe[2,] = cbind(round(x$ranef$mean,5),round(sqrt(x$ranef$var),5))
    colnames(fe) = c("Estimate","StdErr")
    rownames(fe) = c("Fixed Effects","Random Effects")
    print(fe,quote=F)
    cat("\nHeterogeneity:")
    cat(paste("\nTau squared: \t\t",round(x$ranef$tauSq,5),sep=""))
    cat(paste("\nCochran's Q statistic: \t",round(x$Q$Q,5)," (p.value: ",round(x$Q$p.value,5),")",sep=""))
    cat(paste("\nH-square index: \t", round(x$H2$H2,5),sep=""))
    cat(paste("\nI-square index: \t", round(x$I2$I2*100,3)," %\n",sep=""))
}

predict.uvmeta <- function(object, level = 0.95, ...)
{
  alpha = (1-level)/2

  #The correct number of degrees of freedom for this t distribution is complex, and we use a value of k–2 largely for pragmatic reasons. (Riley 2011)
  df = 2 
  
  pred.mean  <- object$ranef$mean
  pred.lower <- object$ranef$mean + qt(alpha,df=(object$numstudies-df))*sqrt(object$ranef$tauSq+object$ranef$var)
  pred.upper <- object$ranef$mean + qt((1-alpha),df=(object$numstudies-df))*sqrt(object$ranef$tauSq+object$ranef$var)
  predint <- c(pred.mean,pred.lower,pred.upper)
  names(predint) <- c("Estimate", paste((alpha*100),"%"),paste(((1-alpha)*100),"%"))
  predint
}

summary.uvmeta <- function(object, level = 0.95, ...)
{
    alpha = (1-level)/2

    fe.lowerconf =  object$fixef$mean + qnorm(alpha)*sqrt(object$fixef$var)
    fe.upperconf =  object$fixef$mean + qnorm(1-alpha)*sqrt(object$fixef$var)
    re.lowerconf =  object$ranef$mean + qnorm(alpha)*sqrt(object$ranef$var)
    re.upperconf =  object$ranef$mean + qnorm(1-alpha)*sqrt(object$ranef$var)
    lnH.lowerconf = log(sqrt(object$H2$H2)) + qnorm(alpha)*object$H2$se_lnH
    lnH.upperconf = log(sqrt(object$H2$H2)) + qnorm(1-alpha)*object$H2$se_lnH
    H2.lowerconf = (exp(lnH.lowerconf))**2 
    H2.upperconf = (exp(lnH.upperconf))**2
    I2.lowerconf = max(c(0,(H2.lowerconf-1)/H2.lowerconf))
    I2.upperconf = if (H2.upperconf > 0) min(c(1,(H2.upperconf-1)/H2.upperconf)) else 0
    Q.lowerconf = object$Q$Q  + qnorm(alpha)*sqrt(object$Q$var)
    Q.upperconf = object$Q$Q  + qnorm(1-alpha)*sqrt(object$Q$var)
    tauSq.lowerconf = object$ranef$tauSq + qnorm(alpha)*sqrt(object$ranef$varTauSq)
    tauSq.upperconf = object$ranef$tauSq + qnorm(1-alpha)*sqrt(object$ranef$varTauSq)
    lower.conf = c(fe.lowerconf,re.lowerconf,tauSq.lowerconf,Q.lowerconf,H2.lowerconf,I2.lowerconf)
    upper.conf = c(fe.upperconf,re.upperconf,tauSq.upperconf,Q.upperconf,H2.upperconf,I2.upperconf)
    
    TAB = cbind(c(object$fixef$mean,object$ranef$mean,object$ranef$tauSq,object$Q$Q,object$H2$H2,object$I2$I2),lower.conf=lower.conf,upper.conf=upper.conf)
    rownames(TAB) = c("mu (fixed)","mu (random)","Tau squared","Cochran Q","H-square index","I-square index")
    colnames(TAB) = c("Estimate",paste((alpha*100),"%"),paste(((1-alpha)*100),"%"))
    res = list(call=object$call,estimates=TAB)
    class(res) = "summary.macc"
    res
}



