logit <- function(x) { if(is.numeric(x))  log(x/(1-x)) else stop("x is not numeric!") }

inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }

generateMCMCinits <- function(n.chains, model.pars)
{
  inits <- list()
  for (i in 1:n.chains) {
    inits.i <- list()
    for (j in 1:length(model.pars)) {
      parname <- model.pars[[j]]$param
      fprior <- model.pars[[j]]$param.f
      fargs <- model.pars[[j]]$param.args
      inits.i[[parname]] = do.call(fprior, fargs)
    }
    inits[[i]] <- inits.i
  }
  return(inits)
}




