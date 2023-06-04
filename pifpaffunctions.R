library(MASS)
library(devtools)

devtools::install_github("colleenchan/pifpaf", force = TRUE)
library(pifpaf)


#' Returns closed-form PAF given a parametric distribution
#'
#' @param beta beta coefficient in exponential relative risk
#' @param p0 proportion of zeroes
#' @param exp assumed parametric distribution of exposure
#' @param param1 parameter 1 of parametric distribution
#' @param param2 parameter 2 of parametric distribution
#' @param trunc.val upper truncation bound
#'
#' @return A list of the true PAF value, mean and variance of the exposure
#'
true.paf <- function(beta, p0, exp, param1, param2, trunc.val){
  if (exp == "weibull"){
    int <- tryCatch({
      integrate(function(x) {dweibull(x, param1, param2) * exp(beta*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dweibull(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
    meanx <- param2 * gamma(1 + 1/param1)
    varx <- param2^2 * (gamma(1 + 2/param1)  - (gamma(1 + 1/param1))^2)
  } else if (exp == "lognormal"){
    int <- tryCatch({
      integrate(function(x) {dlnorm(x, param1, param2) * exp(beta*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dlnorm(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
    meanx <- exp(param1 - param2^2)
    varx <- (exp(param2^2) - 1) * exp(2 * param1 + param2^2)
  } else if (exp == "gamma"){
    int <- tryCatch({
      integrate(function(x) {dgamma(x, shape = param1, scale = param2) * exp(beta*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dgamma(x, shape = param1, scale = param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
    meanx <- param1 * param2
    varx <- param1 * param2^2
  } else if (exp == "normal"){
    int <- tryCatch({
      integrate(function(x) {dnorm(x, param1, param2) * exp(beta*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dnorm(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
    meanx <- param1
    varx <- param2^2
  }
  else{
    stop("exp must be either lognormal, normal, gamma, or weibull")
  }

  # E[X] = p(0) + (1-p)E[Xnonzero]
  meanx <- meanx * (1 - p0)

  # Var(X) = E[X^2] - [EX]^2 where
  # E[X^2] = p(0^2) * (1-p)E[Xnonzero^2] where
  # E[Xnonzero^2] = Var(Xnonzero) + E[Xnonzero]^2
  varx <- (varx + meanx^2) * (1-p0) - meanx^2

  return(list(paf = 1 - 1/(p0 + (1-p0)*int), mean = meanx, var = varx))
}



#' Returns closed-form PIF of form g(x) = b * x given a parametric distribution
#'
#' @param beta beta coefficient in exponential relative risk
#' @param p0 proportion of zeroes
#' @param exp parametric distribution 
#' @param param1 parameter 1 of parametric distribution
#' @param param2 parameter 2 of parametric distribution
#' @param trunc.val upper truncation bound
#' @param b slope of linear counterfactual exposure function
#'
#' @return A list of the true PAF value, mean and variance of the exposure
#'
true.pif <- function(beta, 
                     p0, 
                     exp, 
                     param1, 
                     param2, 
                     trunc.val, 
                     b = 0){
  if (exp == "weibull"){
    int <- tryCatch({
      integrate(function(x) {dweibull(x, param1, param2) * exp(beta*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dweibull(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
    intgx <- tryCatch({
      integrate(function(x) {dweibull(x, param1, param2) * exp(beta*b*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dweibull(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
  } else if (exp == "lognormal"){
    int <- tryCatch({
      integrate(function(x) {dlnorm(x, param1, param2) * exp(beta*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dlnorm(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
    intgx <- tryCatch({
      integrate(function(x) {dlnorm(x, param1, param2) * exp(beta*b*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dlnorm(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
  } else if (exp == "gamma"){
    int <- tryCatch({
      integrate(function(x) {dgamma(x, shape = param1, scale = param2) * 
          exp(beta*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dgamma(x, shape = param1, scale = param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
    intgx <- tryCatch({
      integrate(function(x) {dgamma(x, shape = param1, scale = param2) * 
          exp(beta*b*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dgamma(x, shape = param1, scale = param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
  } else if (exp == "normal"){
    int <- tryCatch({
      integrate(function(x) {dnorm(x, param1, param2) * exp(beta*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dnorm(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
    intgx <- tryCatch({
      integrate(function(x) {dnorm(x, param1, param2) * exp(beta*b*x)},
                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dnorm(x, param1, param2)},
                  lower = 0, upper = trunc.val)$value},
      error = function(e){ return(Inf) })
  }
  else{
    stop("exp must be either lognormal, normal, gamma, or weibull")
  }
  return(pif = (int - intgx) / int)
}



#' Standard approach to estimating the PAF. 
#' Matches mean and variance using method of moments
#'
#' @param meanx mean of the exposure data
#' @param varx variance of the exposure data
#' @param exposure parametric distribution 
#' @param beta beta coefficient in exponential relative risk
#'
#' @return A list of the PAF, mean and variance of the exposure
paf.std <- function(meanx, varx, exposure, beta){
  if (exposure == "gamma"){
    par1 <- meanx^2 / varx # shape k
    par2 <- varx / meanx # scale theta
    int <- tryCatch({
      integrate(function(x) {dgamma(x, shape = par1, scale = par2) * exp(beta*x)},
                lower = 0, upper = Inf)$value},
      error = function(e){ return(Inf) })
  }
  if (exposure == "lognormal"){
    par1 <- log(meanx^2 / sqrt(varx + meanx^2)) # log mu
    par2 <- log(sqrt((varx / meanx^2) + 1)) # log sigma
    int <- tryCatch({
      integrate(function(x) {dlnorm(x, par1, par2) * exp(beta*x)},
                lower = 0, upper = Inf, reltol = 1e-20)$value},
      error = function(e){ return(Inf) })
  }
  if (exposure == "normal"){
    par1 <- meanx
    par2 <- sqrt(varx)
    int <- tryCatch({
      integrate(function(x) {dnorm(x, par1, par2) * exp(beta*x)},
                lower = 0, upper = Inf)$value /
        integrate(function(x) {dnorm(x, par1, par2)},
                  lower = 0, upper = Inf)$value},
      error = function(e){ return(Inf) })
    print(int)
  }
  if (exposure == "weibull"){
    mcf <- function(k, x.bar, sd.x) {
      t1 <- sqrt((gamma((k + 2)/k)/(gamma((k + 1)/k)^2)) - 1)
      ((sd.x/x.bar) - t1)^2
    }
    par1 <- nlminb(start = 1, objective = mcf, lower = .Machine$double.eps,
                   x.bar = meanx, sd.x = sqrt(varx))$par # shape k
    par2 <- meanx / gamma((par1 + 1)/par1) # scale lambda
    int <- tryCatch({
      integrate(function(x) {dweibull(x, par1, par2) * exp(beta*x)},
                lower = 0, upper = Inf)$value},
      error = function(e){ return(Inf) })
  }
  paf <- 1 - 1/int
  return(list(paf = paf, param1 = par1, param2 = par2))
}


# Kehoe mixture method
#'
#' @param x vector of exposure values
#' @param exposure assumed parametric distribution of exposure
#' @param beta beta coefficient in exponential relative risk
#' @param trunc.val upper truncation bound
#'
#' @return PAF and mean and variance of exposure
paf.kehoe <- function(x, exposure, beta, trunc.val = Inf){
  p0 <- mean(x == 0)
  x <- x[x != 0]
  mod <- suppressWarnings(as.numeric(fitdistr(x, exposure)$estimate))
  par1 <- mod[1]
  par2 <- mod[2]
  if (exposure == "gamma"){
    int <- tryCatch( 
      {integrate(function(x) {dgamma(x, shape = par1, scale = par2) * exp(beta*x)},
                     lower = 0, upper = trunc.val)$value /
      integrate(function(x) {dgamma(x, par1, par2)},
                lower = 0, upper = trunc.val)$value},
      error = function(e){ Inf })
  }
  if (exposure == "lognormal"){
    int <- tryCatch( 
      {integrate(function(x) {dlnorm(x, par1, par2) * exp(beta*x)},
                     lower = 0, upper = trunc.val)$value /
      integrate(function(x) {dlnorm(x, par1, par2)},
                lower = 0, upper = trunc.val)$value},
      error = function(e){ Inf })
  }
  if (exposure == "normal"){
    int <- tryCatch( 
      {integrate(function(x) {dnorm(x, par1, par2) * exp(beta*x)},
                                lower = 0, upper = trunc.val)$value /
        integrate(function(x) {dnorm(x, par1, par2)},
                  lower = 0, upper = trunc.val)$value},
        error = function(e){ Inf })
  }
  if (exposure == "weibull"){
    int <- tryCatch( 
      {integrate(function(x) {dweibull(x, par1, par2) * exp(beta*x)},
                     lower = 0, upper = trunc.val)$value /
      integrate(function(x) {dweibull(x, par1, par2)},
                lower = 0, upper = trunc.val)$value},
      error = function(e){ Inf })
  }
  paf <- 1 - 1/(p0 + (1 - p0) * int)
  return(c(paf = paf, param1 = par1, param2 = par2))
}
