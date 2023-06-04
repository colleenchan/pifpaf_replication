# Generates table 1: relative bias of PIF estimators using standard method

library(ReIns) # package for generating trunc Weibull r.v.
library(EnvStats) # for generating trunc lognormal
source("pifpaffunctions.R")


params <- expand.grid(p = 0,
                      true.exp = list(c("Gamma", 1.15, 1.29),
                                      c("Normal", 1.48, 1.38),
                                      c("Weibull", 1.08, 1.53)
                      ))
pafest <- matrix(NA, nrow = nrow(params), ncol = 4)
relbias <- matrix(NA, nrow = nrow(params), ncol = 4)

true.beta <- log(1.27)
varbeta <- (mean(c(log(1.38) - log(1.27), log(1.27) - log(1.16))) / 
              qnorm(0.975))^2
paf <- numeric(nrow(params))
gam.shape <- numeric(nrow(params))
gam.scale <- numeric(nrow(params))
weib.shape <- numeric(nrow(params))
weib.scale <- numeric(nrow(params))
truepaf <- numeric(nrow(params))

for (i in 1:nrow(params)){
  print(paste0("Param setting ", i))
  true.exp <- unlist(params$true.exp[i])
  true.p <- params$p[i]
  truth <- true.paf(beta = true.beta, p0 = true.p, exp = tolower(true.exp[1]),
                    param1 = as.numeric(true.exp[2]),
                    param2 = as.numeric(true.exp[3]),
                    trunc.val = Inf)
  truepaf[i] <- truth$paf
  meanx <- truth$mean
  varx <- truth$var

  gam <- paf.std(meanx = meanx, varx = varx, exposure = "gamma",
                      beta = true.beta)$paf
  lnorm <- paf.std(meanx = meanx, varx = varx, exposure = "lognormal",
                       beta = true.beta)$paf
  norm <- paf.std(meanx = meanx, varx = varx, exposure = "normal",
                     beta = true.beta)$paf
  weib <- paf.std(meanx = meanx, varx = varx, exposure = "weibull",
                        beta = true.beta)$paf
  pafest[i, ] <- c(gamma.std = gam,
                   lnorm.std = lnorm,
                   norm.std = norm,
                   weibull.std = weib)
  relbias[i, ] <- (pafest[i, ] - truepaf[i]) / truepaf[i]

  params$true.exp[i] <- paste0(true.exp[1], "(", true.exp[2], ", ", true.exp[3], ")")
}

tab.bias <- cbind(params, round(truepaf, 4), round(relbias * 100, 1))
colnames(tab.bias) <- c(names(params), "true.paf", "gam.std", "lnorm.std",
                        "norm.std", "weib.std")

library(stargazer)
tab.bias$beta <- NULL
stargazer(as.matrix(tab.bias), header = FALSE)

