# Generates table 3: relative bias and coverage rate of approximate and 
# empirical estimators

library(ReIns) # package for generating trunc Weibull r.v.
library(EnvStats) # for generating trunc lognormal
source("pifpaffunctions.R")


true.beta <- log(1.27)
true.var.beta <- (mean(c(log(1.38) - log(1.27), log(1.27) - log(1.16))) / 
                    qnorm(0.975))^2
trunc.val <- 12
B <- 10000


start <- Sys.time()
params <- expand.grid(N = c(100, 1000, 10000),
                      p = c(0, 0.05, 0.25, 0.5, 0.75))
                     
paf <- numeric(nrow(params))
emp <-  matrix(NA, nrow = nrow(params), ncol = 5)
app <- matrix(NA, nrow = nrow(params), ncol = 5)
true.exp <- c("weibull", 1.2, 1.66)

# weibull
for (i in 1:nrow(params)){
  print(paste("Param setting", i))
  
  # true params and PAF
  n <- params$N[i]
  true.p <- params$p[i]
  var.beta <- true.var.beta * 70000 / (7*n)
  paf[i] <- unlist(true.paf(beta = true.beta, p0 = true.p, true.exp[1], 
                     as.numeric(true.exp[2]), as.numeric(true.exp[3]), 
                     trunc.val = trunc.val)$paf)
  
  paf.a.ci <- matrix(NA, nrow = B, ncol = 2)
  paf.e.ci <- matrix(NA, nrow = B, ncol = 2)
  paf.a <- numeric(B)
  paf.e <- numeric(B)
  
  for (b in 1:B){
    set.seed(b*i) 
    beta <- rnorm(1, mean = true.beta, sd = sqrt(var.beta))  
    
    set.seed(b*i)
    X <- rbinom(n = n, 1, 1 - true.p) *
      rtweibull(n = n, shape = as.numeric(true.exp[2]),
                scale = as.numeric(true.exp[3]),
                endpoint = trunc.val)
    
    tmppaf.a <- pif.app(meanx = mean(X), varx = var(X), n = n,
                        beta = beta, varbeta = var.beta)
    paf.a[b] <- as.numeric(tmppaf.a[1])
    paf.a.ci[b, ] <- as.numeric(unlist(tmppaf.a[2]))
    
    tmppaf.e <- pif.ind(x = X, beta = beta, varbeta = var.beta)
    paf.e[b] <- as.numeric(tmppaf.e[1])
    paf.e.ci[b, ] <- as.numeric(unlist(tmppaf.e[2]))
  }
  emp[i, ] <- c(avgest = mean(paf.e), 
           avgrb = (mean(paf.e) - paf[i]) / paf[i],
           avgse = mean((paf.e.ci[, 2] - paf.e) / qnorm(0.975)), 
           sd = sd(paf.e),
           cov = mean(paf[i] > paf.e.ci[, 1] & paf[i] < paf.e.ci[, 2]))
  app[i, ] <- c(avgest = mean(paf.a), 
           avgrb = (mean(paf.a) - paf[i]) / paf[i],
           avgse = mean((paf.a.ci[, 2] - paf.a) / qnorm(0.975)), 
           sd = sd(paf.a),
           cov = mean(paf[i] > paf.a.ci[, 1] & paf[i] < paf.a.ci[, 2]))
}
weib <- cbind(params[c(1,3,2)], "true PAF " = round(paf, 4), round(emp, 4), 
              round(app, 4))
names(weib)[5:14] <- c("emp.avgest", "emp.avgrelbias", "emp.se", "emp.sd",
                       "emp.cov", "app.avgest", "app.avgrelbias", "app.se", 
                       "app.sd", "app.cov")




## log normal
true.exp <- c("lognormal", 0.05, 0.98) # log sigma = 0.98
paf <- numeric(nrow(params))
emp <-  matrix(NA, nrow = nrow(params), ncol = 5)
app <- matrix(NA, nrow = nrow(params), ncol = 5)

for (i in 1:nrow(params)){
  print(paste0("Param setting ", i))
  
  # true params and PAF
  n <- params$N[i]
  true.p <- params$p[i]
  var.beta <- true.var.beta * 70000 / (7*n)
  paf[i] <- unlist(true.paf(beta = true.beta, p0 = true.p, true.exp[1], 
                            as.numeric(true.exp[2]), as.numeric(true.exp[3]), 
                            trunc.val = trunc.val)$paf)
  
  paf.a.ci <- matrix(NA, nrow = B, ncol = 2)
  paf.e.ci <- matrix(NA, nrow = B, ncol = 2)
  paf.a <- numeric(B)
  paf.e <- numeric(B)
  
  for (b in 1:B){
    set.seed(b*i) 
    beta <- rnorm(1, mean = true.beta, sd = sqrt(var.beta)) 
    
    set.seed(b*i)
    X <- rbinom(n, 1, 1 - true.p) *
      rlnormTrunc(n, meanlog = as.numeric(true.exp[2]), sdlog = as.numeric(true.exp[3]), 
                  max = trunc.val)
    
    tmppaf.a <- pif.app(meanx = mean(X), varx = var(X), n = n, 
                        beta = beta, varbeta = var.beta)
    paf.a[b] <- as.numeric(tmppaf.a[1])
    paf.a.ci[b, ] <- as.numeric(unlist(tmppaf.a[2]))
    
    tmppaf.e <- pif.ind(x = X, beta = beta, varbeta = var.beta)
    paf.e[b] <- as.numeric(tmppaf.e[1])
    paf.e.ci[b, ] <- as.numeric(unlist(tmppaf.e[2]))
   
  }
  emp[i, ] <- c(avgest = mean(paf.e), 
                avgrb = (mean(paf.e) - paf[i]) / paf[i],
                avgse = mean((paf.e.ci[, 2] - paf.e) / qnorm(0.975)), 
                sd = sd(paf.e),
                cov = mean(paf[i] > paf.e.ci[, 1] & paf[i] < paf.e.ci[, 2]))
  app[i, ] <- c(avgest = mean(paf.a), 
                avgrb = (mean(paf.a) - paf[i]) / paf[i],
                avgse = mean((paf.a.ci[, 2] - paf.a) / qnorm(0.975)), 
                sd = sd(paf.a),
                cov = mean(paf[i] > paf.a.ci[, 1] & paf[i] < paf.a.ci[, 2]))
}
lnorm <- cbind(params[c(1,3,2)], "true PAF " = round(paf, 4), round(emp, 4), 
              round(app, 4))
names(lnorm)[5:14] <- c("emp.avgest", "emp.avgrelbias", "emp.se", "emp.sd",
                        "emp.cov", "app.avgest", "app.avgrelbias", "app.se", 
                        "app.sd", "app.cov")




## normal
library(truncnorm)
true.exp <- c("normal", 1.56, 1.37)
paf <- numeric(nrow(params))
emp <-  matrix(NA, nrow = nrow(params), ncol = 5)
app <- matrix(NA, nrow = nrow(params), ncol = 5)

for (i in 1:nrow(params)){
  print(paste0("Param setting ", i))
  
  # true params and PAF
  n <- params$N[i]
  true.p <- params$p[i]
  var.beta <- true.var.beta * 70000 / (7*n)
  
  paf[i] <- unlist(true.paf(beta = true.beta, p0 = true.p, true.exp[1], 
                            as.numeric(true.exp[2]), sqrt(as.numeric(true.exp[3])), 
                            trunc.val = trunc.val)$paf)
  
  paf.a.ci <- matrix(NA, nrow = B, ncol = 2)
  paf.e.ci <- matrix(NA, nrow = B, ncol = 2)
  paf.a <- numeric(B)
  paf.e <- numeric(B)

  for (b in 1:B){
    set.seed(b*i) 
    beta <- rnorm(1, mean = true.beta, sd = sqrt(var.beta))
    
    set.seed(b*i)
    X <- rbinom(n, 1, 1 - true.p) *
      rtruncnorm(n, mean= as.numeric(true.exp[2]), sd = sqrt(as.numeric(true.exp[3])),
                  a = 0, b = trunc.val)
    
    tmppaf.a <- pif.app(meanx = mean(X), varx = var(X), n = n, 
                        beta = beta, varbeta = var.beta)
    paf.a[b] <- as.numeric(tmppaf.a[1])
    paf.a.ci[b, ] <- as.numeric(unlist(tmppaf.a[2]))
    
    tmppaf.e <- pif.ind(x = X, beta = beta, varbeta = var.beta)
    paf.e[b] <- as.numeric(tmppaf.e[1])
    paf.e.ci[b, ] <- as.numeric(unlist(tmppaf.e[2]))
  }
  emp[i, ] <- c(avgest = mean(paf.e), 
                avgrb = (mean(paf.e) - paf[i]) / paf[i],
                avgse = mean((paf.e.ci[, 2] - paf.e) / qnorm(0.975)), 
                sd = sd(paf.e),
                cov = mean(paf[i] > paf.e.ci[, 1] & paf[i] < paf.e.ci[, 2]))
  app[i, ] <- c(avgest = mean(paf.a), 
                avgrb = (mean(paf.a) - paf[i]) / paf[i],
                avgse = mean((paf.a.ci[, 2] - paf.a) / qnorm(0.975)), 
                sd = sd(paf.a),
                cov = mean(paf[i] > paf.a.ci[, 1] & paf[i] < paf.a.ci[, 2]))
}
norm <- cbind(params[c(1,3,2)], "true PAF " = round(paf, 4), round(emp, 4), 
              round(app, 4))
names(norm)[5:14] <-  c("emp.avgest", "emp.avgrelbias", "emp.se", "emp.sd",
                        "emp.cov", "app.avgest", "app.avgrelbias", "app.se", 
                        "app.sd", "app.cov")
end <- Sys.time()
end - start 

tab <- as.matrix(rbind(lnorm, norm, weib))
library(stargazer)
stargazer(tab, digits = 4)



