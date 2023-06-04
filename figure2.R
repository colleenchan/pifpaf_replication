# Generates figure 2 and values in table 2. Plots distribution and estimates
# PAF using standard, mixture, approximate, and empirical methods in SSB data

library(ggplot2)
source("pifpaffunctions.R")

df <- read.csv("data/ENSANUT2016SSBV2.csv", as.is = TRUE)
x <- df$serving # dist of SSB

# From Dalia's paper, multivariate adj HR = 1.27, 95% CI = (1.16, 1.38)
beta <- log(1.27)
varbeta <- (mean(c(log(1.38) - log(1.27), log(1.27) - log(1.16))) / 
              qnorm(0.975))^2


# PAF estimates
paf.std(mean(x), var(x), "gamma", beta)
(gam <- paf.kehoe(x, "gamma", beta))
paf.kehoe(x, "gamma", beta, trunc.val = 12)

paf.std(mean(x), var(x), "lognormal", beta)
(lnorm <- paf.kehoe(x, "lognormal", beta))
paf.kehoe(x, "lognormal", beta, trunc.val = 12)

paf.std(mean(x), var(x), "normal", beta)
(norm <- paf.kehoe(x, "normal", beta))
paf.kehoe(x, "normal", beta, trunc.val = 12)

paf.std(mean(x), var(x), "weibull", beta)
(weib <- paf.kehoe(x, "weibull", beta))
paf.kehoe(x, "weibull", beta, trunc.val = 12)

pif.ind(x, beta, varbeta)
pif.app(mean(x), var(x), length(x), beta, varbeta)

# plot SSB distribution
ggplot(data.frame(x), aes(x = x)) +
  geom_histogram(mapping = aes(x = x, y = ..density..), binwidth = 0.5) +
  geom_function(fun = function(x) dweibull(x, shape = weib[2], scale = weib[3]),
                aes(color = "Weibull"), size = 1, alpha = 0.8) +
  geom_function(fun = function(x) dgamma(x, shape = gam[2], scale = gam[3]),
                aes(color = "Gamma"), size = 1, alpha = 0.75) +
  geom_function(fun = function(x) dlnorm(x, lnorm[2], lnorm[3]),
                aes(color = "Lognormal"), size = 1, alpha = 0.75) +
  geom_function(fun = function(x) dnorm(x, norm[2], norm[3]),
                aes(color = "Normal"), size = 1, alpha = 0.7) +
  xlim(-2, 12) +
  labs(x = "SSB consumption (servings/day)", color = "Function") + theme_bw()
