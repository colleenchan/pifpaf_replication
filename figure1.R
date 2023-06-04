# Generates figure 1: PIF as a function of truncation bound

library(ReIns) # package for generating trunc Weibull r.v.
library(EnvStats) # for generating trunc lognormal
source("pifpaffunctions.R")

M <- seq(1, 100, 1)
paf <- numeric(length(M))
pif25 <- numeric(length(M))
pif50 <- numeric(length(M))
pif75 <- numeric(length(M))

for (i in 1:length(M)){
  paf[i] <- true.pif(beta = log(1.27), p0 = 0, exp = "lognormal",
                    param1 = 0.05,
                    param2 = 0.98,
                    trunc.val = M[i],
                    b = 0)
  pif25[i] <- true.pif(beta = log(1.27), p0 = 0, exp = "lognormal",
                         param1 = 0.05,
                         param2 = 0.98,
                         trunc.val = M[i],
                         b = 0.25)
  pif50[i] <- true.pif(beta = log(1.27), p0 = 0, exp = "lognormal",
                     param1 = 0.05,
                     param2 = 0.98,
                     trunc.val = M[i],
                     b = 0.5)
  pif75[i] <- true.pif(beta = log(1.27), p0 = 0, exp = "lognormal",
                     param1 = 0.05,
                     param2 = 0.98,
                     trunc.val = M[i],
                     b = 0.75)
}


# color blind color palette
cbb <- c("PAF" = "#000000", "PIF g(x) = 0.25x" = "#E69F00",
            "PIF g(x) = 0.5x" = "#56B4E9", "PIF g(x) = 0.75x" = "#009E73")
library(ggplot2)
ggplot(data.frame(M, paf, pif25, pif50, pif75), aes(x = M)) +
  geom_line(aes(y = paf, color = "PAF"), alpha = 0.7, lwd = 1) +
  geom_line(aes(y = pif25, color = "PIF g(x) = 0.25x"), alpha = 0.7, lwd = 1) +
  geom_line(aes(y = pif50, color = "PIF g(x) = 0.5x"), alpha = 0.7, lwd = 1) +
  geom_line(aes(y = pif75, color = "PIF g(x) = 0.75x"), alpha = 0.7, lwd = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Truncation Bound", y = "PIF", color = "") +
  scale_color_manual(values = cbb) + theme_bw()
ggsave("truncboundpif.pdf", width = 5.5, height = 4) # 4 x 4in.

# for presentation version
ggplot(data.frame(M, paf, pif25, pif50, pif75), aes(x = M)) +
  geom_line(aes(y = paf, color = "PAF"), alpha = 0.7, lwd = 1) +
  geom_line(aes(y = pif25, color = "PIF (cft. 75% decrease)"), alpha = 0.7, lwd = 1) +
  geom_line(aes(y = pif50, color = "PIF (cft. 50% decrease)"), alpha = 0.7, lwd = 1) +
  geom_line(aes(y = pif75, color = "PIF (cft. 25% decrease)"), alpha = 0.7, lwd = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Truncation Bound", y = "PIF", color = "") +
  scale_color_manual(values = cbb) + theme_bw()
ggsave("truncboundpif_IBC.pdf", width = 5.5,   height = 3.5)

