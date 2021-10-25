##Set working directory
setwd("C:/Users/ewelti/Desktop/git/BOSCH/")

# load libraries
library(MuMIn)
library(car)
library(nlme)
library(lme4)
library(scales)

# attach data
inter <- read.csv("output_data/CWM.csv", header=T)
head(inter)

#subset by site
Au <- inter[which(inter$site=="Auba"), ]

##scale everything
#function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- plyr::numcolwise(scale)(df)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df, newd)
}
#apply function
inter <- scaleVars(inter)
head(inter)

inter$Lab <- log10(inter$Abundance_per_m2)

#subset by site
Au <- inter[which(inter$site=="Auba"), ]
Bi <- inter[which(inter$site=="Bieb"), ]
KiO3 <- inter[which(inter$site=="KiO3"), ]
KiW1 <- inter[which(inter$site=="KiW1"), ]

#####model all sites
cwm_m <- lmer(inter$CWM ~ inter$syear + inter$sDOY + inter$Lab + (1|inter$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(inter$CWM ~ inter$year)
abline(lm(inter$CWM ~ inter$year))

#####model Aubach

cwm_m <- lmer(Au$CWM ~ Au$syear + Au$sDOY + Au$Lab + (1|Au$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(Au$CWM ~ Au$year)
abline(lm(Au$CWM ~ Au$year))

#####model Bieber
cwm_m <- lmer(Bi$CWM ~ Bi$syear + Bi$sDOY + Bi$Lab + (1|Bi$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(Bi$CWM ~ Bi$year)
abline(lm(Bi$CWM ~ Bi$year))

#####model KiO3
cwm_m <- lmer(KiO3$CWM ~ KiO3$syear + KiO3$sDOY + KiO3$Lab + (1|KiO3$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(KiO3$CWM ~ KiO3$year)
abline(lm(KiO3$CWM ~ KiO3$year))

#####model KiW1
cwm_m <- lmer(KiW1$CWM ~ KiW1$syear + KiW1$sDOY + KiW1$Lab + (1|KiW1$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(KiW1$CWM ~ KiW1$year)
abline(lm(KiW1$CWM ~ KiW1$year))