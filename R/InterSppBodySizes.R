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

#######################################
####################################
#no sci notation
options(scipen = 999)

#####model all sites
cwm_m <- lmer(inter$CWM ~ inter$syear + inter$sDOY + (1|inter$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(inter$CWM ~ inter$year)
abline(lm(inter$CWM ~ inter$year))

plot(inter$Lab ~ inter$year)
#####model Aubach

cwm_m <- lm(Au$CWM ~ Au$syear + Au$sDOY)
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(Au$CWM ~ Au$year)
abline(lm(Au$CWM ~ Au$year))

plot(Au$Lab ~ Au$year)
#####model Bieber
cwm_m <- lm(Bi$CWM ~ Bi$syear + Bi$sDOY)
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(Bi$CWM ~ Bi$year)
abline(lm(Bi$CWM ~ Bi$year))

plot(Bi$Lab ~ Bi$year)
#####model KiO3
cwm_m <- lm(KiO3$CWM ~ KiO3$syear + KiO3$sDOY)
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(KiO3$CWM ~ KiO3$year)
abline(lm(KiO3$CWM ~ KiO3$year))

plot(KiO3$Lab ~ KiO3$year)
#####model KiW1
cwm_m <- lmer(KiW1$CWM ~ KiW1$syear + KiW1$sDOY)
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(KiW1$CWM ~ KiW1$year)
abline(lm(KiW1$CWM ~ KiW1$year))

plot(KiW1$Lab ~ KiW1$year)

########################################