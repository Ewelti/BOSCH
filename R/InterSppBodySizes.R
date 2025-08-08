##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

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
options(na.action = "na.omit")

########################year
#####model all sites

cwm_m <- lmer(inter$CWM ~ inter$syear + poly(inter$sDOY,2) + (1|inter$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

#####model Aubach

cwm_m <- lm(Au$CWM ~ Au$syear + poly(Au$sDOY,2))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(Au$CWM ~ Au$year)
abline(lm(Au$CWM ~ Au$year))

plot(Au$Lab ~ Au$year)

#####model Bieber
cwm_m <- lm(Bi$CWM ~ Bi$syear + poly(Bi$sDOY,2))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(Bi$CWM ~ Bi$year)
abline(lm(Bi$CWM ~ Bi$year))

plot(Bi$Lab ~ Bi$year)

#####model KiO3
cwm_m <- lm(KiO3$CWM ~ KiO3$syear + poly(KiO3$sDOY,2))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(KiO3$CWM ~ KiO3$year)
abline(lm(KiO3$CWM ~ KiO3$year))

plot(KiO3$Lab ~ KiO3$year)

#####model KiW1
cwm_m <- lm(KiW1$CWM ~ KiW1$syear + poly(KiW1$sDOY,2))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(KiW1$CWM ~ KiW1$year)
abline(lm(KiW1$CWM ~ KiW1$year))

plot(KiW1$Lab ~ KiW1$year)

######################temperature
cwm_m <- lmer(inter$CWM ~ inter$sYryly_Temp + poly(inter$sDOY,2) + (1|inter$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(inter$CWM ~ inter$year)
abline(lm(inter$CWM ~ inter$year))

plot(sub$CWM ~ sub$Yryly_Temp)
plot(inter$Lab ~ inter$year)

#####model Aubach

cwm_m <- lm(Au$CWM ~ Au$sYryly_Temp + poly(Au$sDOY,2))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(Au$CWM ~ Au$Yryly_Temp)
abline(lm(Au$CWM ~ Au$Yryly_Temp))

plot(Au$Lab ~ Au$Yryly_Temp)

#####model Bieber
cwm_m <- lm(Bi$CWM ~ Bi$sYryly_Temp + poly(Bi$sDOY,2))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(Bi$CWM ~ Bi$Yryly_Temp)
abline(lm(Bi$CWM ~ Bi$Yryly_Temp))

plot(Bi$Lab ~ Bi$Yryly_Temp)

#####model KiO3
cwm_m <- lm(KiO3$CWM ~ KiO3$sYryly_Temp + poly(KiO3$sDOY,2))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(KiO3$CWM ~ KiO3$Yryly_Temp)
abline(lm(KiO3$CWM ~ KiO3$Yryly_Temp))

plot(KiO3$Lab ~ KiO3$Yryly_Temp)

#####model KiW1
cwm_m <- lm(KiW1$CWM ~ KiW1$sYryly_Temp + poly(KiW1$sDOY,2))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(KiW1$CWM ~ KiW1$Yryly_Temp)
abline(lm(KiW1$CWM ~ KiW1$Yryly_Temp))

plot(KiW1$Lab ~ KiW1$Yryly_Temp)


########################################