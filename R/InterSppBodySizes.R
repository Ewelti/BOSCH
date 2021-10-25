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

#####model
cwm_m <- lmer(inter$CWM ~ inter$syear + inter$sDOY + inter$Lab + (1|inter$site))
summary(cwm_m)
coefs <- data.frame(coef(summary(cwm_m)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

plot(inter$CWM ~ inter$year)
abline(lm(inter$CWM ~ inter$year))