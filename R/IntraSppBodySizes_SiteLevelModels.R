##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

# load libraries
library(MuMIn)
library(car)
library(nlme)
library(lme4)
library(scales)

# attach data
intra <- read.csv("RawData/IntraSppBS.csv", header=T)
head(intra)
nrow(intra)
unique(intra$Species)

##scale everything
#function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- plyr::numcolwise(scale)(df)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df, newd)
}
#apply function
intra <- scaleVars(intra)
head(intra)

#fix some problem where these variables are not numeric
intra$BL <- as.numeric(intra$Body_Length) 
intra$HW <- as.numeric(intra$Head_Width)
intra$BW <- as.numeric(intra$Body_Width)
intra$BH <- as.numeric(intra$Height)
intra$LA <- as.numeric(intra$Length.of.1st.Antennae)

#log densities
intra$dens <- as.numeric(intra$density_per_m2)
intra$Ldens <- log10(intra$dens+1)

#get initial summary data
library(dplyr)
# Summarize the data by 'Species'
summary_intra <- intra %>%
  group_by(Species) %>%
  summarise(
    count = n(),                       # Count of rows per species
    mean_body_length = mean(Body_Length, na.rm = TRUE),     # Mean body length per species
    mean_head_width = mean(Head_Width, na.rm = TRUE),) %>%  # Mean head width per species
  arrange(mean_body_length)
detach("package:dplyr", unload = TRUE)
summary_intra

#######################################
####################################
#no sci notation
options(scipen = 999)
#################################################################

options(na.action = "na.omit")
unique(intra$SppCode)

#################################YEAR MODELS######################################
head(intra)
######sites overall
######################BL
intra$spp_site <- paste(intra$SppCode, intra$SiteShort) #concatinate spp and site

##test model- first spp
PO_W1 <- intra[which(intra$spp_site=="PO W1"), ]
head(PO_W1)
trend.i <- summary(gls(BL ~ syr + sDOY,na.action=na.omit, data = PO_W1))$tTable[2, c(1,2,4)]
trend.i

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "BL"]),]
  	trend.i <- summary(gls(BL ~ syr + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/YearModels/BodyLength_SITELEVEL_modeloutputs_YEAR.csv")

######################HW

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "HW"]),]
  	trend.i <- summary(gls(HW ~ syr + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/YearModels/HeadWidth_SITELEVEL_modeloutputs_YEAR.csv")

#######################BW

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "BW"]),]
  	trend.i <- summary(gls(BW ~ syr + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(est_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/YearModels/BodyWidth_SITELEVEL_modeloutputs_YEAR.csv")

######################BH

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "BH"]),]
  	trend.i <- summary(gls(BH ~ syr + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/YearModels/BodyHeight_SITELEVEL_modeloutputs_YEAR.csv")

######################LA

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "LA"]),]
  	trend.i <- summary(gls(LA ~ syr + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/YearModels/AntennaeLength_SITELEVEL_modeloutputs_YEAR.csv")

##########################################################################################
#############################################################################################

#################################Temperature MODELS######################################

######################BL

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "BL"]),]
	sub<-sub[complete.cases(sub[, "sYryly_Temp"]),]
  	trend.i <- summary(gls(BL ~ sYryly_Temp + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/TemperatureModels/BodyLength_SITELEVEL_modeloutputs_Temperature.csv")

######################HW

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "HW"]),]
	sub<-sub[complete.cases(sub[, "sYryly_Temp"]),]
  	trend.i <- summary(gls(HW ~ sYryly_Temp + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/TemperatureModels/HeadWidth_SITELEVEL_modeloutputs_Temperature.csv")

#######################BW

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "BW"]),]
	sub<-sub[complete.cases(sub[, "sYryly_Temp"]),]
  	trend.i <- summary(gls(BW ~ sYryly_Temp + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/TemperatureModels/BodyWidth_SITELEVEL_modeloutputs_Temperature.csv")

######################BH

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "BH"]),]
	sub<-sub[complete.cases(sub[, "sYryly_Temp"]),]
  	trend.i <- summary(gls(BH ~ sYryly_Temp + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/TemperatureModels/BodyHeight_SITELEVEL_modeloutputs_Temperature.csv")

######################LA

ests <- NULL
for(i in unique(intra$spp_site)){
  tryCatch({
  sub <- intra[intra$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "LA"]),]
	sub<-sub[complete.cases(sub[, "sYryly_Temp"]),]
  	trend.i <- summary(gls(LA ~ sYryly_Temp + sDOY,na.action=na.omit, data = sub))$tTable[2, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    ests <- rbind(ests, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)
ests
ests[c('spp','site')] <- str_split_fixed(ests$spp_site, ' ', 2)
ests_s = subset(ests, select = -c(spp_site))
head(ests_s)

write.csv(ests_s,"output_data/IntraSpp_ModelOutputs/SiteLevel/TemperatureModels/AntennaeLength_SITELEVEL_modeloutputs_Temperature.csv")

#############################################################
##################################################################
######################################################################
##############################################################################