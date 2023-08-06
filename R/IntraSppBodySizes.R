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

#subset by spp
unique(intra$SppCode)
ed <- intra[which(intra$SppCode=="ED"), ]
nrow(ed)
head(ed)
po <- intra[which(intra$SppCode=="PO"), ]
nrow(po)
hs <- intra[which(intra$SppCode=="HS"), ]
nrow(hs)
ov <- intra[which(intra$SppCode=="OV"), ]
nrow(ov)
gr <- intra[which(intra$SppCode=="GR"), ]
nrow(gr)
head(gr)
et <- intra[which(intra$SppCode=="ET"), ]
nrow(et)
af <- intra[which(intra$SppCode=="AF"), ]
nrow(af)
head(af)
br <- intra[which(intra$SppCode=="BR"), ]
nrow(br)
aa <- intra[which(intra$SppCode=="AA"), ]
nrow(aa)

##check response ditributions
hist(ed$BL)
hist(po$BL)
hist(hs$BL)
hist(ov$BL)
hist(gr$BL)
hist(et$BL)
hist(af$BL)
hist(br$BL)
hist(aa$BL)

hist(ed$HW)
hist(po$HW)
hist(hs$HW)
hist(ov$HW)
hist(gr$Length.of.1st.Antennae) ## different part measured
hist(et$HW)
hist(af$BW) ## different part measured
hist(af$BH) ## different part measured
hist(br$HW)
hist(aa$HW)

hist(ed$Ldens)
hist(po$Ldens)
hist(hs$Ldens)
hist(ov$Ldens)
hist(gr$Ldens)
hist(et$Ldens)
hist(af$Ldens)
hist(br$Ldens)
hist(aa$Ldens)
#######################################
####################################
#no sci notation
options(scipen = 999)
#################################################################

options(na.action = "na.omit")
unique(intra$SppCode)

#################################YEAR MODELS######################################

######################BL

ests <- NULL
for(i in unique(intra$SppCode)){
  sub <- intra[intra$SppCode == i, ]
	sub<-sub[complete.cases(sub[, "BL"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub$SiteShort <- as.factor(sub$SiteShort)
  	coefs <- data.frame(coef(summary(lmer(BL ~ syr + sDOY + Ldens + (1|SiteShort), data = sub))))
	coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
	colnames(coefs)[1] ="Est"
	colnames(coefs)[2] ="SE"
	colnames(coefs)[3] ="t"
	colnames(coefs)[4] ="p"
	ests.i <- coefs[2:4,1:4]
  ests.i <- data.frame(SppCode = i, t(ests.i))
  ests <- rbind(ests, ests.i) ; rm(ests.i, sub)
} ; rm(i)
ests

write.csv(ests,"output_data/IntraSpp_ModelOutputs/YearModels/BodyLength_modeloutputs_YEAR.csv")

######################HW

ests <- NULL
for(i in unique(intra$SppCode)){
  tryCatch({
  sub <- intra[intra$SppCode == i, ]
	sub<-sub[complete.cases(sub[, "HW"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub$SiteShort <- as.factor(sub$SiteShort)
  	coefs <- data.frame(coef(summary(lmer(HW ~ syr + sDOY + Ldens + (1|SiteShort), data = sub))))
	coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
	colnames(coefs)[1] ="Est"
	colnames(coefs)[2] ="SE"
	colnames(coefs)[3] ="t"
	colnames(coefs)[4] ="p"
	ests.i <- coefs[2:4,1:4]
  ests.i <- data.frame(SppCode = i, t(ests.i))
  ests <- rbind(ests, ests.i) ; rm(ests.i, sub)
    }, error=function(e){cat(unique(sub$SppCode),conditionMessage(e), "\n")})
} ; rm(i)
ests

write.csv(ests,"output_data/IntraSpp_ModelOutputs/YearModels/HeadWidth_modeloutputs_YEAR.csv")

#######################BW

BW_sub<-af[complete.cases(af[, "BW"]),]
BW_sub<-BW_sub[complete.cases(BW_sub[, "Ldens"]),]
BW_sub$SiteShort <- as.factor(BW_sub$SiteShort)
af_BW <- lmer(BW_sub$BW ~ BW_sub$syr + BW_sub$sDOY + BW_sub$Ldens + (1|BW_sub$SiteShort))
coefs <- data.frame(coef(summary(af_BW)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
colnames(coefs)[1] ="Est"
colnames(coefs)[2] ="SE"
colnames(coefs)[3] ="t"
colnames(coefs)[4] ="p"
tco <- t(coefs)
tco

write.csv(tco,"output_data/IntraSpp_ModelOutputs/YearModels/BodyWidth_modeloutputs_YEAR.csv")

######################BH

BH_sub<-af[complete.cases(af[, "BH"]),]
BH_sub<-BH_sub[complete.cases(BH_sub[, "Ldens"]),]
BH_sub$SiteShort <- as.factor(BH_sub$SiteShort)
af_BH <- lmer(BH_sub$BH ~ BH_sub$syr + BH_sub$sDOY + BH_sub$Ldens + (1|BH_sub$SiteShort))
coefs <- data.frame(coef(summary(af_BH)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
colnames(coefs)[1] ="Est"
colnames(coefs)[2] ="SE"
colnames(coefs)[3] ="t"
colnames(coefs)[4] ="p"
tco <- t(coefs)
tco

write.csv(tco,"output_data/IntraSpp_ModelOutputs/YearModels/BodyHeight_modeloutputs_YEAR.csv")

######################LA

LA_sub<-gr[complete.cases(gr[, "LA"]),]
LA_sub<-LA_sub[complete.cases(LA_sub[, "Ldens"]),]
LA_sub$SiteShort <- as.factor(LA_sub$SiteShort)
gr_LA <- lmer(LA_sub$LA ~ LA_sub$syr + LA_sub$sDOY + LA_sub$Ldens + (1|LA_sub$SiteShort))
coefs <- data.frame(coef(summary(gr_LA)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
colnames(coefs)[1] ="Est"
colnames(coefs)[2] ="SE"
colnames(coefs)[3] ="t"
colnames(coefs)[4] ="p"
tco <- t(coefs)
tco

write.csv(tco,"output_data/IntraSpp_ModelOutputs/YearModels/AntennaeLength_modeloutputs_YEAR.csv")

##########################################################################################
#############################################################################################

#################################Temperature MODELS######################################

######################BL

ests <- NULL
for(i in unique(intra$SppCode)){
  sub <- intra[intra$SppCode == i, ]
	sub<-sub[complete.cases(sub[, "BL"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub<-sub[complete.cases(sub[, "sYryly_Temp"]),]
	sub$SiteShort <- as.factor(sub$SiteShort)
  	coefs <- data.frame(coef(summary(lmer(BL ~ sYryly_Temp + sDOY + Ldens + (1|SiteShort), data = sub))))
	coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
	colnames(coefs)[1] ="Est"
	colnames(coefs)[2] ="SE"
	colnames(coefs)[3] ="t"
	colnames(coefs)[4] ="p"
	ests.i <- coefs[2:4,1:4]
  ests.i <- data.frame(SppCode = i, t(ests.i))
  ests <- rbind(ests, ests.i) ; rm(ests.i, sub)
} ; rm(i)
ests

write.csv(ests,"output_data/IntraSpp_ModelOutputs/TemperatureModels/BodyLength_modeloutputs_Temperature.csv")

######################HW

ests <- NULL
for(i in unique(intra$SppCode)){
  tryCatch({
  sub <- intra[intra$SppCode == i, ]
	sub<-sub[complete.cases(sub[, "HW"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub<-sub[complete.cases(sub[, "sYryly_Temp"]),]
	sub$SiteShort <- as.factor(sub$SiteShort)
  	coefs <- data.frame(coef(summary(lmer(HW ~ sYryly_Temp + sDOY + Ldens + (1|SiteShort), data = sub))))
	coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
	colnames(coefs)[1] ="Est"
	colnames(coefs)[2] ="SE"
	colnames(coefs)[3] ="t"
	colnames(coefs)[4] ="p"
	ests.i <- coefs[2:4,1:4]
  ests.i <- data.frame(SppCode = i, t(ests.i))
  ests <- rbind(ests, ests.i) ; rm(ests.i, sub)
    }, error=function(e){cat(unique(sub$SppCode),conditionMessage(e), "\n")})
} ; rm(i)
ests

write.csv(ests,"output_data/IntraSpp_ModelOutputs/TemperatureModels/HeadWidth_modeloutputs_Temperature.csv")

#######################BW

BW_sub<-af[complete.cases(af[, "BW"]),]
BW_sub<-BW_sub[complete.cases(BW_sub[, "Ldens"]),]
BW_sub<-BW_sub[complete.cases(BW_sub[, "sYryly_Temp"]),]
BW_sub$SiteShort <- as.factor(BW_sub$SiteShort)
af_BW <- lmer(BW_sub$BW ~ BW_sub$sYryly_Temp + BW_sub$sDOY + BW_sub$Ldens + (1|BW_sub$SiteShort))
coefs <- data.frame(coef(summary(af_BW)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
colnames(coefs)[1] ="Est"
colnames(coefs)[2] ="SE"
colnames(coefs)[3] ="t"
colnames(coefs)[4] ="p"
tco <- t(coefs)

write.csv(tco,"output_data/IntraSpp_ModelOutputs/TemperatureModels/BodyWidth_modeloutputs_Temperature.csv")

######################BH

BH_sub<-af[complete.cases(af[, "BH"]),]
BH_sub<-BH_sub[complete.cases(BH_sub[, "Ldens"]),]
BH_sub<-BH_sub[complete.cases(BH_sub[, "sYryly_Temp"]),]
BH_sub$SiteShort <- as.factor(BH_sub$SiteShort)
af_BH <- lmer(BH_sub$BH ~ BH_sub$sYryly_Temp + BH_sub$sDOY + BH_sub$Ldens + (1|BH_sub$SiteShort))
coefs <- data.frame(coef(summary(af_BH)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
colnames(coefs)[1] ="Est"
colnames(coefs)[2] ="SE"
colnames(coefs)[3] ="t"
colnames(coefs)[4] ="p"
tco <- t(coefs)

write.csv(tco,"output_data/IntraSpp_ModelOutputs/TemperatureModels/BodyHeight_modeloutputs_Temperature.csv")

######################LA

LA_sub<-gr[complete.cases(gr[, "LA"]),]
LA_sub<-LA_sub[complete.cases(LA_sub[, "Ldens"]),]
LA_sub<-LA_sub[complete.cases(LA_sub[, "sYryly_Temp"]),]
LA_sub$SiteShort <- as.factor(LA_sub$SiteShort)
gr_LA <- lmer(LA_sub$LA ~ LA_sub$sYryly_Temp + LA_sub$sDOY + LA_sub$Ldens + (1|LA_sub$SiteShort))
coefs <- data.frame(coef(summary(gr_LA)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
colnames(coefs)[1] ="Est"
colnames(coefs)[2] ="SE"
colnames(coefs)[3] ="t"
colnames(coefs)[4] ="p"
tco <- t(coefs)

write.csv(tco,"output_data/IntraSpp_ModelOutputs/TemperatureModels/AntennaeLength_modeloutputs_Temperature.csv")

#############################################################
##################################################################
######################################################################
##############################################################################

####BL YEAR plots
#######################################
####################################
#no sci notation
options(scipen = 999)
options(na.action = "na.omit")

tiff(filename = "plots/Year_BL.tiff", width = 10, height = 6, units = 'in', res = 600, compression = 'lzw')
par(mar=c(2,4,4,0.4),mfrow=c(3,3))


#####Aphelocheirus aestivalis

head(aa)
ad_aa <-aa[which(aa$Adult.=="yes"),]

plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ad_aa$BL),max(ad_aa$BL)), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
title(main="a. Aphelocheirus aestivalis adults", line=0.5,cex.lab=1.5)
box(lwd=3)

points(x=ad_aa$yr[ad_aa$SiteShort=="Auba"], y=ad_aa$BL[ad_aa$SiteShort=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ad_aa$yr[ad_aa$SiteShort=="Bieb"], y=ad_aa$BL[ad_aa$SiteShort=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ad_aa$yr[ad_aa$SiteShort=="O3"], y=ad_aa$BL[ad_aa$SiteShort=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ad_aa$yr[ad_aa$SiteShort=="W1"], y=ad_aa$BL[ad_aa$SiteShort=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
abline(lm(ad_aa$BL ~ ad_aa$yr), lwd=4, lty=2)
#abline(lm(ad_aa$BL[ad_aa$SiteShort=="Auba"] ~ ad_aa$yr[ad_aa$SiteShort=="Auba"]),lwd=2,col=alpha(1,0.6),lty=1)
#abline(lm(ad_aa$BL[ad_aa$SiteShort=="Bieb"] ~ ad_aa$yr[ad_aa$SiteShort=="Bieb"]),lwd=2,col=alpha(2,0.6),lty=1)
#abline(lm(ad_aa$BL[ad_aa$SiteShort=="O3"] ~ ad_aa$yr[ad_aa$SiteShort=="O3"]),lwd=2,col=alpha(3,0.6),lty=1)
#abline(lm(ad_aa$BL[ad_aa$SiteShort=="W1"] ~ ad_aa$yr[ad_aa$SiteShort=="W1"]),lwd=2,col=alpha(4,0.6),lty=1)
#legend("topleft", legend=c("Aubach","Bieber","KiO3","KiW1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2.5, cex=1.5)

ju_aa <-aa[which(aa$Adult.=="no"),]

plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ju_aa$BL),max(ju_aa$BL)), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
title(main="b. Aphelocheirus aestivalis immatures", line=0.5,cex.lab=1.5)
box(lwd=3)

points(x=ju_aa$yr[ju_aa$SiteShort=="Auba"], y=ju_aa$BL[ju_aa$SiteShort=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ju_aa$yr[ju_aa$SiteShort=="Bieb"], y=ju_aa$BL[ju_aa$SiteShort=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ju_aa$yr[ju_aa$SiteShort=="O3"], y=ju_aa$BL[ju_aa$SiteShort=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ju_aa$yr[ju_aa$SiteShort=="W1"], y=ju_aa$BL[ju_aa$SiteShort=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
abline(lm(ju_aa$BL ~ ju_aa$yr), lwd=4, lty=2)
#abline(lm(ju_aa$BL[ju_aa$SiteShort=="Auba"] ~ ju_aa$yr[ju_aa$SiteShort=="Auba"]),lwd=2,col=alpha(1,0.6),lty=1)
#abline(lm(ju_aa$BL[ju_aa$SiteShort=="Bieb"] ~ ju_aa$yr[ju_aa$SiteShort=="Bieb"]),lwd=2,col=alpha(2,0.6),lty=1)
#abline(lm(ju_aa$BL[ju_aa$SiteShort=="O3"] ~ ju_aa$yr[ju_aa$SiteShort=="O3"]),lwd=2,col=alpha(3,0.6),lty=1)
#abline(lm(ju_aa$BL[ju_aa$SiteShort=="W1"] ~ ju_aa$yr[ju_aa$SiteShort=="W1"]),lwd=2,col=alpha(4,0.6),lty=1)
#legend("topleft", legend=c("Aubach","Bieber","KiO3","KiW1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2.5, cex=1.5)

#####Ancylus fluviatilis

plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(af$BL),max(af$BL)), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
title(main="c. Ancylus fluviatilis", line=0.5,cex.lab=1.5)
box(lwd=3)

points(x=af$yr[af$SiteShort=="O3"], y=af$BL[af$SiteShort=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$SiteShort=="W1"], y=af$BL[af$SiteShort=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$SiteShort=="Auba"], y=af$BL[af$SiteShort=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$SiteShort=="Bieb"], y=af$BL[af$SiteShort=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
abline(lm(af$BL ~ af$yr), lwd=4, lty=2)
#abline(lm(af$BL[af$SiteShort=="O3"] ~ af$yr[af$SiteShort=="O3"]),lwd=2,col=alpha(3,0.6),lty=2)
#abline(lm(af$BL[af$SiteShort=="W1"] ~ af$yr[af$SiteShort=="W1"]),lwd=2,col=alpha(4,0.6),lty=2)
#abline(lm(af$BL[af$SiteShort=="Auba"] ~ af$yr[af$SiteShort=="Auba"]),lwd=2,col=alpha(1,0.6),lty=1)
#abline(lm(af$BL[af$SiteShort=="Bieb"] ~ af$yr[af$SiteShort=="Bieb"]),lwd=2,col=alpha(2,0.6),lty=1)

#####Baetis rhodani

head(br)
br_sub<-br[complete.cases(br[, "BL"]),]
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(br_sub$BL),max(br_sub$BL)), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
title(main="d. Baetis rhodani", line=0.5,cex.lab=1.5)
box(lwd=3)

points(jitter(x=br$yr[br$SiteShort=="O3"],3), y=br$BL[br$SiteShort=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(jitter(x=br$yr[br$SiteShort=="W1"],3), y=br$BL[br$SiteShort=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(jitter(x=br$yr[br$SiteShort=="Auba"],3), y=br$BL[br$SiteShort=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(jitter(x=br$yr[br$SiteShort=="Bieb"],3), y=br$BL[br$SiteShort=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
#abline(lm(br$BL ~ br$yr), lwd=4, lty=2)
abline(lm(br$BL[br$SiteShort=="O3"] ~ br$yr[br$SiteShort=="O3"]),lwd=2,col=alpha(3,0.6),lty=2)
abline(lm(br$BL[br$SiteShort=="W1"] ~ br$yr[br$SiteShort=="W1"]),lwd=2,col=alpha(4,0.6),lty=2)
abline(lm(br$BL[br$SiteShort=="Auba"] ~ br$yr[br$SiteShort=="Auba"]),lwd=2,col=alpha(1,0.6),lty=1)
abline(lm(br$BL[br$SiteShort=="Bieb"] ~ br$yr[br$SiteShort=="Bieb"]),lwd=2,col=alpha(2,0.6),lty=1)

dev.off()

##
##















##first visualizations of intraspecific body sizes over years

####ED
##HW
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,3), xlim=c(2000,2020))
title(ylab="Head width (mm)", line=2.5,cex.lab=1.5) ##title(ylab="", line=2.5,cex.lab=1.5)
title(main="Ephemera danica", line=1)
points(x=jitter(ed_HW_sub$yr,3), y=ed_HW_sub$HW, pch=(ed_HW_sub$SeasonCode+20), 
bg=alpha((ed_HW_sub$SiteCode),0.6),col=alpha((ed_HW_sub$SiteCode),0.6),lwd=2,cex=1.4)
legend("bottomleft", legend=c("Auba","Bieb","O3","W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)
legend("bottom", legend=c("early","late"), pch=c(21,22), bty="n", cex=1.5,pt.lwd=2)

##BL
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,34), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.5,cex.lab=1.5)
title(main="Ephemera danica", line=1)
points(x=jitter(ed$yr,3), y=ed$BL, pch=(ed$SeasonCode+20), 
bg=alpha((ed$SiteCode),0.6),col=alpha((ed$SiteCode),0.6),lwd=2,cex=1.4)
legend("bottomleft", legend=c("Auba","Bieb","O3","W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)
legend("bottom", legend=c("early","late"), pch=c(21,22), bty="n", cex=1.5,pt.lwd=2)

####PO
##HW
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0.1,1.2), xlim=c(2000,2020))
title(ylab="Head width (mm)", line=2.5,cex.lab=1.5) ##title(ylab="", line=2.5,cex.lab=1.5)
title(main="Prodiamesa olivacea", line=1)
points(x=jitter(po_HW_sub$yr,3), y=po_HW_sub$HW, pch=(po_HW_sub$SeasonCode+20), 
bg=alpha((po_HW_sub$SiteCode),0.6),col=alpha((po_HW_sub$SiteCode),0.6),lwd=2,cex=1.4)
legend("topleft", legend=c("Auba","Bieb","O3","W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)
legend("top", legend=c("early","late"), pch=c(21,22), bty="n", cex=1.5,pt.lwd=2)

##BL
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,20), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.5,cex.lab=1.5)
title(main="Prodiamesa olivacea", line=1)
points(x=jitter(po$yr,3), y=po$BL, pch=(po$SeasonCode+20), 
bg=alpha((po$SiteCode),0.6),col=alpha((po$SiteCode),0.6),lwd=2,cex=1.4)
legend("topleft", legend=c("Auba","Bieb","O3","W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)
legend("top", legend=c("early","late"), pch=c(21,22), bty="n", cex=1.5,pt.lwd=2)

####HS
##HW
max(hs_HW_sub$HW)
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0.45,1.9), xlim=c(2000,2020))
title(ylab="Head width (mm)", line=2.5,cex.lab=1.5) ##title(ylab="", line=2.5,cex.lab=1.5)
title(main="Hydropsyche siltalai", line=1)
points(x=jitter(hs_HW_sub$yr,3), y=hs_HW_sub$HW, pch=(hs_HW_sub$SeasonCode+20), 
bg=alpha((hs_HW_sub$SiteCode),0.6),col=alpha((hs_HW_sub$SiteCode),0.6),lwd=2,cex=1.4)
legend("topright", legend=c("Auba","Bieb","O3","W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)
legend("top", legend=c("early","late"), pch=c(21,22), bty="n", cex=1.5,pt.lwd=2)

##BL
min(hs_BL_sub$BL)
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(3,20), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.5,cex.lab=1.5)
title(main="Hydropsyche siltalai", line=1)
points(x=jitter(hs_BL_sub$yr,3), y=hs_BL_sub$BL, pch=(hs_BL_sub$SeasonCode+20), 
bg=alpha((hs_BL_sub$SiteCode),0.6),col=alpha((hs_BL_sub$SiteCode),0.6),lwd=2,cex=1.4)
legend("topright", legend=c("Auba","Bieb","O3","W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)
legend("top", legend=c("early","late"), pch=c(21,22), bty="n", cex=1.5,pt.lwd=2)

####aa
##HW
ad<-aa_HW_sub[aa_HW_sub$Adult.=="yes",]
ju<-aa_HW_sub[aa_HW_sub$Adult.=="no",]
max(aa_HW_sub$HW)
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,2.4), xlim=c(1999,2020))
title(ylab="Head width (mm)", line=2.5,cex.lab=1.5) ##title(ylab="", line=2.5,cex.lab=1.5)
title(main="Aphelocheirus aestivalis", line=1)
points(x=jitter(ad$yr,3), y=ad$HW, pch=(ad$SeasonCode+20), 
       bg=alpha((ad$SiteCode),0.6),col=alpha((ad$SiteCode),0.6),lwd=2,cex=1.4)
points(x=jitter(ju$yr,3), y=ju$HW, pch=(ju$SeasonCode+20), 
       bg=alpha((ju$SiteCode+2),0.6),col=alpha((ju$SiteCode+2),0.6),lwd=2,cex=1.4)
legend("topleft", legend=c("O3 adult","W1 adult","O3 juvenile","W1 juvenile"),col=c(3,4,5,6),pt.bg=c(3,4,5,6),pt.lwd=1, pch=c(24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)
legend("left", legend=c("early","late"), pch=c(21,22), bty="n", cex=1.5,pt.lwd=2)

##BL
ad<-aa_BL_sub[aa_BL_sub$Adult.=="yes",]
ju<-aa_BL_sub[aa_BL_sub$Adult.=="no",]
max(aa_BL_sub$BL)
head(aa_BL_sub)
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(1,11), xlim=c(2000,2020))
title(ylab="Body length (mm)", line=2.5,cex.lab=1.5)
title(main="Aphelocheirus aestivalis", line=1)
points(x=jitter(ad$yr,3), y=ad$BL, pch=(ad$SeasonCode+20), 
       bg=alpha((ad$SiteCode),0.6),col=alpha((ad$SiteCode),0.6),lwd=2,cex=1.4)
points(x=jitter(ju$yr,3), y=ju$BL, pch=(ju$SeasonCode+20), 
       bg=alpha((ju$SiteCode+2),0.6),col=alpha((ju$SiteCode+2),0.6),lwd=2,cex=1.4)
legend("topleft", legend=c("O3 adult","W1 adult","O3 juvenile","W1 juvenile"),col=c(3,4,5,6),pt.bg=c(3,4,5,6),pt.lwd=1, pch=c(24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)
legend("left", legend=c("early","late"), pch=c(21,22), bty="n", cex=1.5,pt.lwd=2)
########################################
###########################################################################################################################3

