##Set working directory
setwd("C:/Users/Ellen/Desktop/aquatic_data/git/BOSCH/")

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
intra$HW <- as.numeric(intra$Length.of.1st.Antennae)

#log densities
intra$dens <- as.numeric(intra$density_per_m2)
intra$Ldens <- log10(intra$dens+1)

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

###### trying first with simple linear models #### we can discuss what models would be best 
####also still need to get more environmental data, e.g. at least temperature

options(na.action = "na.omit")

####ED
ed_HW_sub<-ed[complete.cases(ed[, "HW"]),]
nrow(ed_HW_sub)
head(ed_HW_sub)
ed_hw <- lm(ed_HW_sub$HW ~ ed_HW_sub$syr + ed_HW_sub$sDOY + ed_HW_sub$Ldens + ed_HW_sub$SiteShort)
ed_hw <- lmer(ed_HW_sub$HW ~ ed_HW_sub$syr + ed_HW_sub$sDOY + ed_HW_sub$Ldens +  (1 | ed_HW_sub$SiteShort)) ##error :/
ed_hw <- gls(HW ~ syr + sDOY + Ldens + SiteShort,na.action=na.omit, data=ed_HW_sub) # did not include a yr autocorr because year is repeated- need to decide how to deal with repeated measures
summary(ed_hw)
plot(ed_HW_sub$HW ~ ed_HW_sub$yr)
abline(lm(ed_HW_sub$HW ~ ed_HW_sub$yr))
abline(lm(ed_HW_sub$HW[ed_HW_sub$Season=="early"] ~ ed_HW_sub$yr[ed_HW_sub$Season=="early"]),lty=2)
abline(lm(ed_HW_sub$HW[ed_HW_sub$Season=="late"] ~ ed_HW_sub$yr[ed_HW_sub$Season=="late"]),lty=2)


ed_BL_sub<-ed[complete.cases(ed[ , "BL"]),]
ed_bl <- lm(ed_BL_sub$BL ~ ed_BL_sub$syr + ed_BL_sub$sDOY + ed_BL_sub$Ldens + ed_BL_sub$SiteShort)
ed_bl <- lmer(ed_BL_sub$BL ~ ed_BL_sub$syr + ed_BL_sub$sDOY + ed_BL_sub$Ldens +(1 | ed_BL_sub$SiteShort))
summary(ed_bl)
plot(ed_BL_sub$BL ~ ed_BL_sub$yr)
abline(lm(ed_BL_sub$BL ~ ed_BL_sub$yr))
abline(lm(ed_BL_sub$BL[ed_BL_sub$Season=="early"] ~ ed_BL_sub$yr[ed_BL_sub$Season=="early"]),lty=2)
abline(lm(ed_BL_sub$BL[ed_BL_sub$Season=="late"] ~ ed_BL_sub$yr[ed_BL_sub$Season=="late"]),lty=2)

####PO
po_HW_sub<-po[complete.cases(po[, "HW"]),]
nrow(po_HW_sub)
po_hw <- lm(po_HW_sub$HW ~ po_HW_sub$syr + po_HW_sub$sDOY + po_HW_sub$Ldens + po_HW_sub$SiteShort)
summary(po_hw)
plot(po_HW_sub$HW ~ po_HW_sub$yr)
abline(lm(po_HW_sub$HW ~ po_HW_sub$yr))

po_BL_sub<-po[complete.cases(po[ , "BL"]),]
po_bl <- lm(po_BL_sub$BL ~ po_BL_sub$syr + po_BL_sub$sDOY + po_BL_sub$Ldens + po_BL_sub$SiteShort)
summary(po_bl)
plot(po_BL_sub$BL ~ po_BL_sub$yr)
abline(lm(po_BL_sub$BL ~ po_BL_sub$yr))

####hs
hs_HW_sub<-hs[complete.cases(hs[, "HW"]),]
nrow(hs_HW_sub)
head(hs_HW_sub)
hs_hw <- lm(hs_HW_sub$HW ~ hs_HW_sub$syr + hs_HW_sub$sDOY + hs_HW_sub$Ldens + hs_HW_sub$SiteShort)
summary(hs_hw)
plot(hs_HW_sub$HW ~ hs_HW_sub$yr)
abline(lm(hs_HW_sub$HW[hs_HW_sub$Season=="early"] ~ hs_HW_sub$yr[hs_HW_sub$Season=="early"]))
abline(lm(hs_HW_sub$HW[hs_HW_sub$Season=="late"] ~ hs_HW_sub$yr[hs_HW_sub$Season=="late"]))

hs_BL_sub<-hs[complete.cases(hs[ , "BL"]),]
hs_bl <- lm(hs_BL_sub$BL ~ hs_BL_sub$syr + hs_BL_sub$sDOY + hs_BL_sub$Ldens + hs_BL_sub$SiteShort)
summary(hs_bl)
plot(hs_BL_sub$BL ~ hs_BL_sub$yr)
abline(lm(hs_BL_sub$BL[hs_BL_sub$Season=="early"] ~ hs_BL_sub$yr[hs_BL_sub$Season=="early"]))
abline(lm(hs_BL_sub$BL[hs_BL_sub$Season=="late"] ~ hs_BL_sub$yr[hs_BL_sub$Season=="late"]))

####aa
aa_HW_sub<-aa[complete.cases(aa[, "HW"]),]
nrow(aa_HW_sub)
head(aa_HW_sub)
aa_hw <- lm(aa_HW_sub$HW ~ aa_HW_sub$syr + aa_HW_sub$sDOY + aa_HW_sub$Ldens + aa_HW_sub$SiteShort + aa_HW_sub$Adult.)
summary(aa_hw)
plot(aa_HW_sub$HW ~ aa_HW_sub$yr)
abline(lm(aa_HW_sub$HW[aa_HW_sub$Season=="early"] ~ aa_HW_sub$yr[aa_HW_sub$Season=="early"]))
abline(lm(aa_HW_sub$HW[aa_HW_sub$Season=="late"] ~ aa_HW_sub$yr[aa_HW_sub$Season=="late"]))
plot(aa_HW_sub$HW ~ aa_HW_sub$Ldens)

aa_BL_sub<-aa[complete.cases(aa[ , "BL"]),]
aa_bl <- lm(aa_BL_sub$BL ~ aa_BL_sub$syr + aa_BL_sub$sDOY + aa_BL_sub$Ldens + aa_BL_sub$SiteShort + aa_BL_sub$Adult.)
summary(aa_bl)
plot(aa_BL_sub$BL ~ aa_BL_sub$yr)
plot(aa_BL_sub$BL ~ aa_BL_sub$Ldens)
abline(lm(aa_BL_sub$BL[aa_BL_sub$Season=="early"] ~ aa_BL_sub$yr[aa_BL_sub$Season=="early"]))
abline(lm(aa_BL_sub$BL[aa_BL_sub$Season=="late"] ~ aa_BL_sub$yr[aa_BL_sub$Season=="late"]))
plot(aa_BL_sub$BL ~ aa_BL_sub$Ldens)
##########################################
####################################################

####model selection
# Example from Burnham and Anderson (2002), page 100:
#  prevent fitting sub-models to different datasets
options(na.action = "na.fail") ##note that if you run this line, the models above will not run anymore

####ED
dd <- dredge(ed_hw)
subset(dd, delta < 2)
#'Best'model
summary(get.models(dd, 1)[[1]])

dd <- dredge(ed_bl)
subset(dd, delta < 2)
#'Best'model
summary(get.models(dd, 1)[[1]])

####PO
dd <- dredge(po_hw)
subset(dd, delta < 2)
#'Best'model
summary(get.models(dd, 1)[[1]])

dd <- dredge(po_bl)
subset(dd, delta < 2)
#'Best'model
summary(get.models(dd, 1)[[1]])

####HS
dd <- dredge(hs_hw)
subset(dd, delta < 2)
#'Best'model
summary(get.models(dd, 1)[[1]])

dd <- dredge(hs_bl)
subset(dd, delta < 2)
#'Best'model
summary(get.models(dd, 1)[[1]])

####aa
dd <- dredge(aa_hw)
subset(dd, delta < 2)
#'Best'model
summary(get.models(dd, 1)[[1]])

dd <- dredge(aa_bl)
subset(dd, delta < 2)
#'Best'model
summary(get.models(dd, 1)[[1]])

##########################################################################################
#############################################################################################

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
points(x=jitter(ed_BL_sub$yr,3), y=ed_BL_sub$BL, pch=(ed_BL_sub$SeasonCode+20), 
bg=alpha((ed_BL_sub$SiteCode),0.6),col=alpha((ed_BL_sub$SiteCode),0.6),lwd=2,cex=1.4)
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
points(x=jitter(po_BL_sub$yr,3), y=po_BL_sub$BL, pch=(po_BL_sub$SeasonCode+20), 
bg=alpha((po_BL_sub$SiteCode),0.6),col=alpha((po_BL_sub$SiteCode),0.6),lwd=2,cex=1.4)
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

