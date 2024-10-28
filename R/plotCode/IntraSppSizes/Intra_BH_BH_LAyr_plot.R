##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

# load libraries
library(scales)
library(stringr)

# function for adding silhouettes
source("R/Add_silhouette_function.R")

# attach data
intra <- read.csv("RawData/IntraSppBS.csv", header=T)
head(intra)
nrow(intra)
unique(intra$Species)

#fix some problem where these variables are not numeric
intra$BL <- as.numeric(intra$Body_Length) 
intra$HW <- as.numeric(intra$Head_Width)
intra$BW <- as.numeric(intra$Body_Width)
intra$BH <- as.numeric(intra$Height)
intra$LA <- as.numeric(intra$Length.of.1st.Antennae)

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

#log densities
intra$dens <- as.numeric(intra$density_per_m2)
intra$Ldens <- log10(intra$dens+1)

#calculate BW estimates for each spp, site, and year
intra$spp_site_yr <- paste(intra$SppCode, intra$SiteShort,intra$yr)
ests_s <- NULL
for(i in unique(intra$spp_site_yr)){
  tryCatch({
  sub <- intra[intra$spp_site_yr == i, ]
	sub<-sub[complete.cases(sub[, "BW"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
  ests.i <- coef(summary(lm(BW ~ 1 + sDOY + Ldens, data = sub)))[1,1:2]
  ests.i <- data.frame(spp_site_yr = i, t(ests.i))
  ests_s <- rbind(ests_s, ests.i) ; rm(ests.i, sub)
  }, error=function(e){cat(unique(sub$spp_site_yr),conditionMessage(e), "\n")})
} ; rm(i)
ests_s

ests_s[c('spp','site', 'yr')] <- str_split_fixed(ests_s$spp_site_yr, ' ', 3)
colnames(ests_s)[2] ="BW_est"
colnames(ests_s)[3] ="BW_SE"
ests_s$yr <-as.numeric(ests_s$yr)
ests_s = subset(ests_s, select = -c(spp_site_yr))
head(ests_s)
##

#################################
#calculate BH estimates for each spp, site, and year
intra$spp_site_yr <- paste(intra$SppCode, intra$SiteShort,intra$yr)
ests_bh <- NULL
for(i in unique(intra$spp_site_yr)){
  tryCatch({
  sub <- intra[intra$spp_site_yr == i, ]
	sub<-sub[complete.cases(sub[, "BH"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
  ests.i <- coef(summary(lm(BH ~ 1 + sDOY + Ldens, data = sub)))[1,1:2]
  ests.i <- data.frame(spp_site_yr = i, t(ests.i))
  ests_bh <- rbind(ests_bh, ests.i) ; rm(ests.i, sub)
  }, error=function(e){cat(unique(sub$spp_site_yr),conditionMessage(e), "\n")})
} ; rm(i)
ests_bh

ests_bh[c('spp','site', 'yr')] <- str_split_fixed(ests_bh$spp_site_yr, ' ', 3)
colnames(ests_bh)[2] ="BH_est"
colnames(ests_bh)[3] ="BH_SE"
ests_bh$yr <-as.numeric(ests_bh$yr)
ests_bh = subset(ests_bh, select = -c(spp_site_yr))
head(ests_bh)
##

#################################
#calculate LA estimates for each spp, site, and year
intra$spp_site_yr <- paste(intra$SppCode, intra$SiteShort,intra$yr)
ests_la <- NULL
for(i in unique(intra$spp_site_yr)){
  tryCatch({
  sub <- intra[intra$spp_site_yr == i, ]
	sub<-sub[complete.cases(sub[, "LA"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
  ests.i <- coef(summary(lm(LA ~ 1 + sDOY + Ldens, data = sub)))[1,1:2]
  ests.i <- data.frame(spp_site_yr = i, t(ests.i))
  ests_la <- rbind(ests_la, ests.i) ; rm(ests.i, sub)
  }, error=function(e){cat(unique(sub$spp_site_yr),conditionMessage(e), "\n")})
} ; rm(i)
ests_la

ests_la[c('spp','site', 'yr')] <- str_split_fixed(ests_la$spp_site_yr, ' ', 3)
colnames(ests_la)[2] ="LA_est"
colnames(ests_la)[3] ="LA_SE"
ests_la$yr <-as.numeric(ests_la$yr)
ests_la = subset(ests_la, select = -c(spp_site_yr))
head(ests_la)
##
##

#no sci notation
options(scipen = 999)
options(na.action = "na.omit")

##Other measurements plot##
tiff(filename = "plots/Year_BH_BW_LA_wSilhouette.tiff", width = 8, height = 8, units = 'in', res = 600, compression = 'lzw')
par(mar=c(2.5,5,4,0.5),mfrow=c(2,2))


###############Ancylus fluviatilis
ests_s[is.na(ests_s)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ests_s$BW_est-ests_s$BW_SE-0.3),max(ests_s$BW_est + ests_s$BW_SE+0.3)), xlim=c(2000,2020))
title(ylab="Body width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ests_s$yr[ests_s$site=="Auba"]-0.15, y=ests_s$BW_est[ests_s$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ests_s$yr[ests_s$site=="Bieb"]-0.05, y=ests_s$BW_est[ests_s$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ests_s$yr[ests_s$site=="O3"]+0.05, y=ests_s$BW_est[ests_s$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ests_s$yr[ests_s$site=="W1"]+0.15, y=ests_s$BW_est[ests_s$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=ests_s$yr[ests_s$site=="Auba"]-0.15, y=ests_s$BW_est[ests_s$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ests_s$yr[ests_s$site=="Bieb"]-0.05, y=ests_s$BW_est[ests_s$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ests_s$yr[ests_s$site=="O3"]+0.05, y=ests_s$BW_est[ests_s$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ests_s$yr[ests_s$site=="W1"]+0.15, y=ests_s$BW_est[ests_s$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(ests_s$yr[ests_s$site=="Auba"]-0.15, ests_s$BW_est[ests_s$site=="Auba"]-ests_s$BW_SE[ests_s$site=="Auba"], ests_s$yr[ests_s$site=="Auba"]-0.15, ests_s$BW_est[ests_s$site=="Auba"]+ests_s$BW_SE[ests_s$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_s$yr[ests_s$site=="Bieb"]-0.05, ests_s$BW_est[ests_s$site=="Bieb"]-ests_s$BW_SE[ests_s$site=="Bieb"], ests_s$yr[ests_s$site=="Bieb"]-0.05, ests_s$BW_est[ests_s$site=="Bieb"]+ests_s$BW_SE[ests_s$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_s$yr[ests_s$site=="O3"]+0.05, ests_s$BW_est[ests_s$site=="O3"]-ests_s$BW_SE[ests_s$site=="O3"], ests_s$yr[ests_s$site=="O3"]+0.05, ests_s$BW_est[ests_s$site=="O3"]+ests_s$BW_SE[ests_s$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_s$yr[ests_s$site=="W1"]+0.15, ests_s$BW_est[ests_s$site=="W1"]-ests_s$BW_SE[ests_s$site=="W1"], ests_s$yr[ests_s$site=="W1"]+0.15, ests_s$BW_est[ests_s$site=="W1"]+ests_s$BW_SE[ests_s$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/ancylus_fluviatilis.svg",
  x = 2005,
  y = 1.75,
  width = 2.5,
  height = NULL
)
title(main="a. Ancylus fluviatilis",bty="n",cex.main=2)

###############Ancylus fluviatilis
ests_bh[is.na(ests_bh)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ests_bh$BH_est-ests_bh$BH_SE-0.3),max(ests_bh$BH_est + ests_bh$BH_SE+0.3)), xlim=c(2000,2020))
title(ylab="Body height (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ests_bh$yr[ests_bh$site=="Auba"]-0.15, y=ests_bh$BH_est[ests_bh$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ests_bh$yr[ests_bh$site=="Bieb"]-0.05, y=ests_bh$BH_est[ests_bh$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ests_bh$yr[ests_bh$site=="O3"]+0.05, y=ests_bh$BH_est[ests_bh$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ests_bh$yr[ests_bh$site=="W1"]+0.15, y=ests_bh$BH_est[ests_bh$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=ests_bh$yr[ests_bh$site=="Auba"]-0.15, y=ests_bh$BH_est[ests_bh$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ests_bh$yr[ests_bh$site=="Bieb"]-0.05, y=ests_bh$BH_est[ests_bh$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ests_bh$yr[ests_bh$site=="O3"]+0.05, y=ests_bh$BH_est[ests_bh$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ests_bh$yr[ests_bh$site=="W1"]+0.15, y=ests_bh$BH_est[ests_bh$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(ests_bh$yr[ests_bh$site=="Auba"]-0.15, ests_bh$BH_est[ests_bh$site=="Auba"]-ests_bh$BH_SE[ests_bh$site=="Auba"], ests_bh$yr[ests_bh$site=="Auba"]-0.15, ests_bh$BH_est[ests_bh$site=="Auba"]+ests_bh$BH_SE[ests_bh$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_bh$yr[ests_bh$site=="Bieb"]-0.05, ests_bh$BH_est[ests_bh$site=="Bieb"]-ests_bh$BH_SE[ests_bh$site=="Bieb"], ests_bh$yr[ests_bh$site=="Bieb"]-0.05, ests_bh$BH_est[ests_bh$site=="Bieb"]+ests_bh$BH_SE[ests_bh$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_bh$yr[ests_bh$site=="O3"]+0.05, ests_bh$BH_est[ests_bh$site=="O3"]-ests_bh$BH_SE[ests_bh$site=="O3"], ests_bh$yr[ests_bh$site=="O3"]+0.05, ests_bh$BH_est[ests_bh$site=="O3"]+ests_bh$BH_SE[ests_bh$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_bh$yr[ests_bh$site=="W1"]+0.15, ests_bh$BH_est[ests_bh$site=="W1"]-ests_bh$BH_SE[ests_bh$site=="W1"], ests_bh$yr[ests_bh$site=="W1"]+0.15, ests_bh$BH_est[ests_bh$site=="W1"]+ests_bh$BH_SE[ests_bh$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/ancylus_fluviatilis.svg",
  x = 2005,
  y = 0.625,
  width = 2.5,
  height = NULL
)
title(main="b. Ancylus fluviatilis",bty="n",cex.main=2)

###############Gammarus roeselii
ests_la[is.na(ests_la)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ests_la$LA_est-ests_la$LA_SE-1),max(ests_la$LA_est + ests_la$LA_SE+1)), xlim=c(2000,2020))
title(ylab="Antennae length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ests_la$yr[ests_la$site=="Auba"]-0.15, y=ests_la$LA_est[ests_la$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ests_la$yr[ests_la$site=="Bieb"]-0.05, y=ests_la$LA_est[ests_la$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ests_la$yr[ests_la$site=="O3"]+0.05, y=ests_la$LA_est[ests_la$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ests_la$yr[ests_la$site=="W1"]+0.15, y=ests_la$LA_est[ests_la$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=ests_la$yr[ests_la$site=="Auba"]-0.15, y=ests_la$LA_est[ests_la$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ests_la$yr[ests_la$site=="Bieb"]-0.05, y=ests_la$LA_est[ests_la$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ests_la$yr[ests_la$site=="O3"]+0.05, y=ests_la$LA_est[ests_la$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ests_la$yr[ests_la$site=="W1"]+0.15, y=ests_la$LA_est[ests_la$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
arrows(ests_la$yr[ests_la$site=="Auba"]-0.15, ests_la$LA_est[ests_la$site=="Auba"]-ests_la$LA_SE[ests_la$site=="Auba"], ests_la$yr[ests_la$site=="Auba"]-0.15, ests_la$LA_est[ests_la$site=="Auba"]+ests_la$LA_SE[ests_la$site=="Auba"],col=alpha(1,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_la$yr[ests_la$site=="Bieb"]-0.05, ests_la$LA_est[ests_la$site=="Bieb"]-ests_la$LA_SE[ests_la$site=="Bieb"], ests_la$yr[ests_la$site=="Bieb"]-0.05, ests_la$LA_est[ests_la$site=="Bieb"]+ests_la$LA_SE[ests_la$site=="Bieb"],col=alpha(2,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_la$yr[ests_la$site=="O3"]+0.05, ests_la$LA_est[ests_la$site=="O3"]-ests_la$LA_SE[ests_la$site=="O3"], ests_la$yr[ests_la$site=="O3"]+0.05, ests_la$LA_est[ests_la$site=="O3"]+ests_la$LA_SE[ests_la$site=="O3"],col=alpha(3,0.6),lwd=2,length=0.05, angle=90, code=3)
arrows(ests_la$yr[ests_la$site=="W1"]+0.15, ests_la$LA_est[ests_la$site=="W1"]-ests_la$LA_SE[ests_la$site=="W1"], ests_la$yr[ests_la$site=="W1"]+0.15, ests_la$LA_est[ests_la$site=="W1"]+ests_la$LA_SE[ests_la$site=="W1"],col=alpha(4,0.6),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/Gammarus_roselii.svg",
  x = 2002,
  y = 2.5,
  width = 3,
  height = NULL
)
title(main="c. Gammarus roeselii",bty="n",cex.main=2)

plot(1, 1,type= "n",bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n')
legend("topleft", legend=c("Aubach","Bieber","Kinzig O3","Kinzig W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)


dev.off()
##
##


