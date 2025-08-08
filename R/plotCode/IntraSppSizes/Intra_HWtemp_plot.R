#
#install.packages("effects")
library(effects)
library(stringr)
library(tidyr)
library(scales)
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")
intra <- read.csv("RawData/IntraSppBS_updatedR1.csv", header=T)

# function for adding silhouettes
source("R/Add_silhouette_function.R")

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

int <- intra[complete.cases(intra$Yryly_Temp),]
head(int)

########calculate HW estimates for a given Temp (marginal effects) for each spp and site
int$spp_site <- paste(int$SiteShort, int$SppCode)
ests_s <- NULL
for(i in unique(int$spp_site)){
  tryCatch({
  sub <- int[int$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub<-sub[complete.cases(sub[, "HW"]),]
	mod<-lm(HW~ Yryly_Temp + sDOY + Ldens, data=sub)
	model_effects <- effect("Yryly_Temp", mod, xlevels=list(sDOY=mean(sub$sDOY), Ldens=mean(sub$Ldens)))
	suma <- rbind(summary(model_effects)[[3]],summary(model_effects)[[5]],summary(model_effects)[[7]])
	suma <- rbind(colnames(suma),suma)
	rownames(suma) <- c(paste(sub$spp_site[1],"temp"),paste(sub$spp_site[1],"mean"),paste(sub$spp_site[1],"lower95"),paste(sub$spp_site[1],"upper95"))
  	ests_s <- rbind(ests_s, suma) ; rm(suma, sub)
  }, error=function(e){cat(unique(sub$spp_site),conditionMessage(e), "\n")})
} ; rm(i)
ests_s

#############reformating ests
ests_s <- cbind(p=rownames(ests_s), ests_s)
ests_s <-as.data.frame(ests_s)
ests_s[c('site','spp', 'parameter')] <- str_split_fixed(ests_s$p, ' ', 3)
rownames(ests_s) <- NULL
colnames(ests_s) <- NULL
est_s <- cbind(ests_s[,7:9],ests_s[,2:6])
colnames(est_s) <- c("site", "spp", "parameter", "v1", "v2", "v3", "v4", "v5")
est_long <- pivot_longer(est_s, cols = "v1":"v5", names_to = "count",
                                values_to = "values", values_drop_na = TRUE)
el <- as.data.frame(est_long)
ff <- reshape(el, idvar = c("site", "spp", "count"), timevar = "parameter", dir = "wide")
els = subset(ff, select = -c(count) )
head(els)
colnames(els)[3] ="temp"
colnames(els)[4] ="mean"
colnames(els)[5] ="l95"
colnames(els)[6] ="u95"
els$temp <- as.numeric(els$temp)
els$mean <- as.numeric(els$mean)
els$l95 <- as.numeric(els$l95)
els$u95 <- as.numeric(els$u95)
els <- els[complete.cases(els[ , c('l95')]), ]
dim(els)
##############
#subset by spp
ed <- els[which(els$spp=="ED"), ]
po <- els[which(els$spp=="PO"), ]
hs <- els[which(els$spp=="HS"), ]
ov <- els[which(els$spp=="OV"), ]
gr <- els[which(els$spp=="GR"), ]
et <- els[which(els$spp=="ET"), ]
af <- els[which(els$spp=="AF"), ]
br <- els[which(els$spp=="BR"), ]
aa <- els[which(els$spp=="AA"), ]


##################################################################
##############################################################
############################################
#############################

##plot
tiff(filename = "plots/Temp_HW_wSilhouettes.tiff", width = 8, height = 10, units = 'in', res = 600, compression = 'lzw')
par(mar=c(4,5,4,0.5),mfrow=c(4,2))

#####Aphelocheirus aestivalis
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(aa$l95),max(aa$u95)), xlim=c(min(aa$temp),max(aa$temp)))
title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(aa$temp[aa$site=="Auba"],rev(aa$temp[aa$site=="Auba"])), c(aa$l95[aa$site=="Auba"],rev(aa$u95[aa$site=="Auba"])),col=alpha(1,0.2), border=NA)
points(x=aa$temp[aa$site=="Auba"], y=aa$mean[aa$site=="Auba"],type="l",col=alpha(1,0.2), lwd=2)
polygon(c(aa$temp[aa$site=="Bieb"],rev(aa$temp[aa$site=="Bieb"])), c(aa$l95[aa$site=="Bieb"],rev(aa$u95[aa$site=="Bieb"])),col=alpha(2,0.2), border=NA)
points(x=aa$temp[aa$site=="Bieb"], y=aa$mean[aa$site=="Bieb"],type="l",col=alpha(2,0.2), lwd=2)
polygon(c(aa$temp[aa$site=="O3"],rev(aa$temp[aa$site=="O3"])), c(aa$l95[aa$site=="O3"],rev(aa$u95[aa$site=="O3"])),col=alpha(3,0.2), border=NA)
points(x=aa$temp[aa$site=="O3"], y=aa$mean[aa$site=="O3"],type="l",col=alpha(3,0.2), lwd=2)
polygon(c(aa$temp[aa$site=="W1"],rev(aa$temp[aa$site=="W1"])), c(aa$l95[aa$site=="W1"],rev(aa$u95[aa$site=="W1"])),col=alpha(4,0.2), border=NA)
points(x=aa$temp[aa$site=="W1"], y=aa$mean[aa$site=="W1"],type="l",col=alpha(4,0.2), lwd=2)
add_silhouette(
  upload_img = "Silhouettes/aphelocheirus_aestivalis.svg",
  x = 10.5,
  y = 1.37,
  width = 0.15,
  height = NULL
)
title("a. Aphelocheirus aestivalis",bty="n",cex.main=1.5)

###############Baetis rhodani
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(br$l95),max(br$u95)), xlim=c(min(br$temp),max(br$temp)))
#title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(br$temp[br$site=="Auba"],rev(br$temp[br$site=="Auba"])), c(br$l95[br$site=="Auba"],rev(br$u95[br$site=="Auba"])),col=alpha(1,0.8), border=NA)
points(x=br$temp[br$site=="Auba"], y=br$mean[br$site=="Auba"],type="l",col=alpha(1,0.8), lwd=2)
polygon(c(br$temp[br$site=="Bieb"],rev(br$temp[br$site=="Bieb"])), c(br$l95[br$site=="Bieb"],rev(br$u95[br$site=="Bieb"])),col=alpha(2,0.2), border=NA)
points(x=br$temp[br$site=="Bieb"], y=br$mean[br$site=="Bieb"],type="l",col=alpha(2,0.2), lwd=2)
polygon(c(br$temp[br$site=="O3"],rev(br$temp[br$site=="O3"])), c(br$l95[br$site=="O3"],rev(br$u95[br$site=="O3"])),col=alpha(3,0.8), border=NA)
points(x=br$temp[br$site=="O3"], y=br$mean[br$site=="O3"],type="l",col=alpha(3,0.8), lwd=2)
polygon(c(br$temp[br$site=="W1"],rev(br$temp[br$site=="W1"])), c(br$l95[br$site=="W1"],rev(br$u95[br$site=="W1"])),col=alpha(4,0.2), border=NA)
points(x=br$temp[br$site=="W1"], y=br$mean[br$site=="W1"],type="l",col=alpha(4,0.2), lwd=2)
add_silhouette(
  upload_img = "Silhouettes/baetis_rhodani.svg",
  x = 9.5,
  y = 1.25,
  width = 0.5,
  height = NULL
)
title(main="b. Baetis rhodani",bty="n",cex.main=1.5)
         
###############Eiseniella tetraedra
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0.,4), xlim=c(min(et$temp),max(et$temp)))
title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(et$temp[et$site=="Auba"],rev(et$temp[et$site=="Auba"])), c(et$l95[et$site=="Auba"],rev(et$u95[et$site=="Auba"])),col=alpha(1,0.2), border=NA)
points(x=et$temp[et$site=="Auba"], y=et$mean[et$site=="Auba"],type="l",col=alpha(1,0.2), lwd=2)
polygon(c(et$temp[et$site=="Bieb"],rev(et$temp[et$site=="Bieb"])), c(et$l95[et$site=="Bieb"],rev(et$u95[et$site=="Bieb"])),col=alpha(2,0.2), border=NA)
points(x=et$temp[et$site=="Bieb"], y=et$mean[et$site=="Bieb"],type="l",col=alpha(2,0.2), lwd=2)
polygon(c(et$temp[et$site=="O3"],rev(et$temp[et$site=="O3"])), c(et$l95[et$site=="O3"],rev(et$u95[et$site=="O3"])),col=alpha(3,0.8), border=NA)
points(x=et$temp[et$site=="O3"], y=et$mean[et$site=="O3"],type="l",col=alpha(3,0.8), lwd=2)
#polygon(c(et$temp[et$site=="W1"],rev(et$temp[et$site=="W1"])), c(et$l95[et$site=="W1"],rev(et$u95[et$site=="W1"])),col=alpha(4,0.2), border=NA)
#points(x=et$temp[et$site=="W1"], y=et$mean[et$site=="W1"],type="l",col=alpha(4,0.2), lwd=2)
add_silhouette(
  upload_img = "Silhouettes/eiseniella_tetraedra.svg",
  x = 9.6,
  y = 0.5,
  width = 0.3,
  height = NULL
)
title(main="c. Eiseniella tetraedra",bty="n",cex.main=1.5)

###############Ephemera danica
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ed$l95),max(ed$u95)), xlim=c(min(ed$temp),max(ed$temp)))
#title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(ed$temp[ed$site=="Auba"],rev(ed$temp[ed$site=="Auba"])), c(ed$l95[ed$site=="Auba"],rev(ed$u95[ed$site=="Auba"])),col=alpha(1,0.8), border=NA)
points(x=ed$temp[ed$site=="Auba"], y=ed$mean[ed$site=="Auba"],type="l",col=alpha(1,0.8), lwd=2)
polygon(c(ed$temp[ed$site=="Bieb"],rev(ed$temp[ed$site=="Bieb"])), c(ed$l95[ed$site=="Bieb"],rev(ed$u95[ed$site=="Bieb"])),col=alpha(2,0.8), border=NA)
points(x=ed$temp[ed$site=="Bieb"], y=ed$mean[ed$site=="Bieb"],type="l",col=alpha(2,0.8), lwd=2)
polygon(c(ed$temp[ed$site=="O3"],rev(ed$temp[ed$site=="O3"])), c(ed$l95[ed$site=="O3"],rev(ed$u95[ed$site=="O3"])),col=alpha(3,0.2), border=NA)
points(x=ed$temp[ed$site=="O3"], y=ed$mean[ed$site=="O3"],type="l",col=alpha(3,0.2), lwd=2)
polygon(c(ed$temp[ed$site=="W1"],rev(ed$temp[ed$site=="W1"])), c(ed$l95[ed$site=="W1"],rev(ed$u95[ed$site=="W1"])),col=alpha(4,0.8), border=NA)
points(x=ed$temp[ed$site=="W1"], y=ed$mean[ed$site=="W1"],type="l",col=alpha(4,0.8), lwd=2)
add_silhouette(
  upload_img = "Silhouettes/Ephemera_danica.svg",
  x = 9.5,
  y = 2.2,
  width = 0.5,
  height = NULL
)
title(main="d. Ephemera danica ",bty="n",cex.main=1.5)

###############Hydropsyche siltalai
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(hs$l95),max(hs$u95)), xlim=c(min(hs$temp),max(hs$temp)))
title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(hs$temp[hs$site=="Auba"],rev(hs$temp[hs$site=="Auba"])), c(hs$l95[hs$site=="Auba"],rev(hs$u95[hs$site=="Auba"])),col=alpha(1,0.2), border=NA)
points(x=hs$temp[hs$site=="Auba"], y=hs$mean[hs$site=="Auba"],type="l",col=alpha(1,0.2), lwd=2)
polygon(c(hs$temp[hs$site=="Bieb"],rev(hs$temp[hs$site=="Bieb"])), c(hs$l95[hs$site=="Bieb"],rev(hs$u95[hs$site=="Bieb"])),col=alpha(2,0.2), border=NA)
points(x=hs$temp[hs$site=="Bieb"], y=hs$mean[hs$site=="Bieb"],type="l",col=alpha(2,0.2), lwd=2)
polygon(c(hs$temp[hs$site=="O3"],rev(hs$temp[hs$site=="O3"])), c(hs$l95[hs$site=="O3"],rev(hs$u95[hs$site=="O3"])),col=alpha(3,0.2), border=NA)
points(x=hs$temp[hs$site=="O3"], y=hs$mean[hs$site=="O3"],type="l",col=alpha(3,0.2), lwd=2)
polygon(c(hs$temp[hs$site=="W1"],rev(hs$temp[hs$site=="W1"])), c(hs$l95[hs$site=="W1"],rev(hs$u95[hs$site=="W1"])),col=alpha(4,0.2), border=NA)
points(x=hs$temp[hs$site=="W1"], y=hs$mean[hs$site=="W1"],type="l",col=alpha(4,0.2), lwd=2)
add_silhouette(
  upload_img = "Silhouettes/hydropsyche_siltalai.svg",
  x = 9.5,
  y = 1.5,
  width = 0.5,
  height = NULL
)
title(main="e. Hydropsyche siltalai",bty="n",cex.main=1.5)
 
###############Orectochilus villosus
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ov$l95),max(ov$u95)), xlim=c(min(ov$temp),max(ov$temp)))
#title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(ov$temp[ov$site=="Auba"],rev(ov$temp[ov$site=="Auba"])), c(ov$l95[ov$site=="Auba"],rev(ov$u95[ov$site=="Auba"])),col=alpha(1,0.2), border=NA)
points(x=ov$temp[ov$site=="Auba"], y=ov$mean[ov$site=="Auba"],type="l",col=alpha(1,0.2), lwd=2)
polygon(c(ov$temp[ov$site=="Bieb"],rev(ov$temp[ov$site=="Bieb"])), c(ov$l95[ov$site=="Bieb"],rev(ov$u95[ov$site=="Bieb"])),col=alpha(2,0.2), border=NA)
points(x=ov$temp[ov$site=="Bieb"], y=ov$mean[ov$site=="Bieb"],type="l",col=alpha(2,0.2), lwd=2)
polygon(c(ov$temp[ov$site=="O3"],rev(ov$temp[ov$site=="O3"])), c(ov$l95[ov$site=="O3"],rev(ov$u95[ov$site=="O3"])),col=alpha(3,0.2), border=NA)
points(x=ov$temp[ov$site=="O3"], y=ov$mean[ov$site=="O3"],type="l",col=alpha(3,0.2), lwd=2)
polygon(c(ov$temp[ov$site=="W1"],rev(ov$temp[ov$site=="W1"])), c(ov$l95[ov$site=="W1"],rev(ov$u95[ov$site=="W1"])),col=alpha(4,0.2), border=NA)
points(x=ov$temp[ov$site=="W1"], y=ov$mean[ov$site=="W1"],type="l",col=alpha(4,0.2), lwd=2)
add_silhouette(
  upload_img = "Silhouettes/Orectochilus_villosus.svg",
  x = 9.25,
  y = 0.35,
  width = 0.2,
  height = NULL
)
title(main="f. Orectochilus villosus",bty="n",cex.main=1.5)
 
###############Prodiamesa olivacea
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(po$l95),max(po$u95)), xlim=c(min(po$temp),max(po$temp)))
title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(po$temp[po$site=="Auba"],rev(po$temp[po$site=="Auba"])), c(po$l95[po$site=="Auba"],rev(po$u95[po$site=="Auba"])),col=alpha(1,0.2), border=NA)
points(x=po$temp[po$site=="Auba"], y=po$mean[po$site=="Auba"],type="l",col=alpha(1,0.2), lwd=2)
polygon(c(po$temp[po$site=="Bieb"],rev(po$temp[po$site=="Bieb"])), c(po$l95[po$site=="Bieb"],rev(po$u95[po$site=="Bieb"])),col=alpha(2,0.2), border=NA)
points(x=po$temp[po$site=="Bieb"], y=po$mean[po$site=="Bieb"],type="l",col=alpha(2,0.2), lwd=2)
polygon(c(po$temp[po$site=="O3"],rev(po$temp[po$site=="O3"])), c(po$l95[po$site=="O3"],rev(po$u95[po$site=="O3"])),col=alpha(3,0.2), border=NA)
points(x=po$temp[po$site=="O3"], y=po$mean[po$site=="O3"],type="l",col=alpha(3,0.2), lwd=2)
polygon(c(po$temp[po$site=="W1"],rev(po$temp[po$site=="W1"])), c(po$l95[po$site=="W1"],rev(po$u95[po$site=="W1"])),col=alpha(4,0.2), border=NA)
points(x=po$temp[po$site=="W1"], y=po$mean[po$site=="W1"],type="l",col=alpha(4,0.2), lwd=2)
add_silhouette(
  upload_img = "Silhouettes/Prodiamesa_olivacea.svg",
  x = 9,
  y = 0.65,
  width = 0.45,
  height = NULL
)
title(main="g. Prodiamesa olivacea ",bty="n",cex.main=1.5)

##
plot(1, 1,type= "n",bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n')
legend("top", legend=c("Aubach","Bieber","Kinzig O3","Kinzig W1"),col=c(1,2,3,4),pt.lwd=1,lty=1,lwd=8,cex=2,bty="n")
  
##
dev.off()
##

##
