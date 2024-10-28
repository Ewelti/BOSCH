#
#install.packages("effects")
library(effects)
library(stringr)
library(tidyr)
library(scales)
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")
intra <- read.csv("RawData/IntraSppBS.csv", header=T)

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

########calculate BW estimates for a given Temp (marginal effects) for each spp and site
int$spp_site <- paste(int$SiteShort, int$SppCode)
ests_s <- NULL
for(i in unique(int$spp_site)){
  tryCatch({
  sub <- int[int$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub<-sub[complete.cases(sub[, "BW"]),]
	mod<-lm(BW~ Yryly_Temp + sDOY + Ldens, data=sub)
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
gr <- els[which(els$spp=="GR"), ]
af <- els[which(els$spp=="AF"), ]

##plot
tiff(filename = "plots/Temp_BW_BH_LA_wSilhouette.tiff", width = 8, height = 7, units = 'in', res = 600, compression = 'lzw')
par(mar=c(4,5,4,0.5),mfrow=c(2,2))

###############Ancylus fluviatilis
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(af$l95),max(af$u95)), xlim=c(min(af$temp),max(af$temp)))
title(ylab="Body width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(af$temp[af$site=="Auba"],rev(af$temp[af$site=="Auba"])), c(af$l95[af$site=="Auba"],rev(af$u95[af$site=="Auba"])),col=alpha(1,0.4),border=NA)
points(x=af$temp[af$site=="Auba"], y=af$mean[af$site=="Auba"],type="l",col=alpha(1,0.6),lwd=2)
polygon(c(af$temp[af$site=="Bieb"],rev(af$temp[af$site=="Bieb"])), c(af$l95[af$site=="Bieb"],rev(af$u95[af$site=="Bieb"])),col=alpha(2,0.4),border=NA)
points(x=af$temp[af$site=="Bieb"], y=af$mean[af$site=="Bieb"],type="l",col=alpha(2,0.6),lwd=2)
polygon(c(af$temp[af$site=="O3"],rev(af$temp[af$site=="O3"])), c(af$l95[af$site=="O3"],rev(af$u95[af$site=="O3"])),col=alpha(3,0.4),border=NA)
points(x=af$temp[af$site=="O3"], y=af$mean[af$site=="O3"],type="l",col=alpha(3,0.6),lwd=2)
polygon(c(af$temp[af$site=="W1"],rev(af$temp[af$site=="W1"])), c(af$l95[af$site=="W1"],rev(af$u95[af$site=="W1"])),col=alpha(4,0.4),border=NA)
points(x=af$temp[af$site=="W1"], y=af$mean[af$site=="W1"],type="l",col=alpha(4,0.6),lwd=2)
add_silhouette(
  upload_img = "Silhouettes/ancylus_fluviatilis.svg",
  x = 12,
  y = 6.5,
  width = 0.5,
  height = NULL
)
title(main="a. Ancylus fluviatilis",bty="n",cex.main=1.5)

########calculate BH estimates for a given Temp (marginal effects) for each spp and site
ests_s <- NULL
for(i in unique(int$spp_site)){
  tryCatch({
  sub <- int[int$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub<-sub[complete.cases(sub[, "BH"]),]
	mod<-lm(BH~ Yryly_Temp + sDOY + Ldens, data=sub)
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
af <- els[which(els$spp=="AF"), ]

###############Ancylus fluviatilis
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(af$l95),max(af$u95)), xlim=c(min(af$temp),max(af$temp)))
title(ylab="Body height (mm)", line=2.7,cex.lab=1.5)
title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(af$temp[af$site=="Auba"],rev(af$temp[af$site=="Auba"])), c(af$l95[af$site=="Auba"],rev(af$u95[af$site=="Auba"])),col=alpha(1,0.4),border=NA)
points(x=af$temp[af$site=="Auba"], y=af$mean[af$site=="Auba"],type="l",col=alpha(1,0.6),lwd=2)
polygon(c(af$temp[af$site=="Bieb"],rev(af$temp[af$site=="Bieb"])), c(af$l95[af$site=="Bieb"],rev(af$u95[af$site=="Bieb"])),col=alpha(2,0.4),border=NA)
points(x=af$temp[af$site=="Bieb"], y=af$mean[af$site=="Bieb"],type="l",col=alpha(2,0.6),lwd=2)
polygon(c(af$temp[af$site=="O3"],rev(af$temp[af$site=="O3"])), c(af$l95[af$site=="O3"],rev(af$u95[af$site=="O3"])),col=alpha(3,0.4),border=NA)
points(x=af$temp[af$site=="O3"], y=af$mean[af$site=="O3"],type="l",col=alpha(3,0.6),lwd=2)
polygon(c(af$temp[af$site=="W1"],rev(af$temp[af$site=="W1"])), c(af$l95[af$site=="W1"],rev(af$u95[af$site=="W1"])),col=alpha(4,0.4),border=NA)
points(x=af$temp[af$site=="W1"], y=af$mean[af$site=="W1"],type="l",col=alpha(4,0.6),lwd=2)
add_silhouette(
  upload_img = "Silhouettes/ancylus_fluviatilis.svg",
  x = 12,
  y = 3.5,
  width = 0.5,
  height = NULL
)
title(main="b. Ancylus fluviatilis",bty="n",cex.main=1.5)

########calculate LA estimates for a given Temp (marginal effects) for each spp and site
ests_s <- NULL
for(i in unique(int$spp_site)){
  tryCatch({
  sub <- int[int$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub<-sub[complete.cases(sub[, "LA"]),]
	mod<-lm(LA~ Yryly_Temp + sDOY + Ldens, data=sub)
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
gr <- els[which(els$spp=="GR"), ]

###############Gammarus roeselii
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(gr$l95),max(gr$u95)), xlim=c(min(gr$temp),max(gr$temp)))
title(ylab="Antennae length (mm)", line=2.7,cex.lab=1.5)
title(xlab="Temperature (\u00B0C)", line=2.5,cex.lab=1.5)
box(lwd=3)
polygon(c(gr$temp[gr$site=="Auba"],rev(gr$temp[gr$site=="Auba"])), c(gr$l95[gr$site=="Auba"],rev(gr$u95[gr$site=="Auba"])),col=alpha(1,0.4),border=NA)
points(x=gr$temp[gr$site=="Auba"], y=gr$mean[gr$site=="Auba"],type="l",col=alpha(1,0.6),lwd=2)
polygon(c(gr$temp[gr$site=="Bieb"],rev(gr$temp[gr$site=="Bieb"])), c(gr$l95[gr$site=="Bieb"],rev(gr$u95[gr$site=="Bieb"])),col=alpha(2,0.4),border=NA)
points(x=gr$temp[gr$site=="Bieb"], y=gr$mean[gr$site=="Bieb"],type="l",col=alpha(2,0.6),lwd=2)
polygon(c(gr$temp[gr$site=="O3"],rev(gr$temp[gr$site=="O3"])), c(gr$l95[gr$site=="O3"],rev(gr$u95[gr$site=="O3"])),col=alpha(3,0.4),border=NA)
points(x=gr$temp[gr$site=="O3"], y=gr$mean[gr$site=="O3"],type="l",col=alpha(3,0.6),lwd=2)
polygon(c(gr$temp[gr$site=="W1"],rev(gr$temp[gr$site=="W1"])), c(gr$l95[gr$site=="W1"],rev(gr$u95[gr$site=="W1"])),col=alpha(4,0.4),border=NA)
points(x=gr$temp[gr$site=="W1"], y=gr$mean[gr$site=="W1"],type="l",col=alpha(4,0.6),lwd=2)
add_silhouette(
  upload_img = "Silhouettes/Gammarus_roselii.svg",
  x = 12,
  y = 1.5,
  width = 0.4,
  height = NULL
)
title(main="c. Gammarus roeselii",bty="n",cex.main=1.5)

##
plot(1, 1,type= "n",bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n')
legend("top", legend=c("Aubach","Kinzig O3","Kinzig W1"),col=c(1,3,4),pt.lwd=1,lty=1,lwd=8,cex=2,bty="n")

##
dev.off()
##

##
