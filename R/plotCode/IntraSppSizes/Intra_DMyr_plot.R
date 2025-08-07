##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

# load libraries
library(scales)
library(stringr)
library(vegan)

# function for adding silhouettes
source("R/Add_silhouette_function.R")

# attach data
intra <- read.csv("RawData/IntraSppBS_updatedR1.csv", header=T)
head(intra)
nrow(intra)
unique(intra$Species)

#fix some problem where these variables are not numeric
intra$DM <- as.numeric(intra$dry_mass_mg)

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

#calculate DM estimates for each spp, site, and year
intra$spp_site_yr <- paste(intra$SppCode, intra$SiteShort,intra$yr)
ests_s <- NULL
for(i in unique(intra$spp_site_yr)){
  tryCatch({
  sub <- intra[intra$spp_site_yr == i, ]
	sub<-sub[complete.cases(sub[, "DM"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
  ests.i <- coef(summary(lm(DM ~ 1 + sDOY + Ldens, data = sub)))[1,1:2]
  ests.i <- data.frame(spp_site_yr = i, t(ests.i))
  ests_s <- rbind(ests_s, ests.i) ; rm(ests.i, sub)
  }, error=function(e){cat(unique(sub$spp_site_yr),conditionMessage(e), "\n")})
} ; rm(i)
ests_s

ests_s[c('spp','site', 'yr')] <- str_split_fixed(ests_s$spp_site_yr, ' ', 3)
colnames(ests_s)[2] ="DM_est"
colnames(ests_s)[3] ="DM_SE"
ests_s$yr <-as.numeric(ests_s$yr)
ests_s = subset(ests_s, select = -c(spp_site_yr))
head(ests_s)
##

#subset by spp
ed <- ests_s[which(ests_s$spp=="ED"), ]
po <- ests_s[which(ests_s$spp=="PO"), ]
hs <- ests_s[which(ests_s$spp=="HS"), ]
ov <- ests_s[which(ests_s$spp=="OV"), ]
gr <- ests_s[which(ests_s$spp=="GR"), ]
et <- ests_s[which(ests_s$spp=="ET"), ]
af <- ests_s[which(ests_s$spp=="AF"), ]
br <- ests_s[which(ests_s$spp=="BR"), ]
aa <- ests_s[which(ests_s$spp=="AA"), ]


#no sci notation
options(scipen = 999)
options(na.action = "na.omit")

##BL plot##
tiff(filename = "plots/Year_DM_wSilhouettes.tiff", width = 12, height = 12, units = 'in', res = 600, compression = 'lzw')

##layout
layout_mat <- matrix(c(1:16), nrow = 4, ncol = 4, byrow = TRUE)
my_lay <- layout(mat = layout_mat, 
                 heights = c(2,2,2,0.3),
                 widths = c(0.3,2,2,2),respect = TRUE)      

#empty plot for where the labels will go
par(mar=c(0.4,0.4,0.4,0.4))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

par(mar=c(2,2.5,2.5,0.4))

###############Ancylus fluviatilis
af[is.na(af)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(af$DM_est + af$DM_SE)), xlim=c(2000,2020))
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=af$yr[af$site=="Auba"]-0.15, y=af$DM_est[af$site=="Auba"], pch=21, bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=af$yr[af$site=="Bieb"]-0.05, y=af$DM_est[af$site=="Bieb"], pch=22, bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=af$yr[af$site=="O3"]+0.05, y=af$DM_est[af$site=="O3"], pch=23, bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=af$yr[af$site=="W1"]+0.15, y=af$DM_est[af$site=="W1"], pch=24, bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
points(x=af$yr[af$site=="Auba"]-0.15, y=af$DM_est[af$site=="Auba"], type="l", bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=af$yr[af$site=="Bieb"]-0.05, y=af$DM_est[af$site=="Bieb"], type="l", bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=af$yr[af$site=="O3"]+0.05, y=af$DM_est[af$site=="O3"], type="l", bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=af$yr[af$site=="W1"]+0.15, y=af$DM_est[af$site=="W1"], type="l", bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
arrows(af$yr[af$site=="Auba"]-0.15, af$DM_est[af$site=="Auba"]-af$DM_SE[af$site=="Auba"], af$yr[af$site=="Auba"]-0.15, af$DM_est[af$site=="Auba"]+af$DM_SE[af$site=="Auba"],col=alpha(1,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(af$yr[af$site=="Bieb"]-0.05, af$DM_est[af$site=="Bieb"]-af$DM_SE[af$site=="Bieb"], af$yr[af$site=="Bieb"]-0.05, af$DM_est[af$site=="Bieb"]+af$DM_SE[af$site=="Bieb"],col=alpha(2,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(af$yr[af$site=="O3"]+0.05, af$DM_est[af$site=="O3"]-af$DM_SE[af$site=="O3"], af$yr[af$site=="O3"]+0.05, af$DM_est[af$site=="O3"]+af$DM_SE[af$site=="O3"],col=alpha(3,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(af$yr[af$site=="W1"]+0.15, af$DM_est[af$site=="W1"]-af$DM_SE[af$site=="W1"], af$yr[af$site=="W1"]+0.15, af$DM_est[af$site=="W1"]+af$DM_SE[af$site=="W1"],col=alpha(4,0.2), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/ancylus_fluviatilis.svg",
  x = 2002,
  y = min(af$DM_est-af$DM_SE-1) + 1.3,
  width = 2.5,
  height = NULL
)
title(main="a. Ancylus fluviatilis",bty="n",cex.main=2)

#####Aphelocheirus aestivalis
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(aa$DM_est + aa$DM_SE)), xlim=c(2000,2020))
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=aa$yr[aa$site=="Auba"]-0.15, y=aa$DM_est[aa$site=="Auba"], pch=21, bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Bieb"]-0.05, y=aa$DM_est[aa$site=="Bieb"], pch=22, bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="O3"]+0.05, y=aa$DM_est[aa$site=="O3"], pch=23, bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="W1"]+0.15, y=aa$DM_est[aa$site=="W1"], pch=24, bg=alpha(4,0.8), col=alpha(4,0.8), lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Auba"]-0.15, y=aa$DM_est[aa$site=="Auba"], type="l", bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Bieb"]-0.05, y=aa$DM_est[aa$site=="Bieb"], type="l", bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="O3"]+0.05, y=aa$DM_est[aa$site=="O3"], type="l", bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="W1"]+0.15, y=aa$DM_est[aa$site=="W1"], type="l", bg=alpha(4,0.8), col=alpha(4,0.8), lwd=2,cex=2.5)
arrows(aa$yr[aa$site=="Auba"]-0.15, aa$DM_est[aa$site=="Auba"]-aa$DM_SE[aa$site=="Auba"], aa$yr[aa$site=="Auba"]-0.15, aa$DM_est[aa$site=="Auba"]+aa$DM_SE[aa$site=="Auba"],col=alpha(1,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="Bieb"]-0.05, aa$DM_est[aa$site=="Bieb"]-aa$DM_SE[aa$site=="Bieb"], aa$yr[aa$site=="Bieb"]-0.05, aa$DM_est[aa$site=="Bieb"]+aa$DM_SE[aa$site=="Bieb"],col=alpha(2,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="O3"]+0.05, aa$DM_est[aa$site=="O3"]-aa$DM_SE[aa$site=="O3"], aa$yr[aa$site=="O3"]+0.05, aa$DM_est[aa$site=="O3"]+aa$DM_SE[aa$site=="O3"],col=alpha(3,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="W1"]+0.15, aa$DM_est[aa$site=="W1"]-aa$DM_SE[aa$site=="W1"], aa$yr[aa$site=="W1"]+0.15, aa$DM_est[aa$site=="W1"]+aa$DM_SE[aa$site=="W1"],col=alpha(4,0.8), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/aphelocheirus_aestivalis.svg",
  x = 2002,
  y = min(aa$DM_est - aa$DM_SE - 1) + 2,
  width = 2.5,
  height = NULL
)
title("b. Aphelocheirus aestivalis",bty="n",cex.main=2)

###############Baetis rhodani  
br[is.na(br)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(br$DM_est + br$DM_SE)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=br$yr[br$site=="Auba"]-0.15, y=br$DM_est[br$site=="Auba"], pch=21, bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=br$yr[br$site=="Bieb"]-0.05, y=br$DM_est[br$site=="Bieb"], pch=22, bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=br$yr[br$site=="O3"]+0.05, y=br$DM_est[br$site=="O3"], pch=23, bg=alpha(3,0.8), col=alpha(3,0.8), lwd=2,cex=2.5)
points(x=br$yr[br$site=="W1"]+0.15, y=br$DM_est[br$site=="W1"], pch=24, bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
points(x=br$yr[br$site=="Auba"]-0.15, y=br$DM_est[br$site=="Auba"], type="l", bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=br$yr[br$site=="Bieb"]-0.05, y=br$DM_est[br$site=="Bieb"], type="l", bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=br$yr[br$site=="O3"]+0.05, y=br$DM_est[br$site=="O3"], type="l", bg=alpha(3,0.8), col=alpha(3,0.8), lwd=2,cex=2.5)
points(x=br$yr[br$site=="W1"]+0.15, y=br$DM_est[br$site=="W1"], type="l", bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
arrows(br$yr[br$site=="Auba"]-0.15, br$DM_est[br$site=="Auba"]-br$DM_SE[br$site=="Auba"], br$yr[br$site=="Auba"]-0.15, br$DM_est[br$site=="Auba"]+br$DM_SE[br$site=="Auba"],col=alpha(1,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="Bieb"]-0.05, br$DM_est[br$site=="Bieb"]-br$DM_SE[br$site=="Bieb"], br$yr[br$site=="Bieb"]-0.05, br$DM_est[br$site=="Bieb"]+br$DM_SE[br$site=="Bieb"],col=alpha(2,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="O3"]+0.05, br$DM_est[br$site=="O3"]-br$DM_SE[br$site=="O3"], br$yr[br$site=="O3"]+0.05, br$DM_est[br$site=="O3"]+br$DM_SE[br$site=="O3"],col=alpha(3,0.8), lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="W1"]+0.15, br$DM_est[br$site=="W1"]-br$DM_SE[br$site=="W1"], br$yr[br$site=="W1"]+0.15, br$DM_est[br$site=="W1"]+br$DM_SE[br$site=="W1"],col=alpha(4,0.2), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/baetis_rhodani.svg",
  x = 2002,
  y = min(br$DM_est-br$DM_SE-1)+ 4,
  width = 4,
  height = NULL
)
title(main="c. Baetis rhodani",bty="n",cex.main=2)

#empty plot for where the labels will go
par(mar=c(0.4,0.4,0.4,0.4))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab="Body size (mg)", line=-2,cex.lab=2)
par(mar=c(2,2.5,2.5,0.4))
    
###############Eiseniella tetraedra 
et[is.na(et)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(et$DM_est + et$DM_SE)), xlim=c(2000,2020))
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=et$yr[et$site=="Auba"]-0.15, y=et$DM_est[et$site=="Auba"], pch=21, bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=et$yr[et$site=="Bieb"]-0.05, y=et$DM_est[et$site=="Bieb"], pch=22, bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=et$yr[et$site=="O3"]+0.05, y=et$DM_est[et$site=="O3"], pch=23, bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=et$yr[et$site=="W1"]+0.15, y=et$DM_est[et$site=="W1"], pch=24, bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
points(x=et$yr[et$site=="Auba"]-0.15, y=et$DM_est[et$site=="Auba"], type="l", bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=et$yr[et$site=="Bieb"]-0.05, y=et$DM_est[et$site=="Bieb"], type="l", bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=et$yr[et$site=="O3"]+0.05, y=et$DM_est[et$site=="O3"], type="l", bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=et$yr[et$site=="W1"]+0.15, y=et$DM_est[et$site=="W1"], type="l", bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
arrows(et$yr[et$site=="Auba"]-0.15, et$DM_est[et$site=="Auba"]-et$DM_SE[et$site=="Auba"], et$yr[et$site=="Auba"]-0.15, et$DM_est[et$site=="Auba"]+et$DM_SE[et$site=="Auba"],col=alpha(1,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="Bieb"]-0.05, et$DM_est[et$site=="Bieb"]-et$DM_SE[et$site=="Bieb"], et$yr[et$site=="Bieb"]-0.05, et$DM_est[et$site=="Bieb"]+et$DM_SE[et$site=="Bieb"],col=alpha(2,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="O3"]+0.05, et$DM_est[et$site=="O3"]-et$DM_SE[et$site=="O3"], et$yr[et$site=="O3"]+0.05, et$DM_est[et$site=="O3"]+et$DM_SE[et$site=="O3"],col=alpha(3,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="W1"]+0.15, et$DM_est[et$site=="W1"]-et$DM_SE[et$site=="W1"], et$yr[et$site=="W1"]+0.15, et$DM_est[et$site=="W1"]+et$DM_SE[et$site=="W1"],col=alpha(4,0.2), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/eiseniella_tetraedra.svg",
  x = 2002,
  y = min(et$DM_est-et$DM_SE-1) + 1,
  width = 3,
  height = NULL
)
title(main="d. Eiseniella tetraedra",bty="n",cex.main=2)
legend("topleft", legend=c("Aubach","Bieber","Kinzig O3","Kinzig W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)

###############Ephemera danica 
ed[is.na(ed)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(ed$DM_est + ed$DM_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ed$yr[ed$site=="Auba"]-0.15, y=ed$DM_est[ed$site=="Auba"], pch=21, bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Bieb"]-0.05, y=ed$DM_est[ed$site=="Bieb"], pch=22, bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="O3"]+0.05, y=ed$DM_est[ed$site=="O3"], pch=23, bg=alpha(3,0.8), col=alpha(3,0.8), lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="W1"]+0.15, y=ed$DM_est[ed$site=="W1"], pch=24, bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Auba"]-0.15, y=ed$DM_est[ed$site=="Auba"], type="l", bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Bieb"]-0.05, y=ed$DM_est[ed$site=="Bieb"], type="l", bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="O3"]+0.05, y=ed$DM_est[ed$site=="O3"], type="l", bg=alpha(3,0.8), col=alpha(3,0.8), lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="W1"]+0.15, y=ed$DM_est[ed$site=="W1"], type="l", bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
arrows(ed$yr[ed$site=="Auba"]-0.15, ed$DM_est[ed$site=="Auba"]-ed$DM_SE[ed$site=="Auba"], ed$yr[ed$site=="Auba"]-0.15, ed$DM_est[ed$site=="Auba"]+ed$DM_SE[ed$site=="Auba"],col=alpha(1,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="Bieb"]-0.05, ed$DM_est[ed$site=="Bieb"]-ed$DM_SE[ed$site=="Bieb"], ed$yr[ed$site=="Bieb"]-0.05, ed$DM_est[ed$site=="Bieb"]+ed$DM_SE[ed$site=="Bieb"],col=alpha(2,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="O3"]+0.05, ed$DM_est[ed$site=="O3"]-ed$DM_SE[ed$site=="O3"], ed$yr[ed$site=="O3"]+0.05, ed$DM_est[ed$site=="O3"]+ed$DM_SE[ed$site=="O3"],col=alpha(3,0.8), lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="W1"]+0.15, ed$DM_est[ed$site=="W1"]-ed$DM_SE[ed$site=="W1"], ed$yr[ed$site=="W1"]+0.15, ed$DM_est[ed$site=="W1"]+ed$DM_SE[ed$site=="W1"],col=alpha(4,0.2), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/Ephemera_danica.svg",
  x = 2002,
  y = min(ed$DM_est-ed$DM_SE-1) + 2.75,
  width = 3.5,
  height = NULL
)
title(main="e. Ephemera danica ",bty="n",cex.main=2)

###############Gammarus roeselii
gr[is.na(gr)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(gr$DM_est + gr$DM_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=gr$yr[gr$site=="Auba"]-0.15, y=gr$DM_est[gr$site=="Auba"], pch=21, bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Bieb"]-0.05, y=gr$DM_est[gr$site=="Bieb"], pch=22, bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="O3"]+0.05, y=gr$DM_est[gr$site=="O3"], pch=23, bg=alpha(3,0.8), col=alpha(3,0.8), lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="W1"]+0.15, y=gr$DM_est[gr$site=="W1"], pch=24, bg=alpha(4,0.8), col=alpha(4,0.8), lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Auba"]-0.15, y=gr$DM_est[gr$site=="Auba"], type="l", bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Bieb"]-0.05, y=gr$DM_est[gr$site=="Bieb"], type="l", bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="O3"]+0.05, y=gr$DM_est[gr$site=="O3"], type="l", bg=alpha(3,0.8), col=alpha(3,0.8), lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="W1"]+0.15, y=gr$DM_est[gr$site=="W1"], type="l", bg=alpha(4,0.8), col=alpha(4,0.8), lwd=2,cex=2.5)
arrows(gr$yr[gr$site=="Auba"]-0.15, gr$DM_est[gr$site=="Auba"]-gr$DM_SE[gr$site=="Auba"], gr$yr[gr$site=="Auba"]-0.15, gr$DM_est[gr$site=="Auba"]+gr$DM_SE[gr$site=="Auba"],col=alpha(1,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(gr$yr[gr$site=="Bieb"]-0.05, gr$DM_est[gr$site=="Bieb"]-gr$DM_SE[gr$site=="Bieb"], gr$yr[gr$site=="Bieb"]-0.05, gr$DM_est[gr$site=="Bieb"]+gr$DM_SE[gr$site=="Bieb"],col=alpha(2,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(gr$yr[gr$site=="O3"]+0.05, gr$DM_est[gr$site=="O3"]-gr$DM_SE[gr$site=="O3"], gr$yr[gr$site=="O3"]+0.05, gr$DM_est[gr$site=="O3"]+gr$DM_SE[gr$site=="O3"],col=alpha(3,0.8), lwd=2,length=0.05, angle=90, code=3)
arrows(gr$yr[gr$site=="W1"]+0.15, gr$DM_est[gr$site=="W1"]-gr$DM_SE[gr$site=="W1"], gr$yr[gr$site=="W1"]+0.15, gr$DM_est[gr$site=="W1"]+gr$DM_SE[gr$site=="W1"],col=alpha(4,0.8), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/Gammarus_roselii.svg",
  x = 2002,
  y = min(gr$DM_est-gr$DM_SE-1) + 4,
  width = 3,
  height = NULL
)
title(main="f. Gammarus roeselii",bty="n",cex.main=2)

#empty plot for where the labels will go
par(mar=c(0.4,0.4,0.4,0.4))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
par(mar=c(2,2.5,2.5,0.4))

###############Hydropsyche siltalai 
hs[is.na(hs)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(hs$DM_est + hs$DM_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=hs$yr[hs$site=="Auba"]-0.15, y=hs$DM_est[hs$site=="Auba"], pch=21, bg=alpha(1,0.8), col=alpha(1,0.8), lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Bieb"]-0.05, y=hs$DM_est[hs$site=="Bieb"], pch=22, bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="O3"]+0.05, y=hs$DM_est[hs$site=="O3"], pch=23, bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="W1"]+0.15, y=hs$DM_est[hs$site=="W1"], pch=24, bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Auba"]-0.15, y=hs$DM_est[hs$site=="Auba"], type="l", bg=alpha(1,0.8), col=alpha(1,0.8), lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Bieb"]-0.05, y=hs$DM_est[hs$site=="Bieb"], type="l", bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="O3"]+0.05, y=hs$DM_est[hs$site=="O3"], type="l", bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="W1"]+0.15, y=hs$DM_est[hs$site=="W1"], type="l", bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
arrows(hs$yr[hs$site=="Auba"]-0.15, hs$DM_est[hs$site=="Auba"]-hs$DM_SE[hs$site=="Auba"], hs$yr[hs$site=="Auba"]-0.15, hs$DM_est[hs$site=="Auba"]+hs$DM_SE[hs$site=="Auba"],col=alpha(1,0.8), lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="Bieb"]-0.05, hs$DM_est[hs$site=="Bieb"]-hs$DM_SE[hs$site=="Bieb"], hs$yr[hs$site=="Bieb"]-0.05, hs$DM_est[hs$site=="Bieb"]+hs$DM_SE[hs$site=="Bieb"],col=alpha(2,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="O3"]+0.05, hs$DM_est[hs$site=="O3"]-hs$DM_SE[hs$site=="O3"], hs$yr[hs$site=="O3"]+0.05, hs$DM_est[hs$site=="O3"]+hs$DM_SE[hs$site=="O3"],col=alpha(3,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="W1"]+0.15, hs$DM_est[hs$site=="W1"]-hs$DM_SE[hs$site=="W1"], hs$yr[hs$site=="W1"]+0.15, hs$DM_est[hs$site=="W1"]+hs$DM_SE[hs$site=="W1"],col=alpha(4,0.2), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/hydropsyche_siltalai.svg",
  x = 2002,
  y = min(hs$DM_est-hs$DM_SE-1) + 1,
  width = 4,
  height = NULL
)
title(main="g. Hydropsyche siltalai",bty="n",cex.main=2)

###############Orectochilus villosus
ov[is.na(ov)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(ov$DM_est + ov$DM_SE+1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ov$yr[ov$site=="Auba"]-0.15, y=ov$DM_est[ov$site=="Auba"], pch=21, bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Bieb"]-0.05, y=ov$DM_est[ov$site=="Bieb"], pch=22, bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="O3"]+0.05, y=ov$DM_est[ov$site=="O3"], pch=23, bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="W1"]+0.15, y=ov$DM_est[ov$site=="W1"], pch=24, bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Auba"]-0.15, y=ov$DM_est[ov$site=="Auba"], type="l", bg=alpha(1,0.2), col=alpha(1,0.2), lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Bieb"]-0.05, y=ov$DM_est[ov$site=="Bieb"], type="l", bg=alpha(2,0.2), col=alpha(2,0.2), lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="O3"]+0.05, y=ov$DM_est[ov$site=="O3"], type="l", bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="W1"]+0.15, y=ov$DM_est[ov$site=="W1"], type="l", bg=alpha(4,0.2), col=alpha(4,0.2), lwd=2,cex=2.5)
arrows(ov$yr[ov$site=="Auba"]-0.15, ov$DM_est[ov$site=="Auba"]-ov$DM_SE[ov$site=="Auba"], ov$yr[ov$site=="Auba"]-0.15, ov$DM_est[ov$site=="Auba"]+ov$DM_SE[ov$site=="Auba"],col=alpha(1,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="Bieb"]-0.05, ov$DM_est[ov$site=="Bieb"]-ov$DM_SE[ov$site=="Bieb"], ov$yr[ov$site=="Bieb"]-0.05, ov$DM_est[ov$site=="Bieb"]+ov$DM_SE[ov$site=="Bieb"],col=alpha(2,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="O3"]+0.05, ov$DM_est[ov$site=="O3"]-ov$DM_SE[ov$site=="O3"], ov$yr[ov$site=="O3"]+0.05, ov$DM_est[ov$site=="O3"]+ov$DM_SE[ov$site=="O3"],col=alpha(3,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="W1"]+0.15, ov$DM_est[ov$site=="W1"]-ov$DM_SE[ov$site=="W1"], ov$yr[ov$site=="W1"]+0.15, ov$DM_est[ov$site=="W1"]+ov$DM_SE[ov$site=="W1"],col=alpha(4,0.2), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/Orectochilus_villosus.svg",
  x = 2002,
  y = min(ov$DM_est-ov$DM_SE-1) + 40,
  width = 2.5,
  height = NULL
)
title(main="h. Orectochilus villosus",bty="n",cex.main=2)

###############Prodiamesa olivacea  
po[is.na(po)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(0,max(po$DM_est + po$DM_SE)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=po$yr[po$site=="Auba"]-0.15, y=po$DM_est[po$site=="Auba"], pch=21, bg=alpha(1,0.8), col=alpha(1,0.8), lwd=2,cex=2.5)
points(x=po$yr[po$site=="Bieb"]-0.05, y=po$DM_est[po$site=="Bieb"], pch=22, bg=alpha(2,0.8), col=alpha(2,0.8), lwd=2,cex=2.5)
points(x=po$yr[po$site=="O3"]+0.05, y=po$DM_est[po$site=="O3"], pch=23, bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=po$yr[po$site=="W1"]+0.15, y=po$DM_est[po$site=="W1"], pch=24, bg=alpha(4,0.8), col=alpha(4,0.8), lwd=2,cex=2.5)
points(x=po$yr[po$site=="Auba"]-0.15, y=po$DM_est[po$site=="Auba"], type="l", bg=alpha(1,0.8), col=alpha(1,0.8), lwd=2,cex=2.5)
points(x=po$yr[po$site=="Bieb"]-0.05, y=po$DM_est[po$site=="Bieb"], type="l", bg=alpha(2,0.8), col=alpha(2,0.8), lwd=2,cex=2.5)
points(x=po$yr[po$site=="O3"]+0.05, y=po$DM_est[po$site=="O3"], type="l", bg=alpha(3,0.2), col=alpha(3,0.2), lwd=2,cex=2.5)
points(x=po$yr[po$site=="W1"]+0.15, y=po$DM_est[po$site=="W1"], type="l", bg=alpha(4,0.8), col=alpha(4,0.8), lwd=2,cex=2.5)
arrows(po$yr[po$site=="Auba"]-0.15, po$DM_est[po$site=="Auba"]-po$DM_SE[po$site=="Auba"], po$yr[po$site=="Auba"]-0.15, po$DM_est[po$site=="Auba"]+po$DM_SE[po$site=="Auba"],col=alpha(1,0.8), lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="Bieb"]-0.05, po$DM_est[po$site=="Bieb"]-po$DM_SE[po$site=="Bieb"], po$yr[po$site=="Bieb"]-0.05, po$DM_est[po$site=="Bieb"]+po$DM_SE[po$site=="Bieb"],col=alpha(2,0.8), lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="O3"]+0.05, po$DM_est[po$site=="O3"]-po$DM_SE[po$site=="O3"], po$yr[po$site=="O3"]+0.05, po$DM_est[po$site=="O3"]+po$DM_SE[po$site=="O3"],col=alpha(3,0.2), lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="W1"]+0.15, po$DM_est[po$site=="W1"]-po$DM_SE[po$site=="W1"], po$yr[po$site=="W1"]+0.15, po$DM_est[po$site=="W1"]+po$DM_SE[po$site=="W1"],col=alpha(4,0.8), lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/Prodiamesa_olivacea.svg",
  x = 2002,
  y = min(po$DM_est-po$DM_SE-1) + 2,
  width = 3.5,
  height = NULL
)
title(main="i. Prodiamesa olivacea ",bty="n",cex.main=2)

#empty plot for where the labels will go
par(mar=c(0.4,0.4,0.4,0.4))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(xlab="Year", line=-2,cex.lab=2)
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
 

dev.off()
##
##


