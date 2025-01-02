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

#calculate BL estimates for each spp, site, and year
intra$spp_site_yr <- paste(intra$SppCode, intra$SiteShort,intra$yr)
ests_s <- NULL
for(i in unique(intra$spp_site_yr)){
  tryCatch({
  sub <- intra[intra$spp_site_yr == i, ]
	sub<-sub[complete.cases(sub[, "HW"]),]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
  ests.i <- coef(summary(lm(HW ~ 1 + sDOY + Ldens, data = sub)))[1,1:2]
  ests.i <- data.frame(spp_site_yr = i, t(ests.i))
  ests_s <- rbind(ests_s, ests.i) ; rm(ests.i, sub)
  }, error=function(e){cat(unique(sub$spp_site_yr),conditionMessage(e), "\n")})
} ; rm(i)
ests_s

ests_s[c('spp','site', 'yr')] <- str_split_fixed(ests_s$spp_site_yr, ' ', 3)
colnames(ests_s)[2] ="HW_est"
colnames(ests_s)[3] ="HW_SE"
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

##HW plot##
tiff(filename = "plots/Year_HW_wSilhouette.tiff", width = 8, height = 12, units = 'in', res = 600, compression = 'lzw')
par(mar=c(3.6,5,4,0.5),mfrow=c(4,2))

#####Aphelocheirus aestivalis

plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(aa$HW_est-aa$HW_SE-.1),max(aa$HW_est + aa$HW_SE+.1)), xlim=c(2000,2020))
title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=aa$yr[aa$site=="Auba"]-0.15, y=aa$HW_est[aa$site=="Auba"], pch=21, bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Bieb"]-0.05, y=aa$HW_est[aa$site=="Bieb"], pch=22, bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="O3"]+0.05, y=aa$HW_est[aa$site=="O3"], pch=23, bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="W1"]+0.15, y=aa$HW_est[aa$site=="W1"], pch=24, bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Auba"]-0.15, y=aa$HW_est[aa$site=="Auba"], type="l", bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Bieb"]-0.05, y=aa$HW_est[aa$site=="Bieb"], type="l", bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="O3"]+0.05, y=aa$HW_est[aa$site=="O3"], type="l", bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="W1"]+0.15, y=aa$HW_est[aa$site=="W1"], type="l", bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
arrows(aa$yr[aa$site=="Auba"]-0.15, aa$HW_est[aa$site=="Auba"]-aa$HW_SE[aa$site=="Auba"], aa$yr[aa$site=="Auba"]-0.15, aa$HW_est[aa$site=="Auba"]+aa$HW_SE[aa$site=="Auba"],col=alpha(1,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="Bieb"]-0.05, aa$HW_est[aa$site=="Bieb"]-aa$HW_SE[aa$site=="Bieb"], aa$yr[aa$site=="Bieb"]-0.05, aa$HW_est[aa$site=="Bieb"]+aa$HW_SE[aa$site=="Bieb"],col=alpha(2,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="O3"]+0.05, aa$HW_est[aa$site=="O3"]-aa$HW_SE[aa$site=="O3"], aa$yr[aa$site=="O3"]+0.05, aa$HW_est[aa$site=="O3"]+aa$HW_SE[aa$site=="O3"],col=alpha(3,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(aa$yr[aa$site=="W1"]+0.15, aa$HW_est[aa$site=="W1"]-aa$HW_SE[aa$site=="W1"], aa$yr[aa$site=="W1"]+0.15, aa$HW_est[aa$site=="W1"]+aa$HW_SE[aa$site=="W1"],col=alpha(4,0.2),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/aphelocheirus_aestivalis.svg",
  x = 2001,
  y = 1,
  width = 1.5,
  height = NULL
)
title("a. Aphelocheirus aestivalis",bty="n",cex.main=2)
 
###############Baetis rhodani  
br[is.na(br)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(br$HW_est-br$HW_SE-.1),max(br$HW_est + br$HW_SE+.1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=br$yr[br$site=="Auba"]-0.15, y=br$HW_est[br$site=="Auba"], pch=21, bg=alpha(1,0.8),col=alpha(1,0.8),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Bieb"]-0.05, y=br$HW_est[br$site=="Bieb"], pch=22, bg=alpha(2,0.8),col=alpha(2,0.8),lwd=2,cex=2.5)
points(x=br$yr[br$site=="O3"]+0.05, y=br$HW_est[br$site=="O3"], pch=23, bg=alpha(3,0.8),col=alpha(3,0.8),lwd=2,cex=2.5)
points(x=br$yr[br$site=="W1"]+0.15, y=br$HW_est[br$site=="W1"], pch=24, bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Auba"]-0.15, y=br$HW_est[br$site=="Auba"], type="l", bg=alpha(1,0.8),col=alpha(1,0.8),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Bieb"]-0.05, y=br$HW_est[br$site=="Bieb"], type="l", bg=alpha(2,0.8),col=alpha(2,0.8),lwd=2,cex=2.5)
points(x=br$yr[br$site=="O3"]+0.05, y=br$HW_est[br$site=="O3"], type="l", bg=alpha(3,0.8),col=alpha(3,0.8),lwd=2,cex=2.5)
points(x=br$yr[br$site=="W1"]+0.15, y=br$HW_est[br$site=="W1"], type="l", bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
arrows(br$yr[br$site=="Auba"]-0.15, br$HW_est[br$site=="Auba"]-br$HW_SE[br$site=="Auba"], br$yr[br$site=="Auba"]-0.15, br$HW_est[br$site=="Auba"]+br$HW_SE[br$site=="Auba"],col=alpha(1,0.8),lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="Bieb"]-0.05, br$HW_est[br$site=="Bieb"]-br$HW_SE[br$site=="Bieb"], br$yr[br$site=="Bieb"]-0.05, br$HW_est[br$site=="Bieb"]+br$HW_SE[br$site=="Bieb"],col=alpha(2,0.8),lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="O3"]+0.05, br$HW_est[br$site=="O3"]-br$HW_SE[br$site=="O3"], br$yr[br$site=="O3"]+0.05, br$HW_est[br$site=="O3"]+br$HW_SE[br$site=="O3"],col=alpha(3,0.8),lwd=2,length=0.05, angle=90, code=3)
arrows(br$yr[br$site=="W1"]+0.15, br$HW_est[br$site=="W1"]-br$HW_SE[br$site=="W1"], br$yr[br$site=="W1"]+0.15, br$HW_est[br$site=="W1"]+br$HW_SE[br$site=="W1"],col=alpha(4,0.2),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/baetis_rhodani.svg",
  x = 2001,
  y = 1.3,
  width = 2.75,
  height = NULL
)
title(main="b. Baetis rhodani",bty="n",cex.main=2)
            
###############Eiseniella tetraedra 
et[is.na(et)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(et$HW_est-et$HW_SE-.1),max(et$HW_est + et$HW_SE+.1)), xlim=c(2000,2020))
title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=et$yr[et$site=="Auba"]-0.15, y=et$HW_est[et$site=="Auba"], pch=21, bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Bieb"]-0.05, y=et$HW_est[et$site=="Bieb"], pch=22, bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=et$yr[et$site=="O3"]+0.05, y=et$HW_est[et$site=="O3"], pch=23, bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=et$yr[et$site=="W1"]+0.15, y=et$HW_est[et$site=="W1"], pch=24, bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Auba"]-0.15, y=et$HW_est[et$site=="Auba"], type="l", bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Bieb"]-0.05, y=et$HW_est[et$site=="Bieb"], type="l", bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=et$yr[et$site=="O3"]+0.05, y=et$HW_est[et$site=="O3"], type="l", bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=et$yr[et$site=="W1"]+0.15, y=et$HW_est[et$site=="W1"], type="l", bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
arrows(et$yr[et$site=="Auba"]-0.15, et$HW_est[et$site=="Auba"]-et$HW_SE[et$site=="Auba"], et$yr[et$site=="Auba"]-0.15, et$HW_est[et$site=="Auba"]+et$HW_SE[et$site=="Auba"],col=alpha(1,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="Bieb"]-0.05, et$HW_est[et$site=="Bieb"]-et$HW_SE[et$site=="Bieb"], et$yr[et$site=="Bieb"]-0.05, et$HW_est[et$site=="Bieb"]+et$HW_SE[et$site=="Bieb"],col=alpha(2,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="O3"]+0.05, et$HW_est[et$site=="O3"]-et$HW_SE[et$site=="O3"], et$yr[et$site=="O3"]+0.05, et$HW_est[et$site=="O3"]+et$HW_SE[et$site=="O3"],col=alpha(3,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(et$yr[et$site=="W1"]+0.15, et$HW_est[et$site=="W1"]-et$HW_SE[et$site=="W1"], et$yr[et$site=="W1"]+0.15, et$HW_est[et$site=="W1"]+et$HW_SE[et$site=="W1"],col=alpha(4,0.2),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/eiseniella_tetraedra.svg",
  x = 2001,
  y = 1,
  width = 2.5,
  height = NULL
)
title(main="c. Eiseniella tetraedra",bty="n",cex.main=2)

###############Ephemera danica 
ed[is.na(ed)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ed$HW_est-ed$HW_SE-.1),max(ed$HW_est + ed$HW_SE+.1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ed$yr[ed$site=="Auba"]-0.15, y=ed$HW_est[ed$site=="Auba"], pch=21, bg=alpha(1,0.8),col=alpha(1,0.8),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Bieb"]-0.05, y=ed$HW_est[ed$site=="Bieb"], pch=22, bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="O3"]+0.05, y=ed$HW_est[ed$site=="O3"], pch=23, bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="W1"]+0.15, y=ed$HW_est[ed$site=="W1"], pch=24, bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Auba"]-0.15, y=ed$HW_est[ed$site=="Auba"], type="l", bg=alpha(1,0.8),col=alpha(1,0.8),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Bieb"]-0.05, y=ed$HW_est[ed$site=="Bieb"], type="l", bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="O3"]+0.05, y=ed$HW_est[ed$site=="O3"], type="l", bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="W1"]+0.15, y=ed$HW_est[ed$site=="W1"], type="l", bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
arrows(ed$yr[ed$site=="Auba"]-0.15, ed$HW_est[ed$site=="Auba"]-ed$HW_SE[ed$site=="Auba"], ed$yr[ed$site=="Auba"]-0.15, ed$HW_est[ed$site=="Auba"]+ed$HW_SE[ed$site=="Auba"],col=alpha(1,0.8),lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="Bieb"]-0.05, ed$HW_est[ed$site=="Bieb"]-ed$HW_SE[ed$site=="Bieb"], ed$yr[ed$site=="Bieb"]-0.05, ed$HW_est[ed$site=="Bieb"]+ed$HW_SE[ed$site=="Bieb"],col=alpha(2,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="O3"]+0.05, ed$HW_est[ed$site=="O3"]-ed$HW_SE[ed$site=="O3"], ed$yr[ed$site=="O3"]+0.05, ed$HW_est[ed$site=="O3"]+ed$HW_SE[ed$site=="O3"],col=alpha(3,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(ed$yr[ed$site=="W1"]+0.15, ed$HW_est[ed$site=="W1"]-ed$HW_SE[ed$site=="W1"], ed$yr[ed$site=="W1"]+0.15, ed$HW_est[ed$site=="W1"]+ed$HW_SE[ed$site=="W1"],col=alpha(4,0.2),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/Ephemera_danica.svg",
  x = 2001,
  y = 0.6,
  width = 2.5,
  height = NULL
)
title(main="d. Ephemera danica ",bty="n",cex.main=2)
 
###############Hydropsyche siltalai 
hs[is.na(hs)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(hs$HW_est-hs$HW_SE-.1),max(hs$HW_est + hs$HW_SE+.1)), xlim=c(2000,2020))
title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=hs$yr[hs$site=="Auba"]-0.15, y=hs$HW_est[hs$site=="Auba"], pch=21, bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Bieb"]-0.05, y=hs$HW_est[hs$site=="Bieb"], pch=22, bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="O3"]+0.05, y=hs$HW_est[hs$site=="O3"], pch=23, bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="W1"]+0.15, y=hs$HW_est[hs$site=="W1"], pch=24, bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Auba"]-0.15, y=hs$HW_est[hs$site=="Auba"], type="l", bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Bieb"]-0.05, y=hs$HW_est[hs$site=="Bieb"], type="l", bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="O3"]+0.05, y=hs$HW_est[hs$site=="O3"], type="l", bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="W1"]+0.15, y=hs$HW_est[hs$site=="W1"], type="l", bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
arrows(hs$yr[hs$site=="Auba"]-0.15, hs$HW_est[hs$site=="Auba"]-hs$HW_SE[hs$site=="Auba"], hs$yr[hs$site=="Auba"]-0.15, hs$HW_est[hs$site=="Auba"]+hs$HW_SE[hs$site=="Auba"],col=alpha(1,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="Bieb"]-0.05, hs$HW_est[hs$site=="Bieb"]-hs$HW_SE[hs$site=="Bieb"], hs$yr[hs$site=="Bieb"]-0.05, hs$HW_est[hs$site=="Bieb"]+hs$HW_SE[hs$site=="Bieb"],col=alpha(2,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="O3"]+0.05, hs$HW_est[hs$site=="O3"]-hs$HW_SE[hs$site=="O3"], hs$yr[hs$site=="O3"]+0.05, hs$HW_est[hs$site=="O3"]+hs$HW_SE[hs$site=="O3"],col=alpha(3,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(hs$yr[hs$site=="W1"]+0.15, hs$HW_est[hs$site=="W1"]-hs$HW_SE[hs$site=="W1"], hs$yr[hs$site=="W1"]+0.15, hs$HW_est[hs$site=="W1"]+hs$HW_SE[hs$site=="W1"],col=alpha(4,0.2),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/hydropsyche_siltalai.svg",
  x = 2001,
  y = 0.6,
  width = 3,
  height = NULL
)
title(main="e. Hydropsyche siltalai",bty="n",cex.main=2)
 
###############Orectochilus villosus
ov[is.na(ov)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ov$HW_est-ov$HW_SE-.1),max(ov$HW_est + ov$HW_SE+.1)), xlim=c(2000,2020))
#title(ylab="Body length (mm)", line=2.7,cex.lab=1.5)
title(xlab="Year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ov$yr[ov$site=="Auba"]-0.15, y=ov$HW_est[ov$site=="Auba"], pch=21, bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Bieb"]-0.05, y=ov$HW_est[ov$site=="Bieb"], pch=22, bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="O3"]+0.05, y=ov$HW_est[ov$site=="O3"], pch=23, bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="W1"]+0.15, y=ov$HW_est[ov$site=="W1"], pch=24, bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Auba"]-0.15, y=ov$HW_est[ov$site=="Auba"], type="l", bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Bieb"]-0.05, y=ov$HW_est[ov$site=="Bieb"], type="l", bg=alpha(2,0.2),col=alpha(2,0.2),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="O3"]+0.05, y=ov$HW_est[ov$site=="O3"], type="l", bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="W1"]+0.15, y=ov$HW_est[ov$site=="W1"], type="l", bg=alpha(4,0.2),col=alpha(4,0.2),lwd=2,cex=2.5)
arrows(ov$yr[ov$site=="Auba"]-0.15, ov$HW_est[ov$site=="Auba"]-ov$HW_SE[ov$site=="Auba"], ov$yr[ov$site=="Auba"]-0.15, ov$HW_est[ov$site=="Auba"]+ov$HW_SE[ov$site=="Auba"],col=alpha(1,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="Bieb"]-0.05, ov$HW_est[ov$site=="Bieb"]-ov$HW_SE[ov$site=="Bieb"], ov$yr[ov$site=="Bieb"]-0.05, ov$HW_est[ov$site=="Bieb"]+ov$HW_SE[ov$site=="Bieb"],col=alpha(2,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="O3"]+0.05, ov$HW_est[ov$site=="O3"]-ov$HW_SE[ov$site=="O3"], ov$yr[ov$site=="O3"]+0.05, ov$HW_est[ov$site=="O3"]+ov$HW_SE[ov$site=="O3"],col=alpha(3,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(ov$yr[ov$site=="W1"]+0.15, ov$HW_est[ov$site=="W1"]-ov$HW_SE[ov$site=="W1"], ov$yr[ov$site=="W1"]+0.15, ov$HW_est[ov$site=="W1"]+ov$HW_SE[ov$site=="W1"],col=alpha(4,0.2),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/Orectochilus_villosus.svg",
  x = 2001,
  y = 1.75,
  width = 1.5,
  height = NULL
)
title(main="f. Orectochilus villosus",bty="n",cex.main=2)
 
###############Prodiamesa olivacea  
po[is.na(po)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(po$HW_est-po$HW_SE-.1),max(po$HW_est + po$HW_SE+.1)), xlim=c(2000,2020))
title(ylab="Head width (mm)", line=2.7,cex.lab=1.5)
title(xlab="Year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=po$yr[po$site=="Auba"]-0.15, y=po$HW_est[po$site=="Auba"], pch=21, bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Bieb"]-0.05, y=po$HW_est[po$site=="Bieb"], pch=22, bg=alpha(2,0.8),col=alpha(2,0.8),lwd=2,cex=2.5)
points(x=po$yr[po$site=="O3"]+0.05, y=po$HW_est[po$site=="O3"], pch=23, bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=po$yr[po$site=="W1"]+0.15, y=po$HW_est[po$site=="W1"], pch=24, bg=alpha(4,0.8),col=alpha(4,0.8),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Auba"]-0.15, y=po$HW_est[po$site=="Auba"], type="l", bg=alpha(1,0.2),col=alpha(1,0.2),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Bieb"]-0.05, y=po$HW_est[po$site=="Bieb"], type="l", bg=alpha(2,0.8),col=alpha(2,0.8),lwd=2,cex=2.5)
points(x=po$yr[po$site=="O3"]+0.05, y=po$HW_est[po$site=="O3"], type="l", bg=alpha(3,0.2),col=alpha(3,0.2),lwd=2,cex=2.5)
points(x=po$yr[po$site=="W1"]+0.15, y=po$HW_est[po$site=="W1"], type="l", bg=alpha(4,0.8),col=alpha(4,0.8),lwd=2,cex=2.5)
arrows(po$yr[po$site=="Auba"]-0.15, po$HW_est[po$site=="Auba"]-po$HW_SE[po$site=="Auba"], po$yr[po$site=="Auba"]-0.15, po$HW_est[po$site=="Auba"]+po$HW_SE[po$site=="Auba"],col=alpha(1,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="Bieb"]-0.05, po$HW_est[po$site=="Bieb"]-po$HW_SE[po$site=="Bieb"], po$yr[po$site=="Bieb"]-0.05, po$HW_est[po$site=="Bieb"]+po$HW_SE[po$site=="Bieb"],col=alpha(2,0.8),lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="O3"]+0.05, po$HW_est[po$site=="O3"]-po$HW_SE[po$site=="O3"], po$yr[po$site=="O3"]+0.05, po$HW_est[po$site=="O3"]+po$HW_SE[po$site=="O3"],col=alpha(3,0.2),lwd=2,length=0.05, angle=90, code=3)
arrows(po$yr[po$site=="W1"]+0.15, po$HW_est[po$site=="W1"]-po$HW_SE[po$site=="W1"], po$yr[po$site=="W1"]+0.15, po$HW_est[po$site=="W1"]+po$HW_SE[po$site=="W1"],col=alpha(4,0.8),lwd=2,length=0.05, angle=90, code=3)
add_silhouette(
  upload_img = "Silhouettes/Prodiamesa_olivacea.svg",
  x = 2001,
  y = 0.25,
  width = 2.5,
  height = NULL
)
title(main="g. Prodiamesa olivacea ",bty="n",cex.main=2)

plot(1, 1,type= "n",bty ="n",axes=F,frame.plot=F, xaxt='n', ann=FALSE, yaxt='n')
legend("topleft", legend=c("Aubach","Bieber","Kinzig O3","Kinzig W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=3, cex=2.5)

dev.off()
##
##


