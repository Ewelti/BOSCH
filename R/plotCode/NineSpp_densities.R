##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

# load libraries
library(scales)
library(stringr)

# function for adding silhouettes
source("R/Add_silhouette_function.R")

# attach data
intra <- read.csv("RawData/NineSppDensities.csv", header=T)
head(intra)
nrow(intra)
unique(intra$Species)

#fix some problem where these variables are not numeric
intra$BL <- as.numeric(intra$density_per_m2) 


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

#calculate BL estimates for each spp, site, and year
intra$spp_site_yr <- paste(intra$SppCode, intra$SiteShort,intra$yr)
ests_s <- NULL
for(i in unique(intra$spp_site_yr)){
  tryCatch({
  sub <- intra[intra$spp_site_yr == i, ]
	sub<-sub[complete.cases(sub[, "density_per_m2"]),]
  ests.i <- coef(summary(lm(density_per_m2 ~ 1 + sDOY, data = sub)))[1,1:2]
  ests.i <- data.frame(spp_site_yr = i, t(ests.i))
  ests_s <- rbind(ests_s, ests.i) ; rm(ests.i, sub)
  }, error=function(e){cat(unique(sub$spp_site_yr),conditionMessage(e), "\n")})
} ; rm(i)
ests_s

ests_s[c('spp','site', 'yr')] <- str_split_fixed(ests_s$spp_site_yr, ' ', 3)
colnames(ests_s)[2] ="Dens_est"
colnames(ests_s)[3] ="Dens_SE"
ests_s$yr <-as.numeric(ests_s$yr)
ests_s = subset(ests_s, select = -c(spp_site_yr))
ests_s = subset(ests_s, select = -c(Dens_SE))
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

##Density plot##
tiff(filename = "plots/Year_Dens_wSilhouette.tiff", width = 12, height = 12, units = 'in', res = 600, compression = 'lzw')

##layout
layout_mat <- matrix(c(1:16), nrow = 4, ncol = 4, byrow = TRUE)
my_lay <- layout(mat = layout_mat, 
                 heights = c(2,2,2,0.3),
                 widths = c(0.3,2,2,2),respect = TRUE)      
layout.show(16)
par(mar=c(2,2.5,2.5,0.4))

#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

#####Aphelocheirus aestivalis

plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(aa$Dens_est),max(aa$Dens_est)), xlim=c(2000,2020))
#title(ylab="Density (XX)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=aa$yr[aa$site=="Auba"]-0.15, y=aa$Dens_est[aa$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Bieb"]-0.05, y=aa$Dens_est[aa$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="O3"]+0.05, y=aa$Dens_est[aa$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="W1"]+0.15, y=aa$Dens_est[aa$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Auba"]-0.15, y=aa$Dens_est[aa$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="Bieb"]-0.05, y=aa$Dens_est[aa$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="O3"]+0.05, y=aa$Dens_est[aa$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=aa$yr[aa$site=="W1"]+0.15, y=aa$Dens_est[aa$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/aphelocheirus_aestivalis.svg",
  x = 2002.5,
  y = min(aa$Dens_est) * 7.5,
  width = 2.5,
  height = NULL
)
title("a. Aphelocheirus aestivalis",bty="n",cex.main=2)

###############Ancylus fluviatilis
af[is.na(af)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(af$Dens_est),max(af$Dens_est)), xlim=c(2000,2020))
#title(ylab="Density (XX)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=af$yr[af$site=="Auba"]-0.15, y=af$Dens_est[af$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="Bieb"]-0.05, y=af$Dens_est[af$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="O3"]+0.05, y=af$Dens_est[af$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="W1"]+0.15, y=af$Dens_est[af$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="Auba"]-0.15, y=af$Dens_est[af$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="Bieb"]-0.05, y=af$Dens_est[af$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="O3"]+0.05, y=af$Dens_est[af$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=af$yr[af$site=="W1"]+0.15, y=af$Dens_est[af$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/ancylus_fluviatilis.svg",
  x = 2002.5,
  y = min(af$Dens_est) * 63.5,
  width = 2.5,
  height = NULL
)
title(main="b. Ancylus fluviatilis",bty="n",cex.main=2)
 
###############Baetis rhodani  
br[is.na(br)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(br$Dens_est),max(br$Dens_est)), xlim=c(2000,2020))
#title(ylab="Density (XX)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=br$yr[br$site=="Auba"]-0.15, y=br$Dens_est[br$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Bieb"]-0.05, y=br$Dens_est[br$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="O3"]+0.05, y=br$Dens_est[br$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="W1"]+0.15, y=br$Dens_est[br$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Auba"]-0.15, y=br$Dens_est[br$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="Bieb"]-0.05, y=br$Dens_est[br$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="O3"]+0.05, y=br$Dens_est[br$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=br$yr[br$site=="W1"]+0.15, y=br$Dens_est[br$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/baetis_rhodani.svg",
  x = 2002.5,
  y = min(br$Dens_est) * 150,
  width = 4,
  height = NULL
)
title(main="c. Baetis rhodani",bty="n",cex.main=2)

#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
title(ylab=parse(text='Density/m^2'), line=-2,cex.lab=2)
         
###############Eiseniella tetraedra 
et[is.na(et)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(et$Dens_est),max(et$Dens_est)), xlim=c(2000,2020))
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=et$yr[et$site=="Auba"]-0.15, y=et$Dens_est[et$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Bieb"]-0.05, y=et$Dens_est[et$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="O3"]+0.05, y=et$Dens_est[et$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="W1"]+0.15, y=et$Dens_est[et$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Auba"]-0.15, y=et$Dens_est[et$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="Bieb"]-0.05, y=et$Dens_est[et$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="O3"]+0.05, y=et$Dens_est[et$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=et$yr[et$site=="W1"]+0.15, y=et$Dens_est[et$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/eiseniella_tetraedra.svg",
  x = 2002.5,
  y = min(et$Dens_est) *3,
  width = 3,
  height = NULL
)
title(main="d. Eiseniella tetraedra",bty="n",cex.main=2)
legend("topleft", legend=c("Aubach","Bieber","Kinzig O3","Kinzig W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2, cex=1.5)

###############Ephemera danica 
ed[is.na(ed)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ed$Dens_est),max(ed$Dens_est)), xlim=c(2000,2020))
#title(ylab="Density (XX)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ed$yr[ed$site=="Auba"]-0.15, y=ed$Dens_est[ed$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Bieb"]-0.05, y=ed$Dens_est[ed$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="O3"]+0.05, y=ed$Dens_est[ed$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="W1"]+0.15, y=ed$Dens_est[ed$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Auba"]-0.15, y=ed$Dens_est[ed$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="Bieb"]-0.05, y=ed$Dens_est[ed$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="O3"]+0.05, y=ed$Dens_est[ed$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ed$yr[ed$site=="W1"]+0.15, y=ed$Dens_est[ed$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/Ephemera_danica.svg",
  x = 2002.5,
  y = min(ed$Dens_est) + 140,
  width = 3.5,
  height = NULL
)
title(main="e. Ephemera danica ",bty="n",cex.main=2)
 
###############Gammarus roeselii
gr[is.na(gr)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(gr$Dens_est),max(gr$Dens_est)), xlim=c(2000,2020))
#title(ylab="Density (XX)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=gr$yr[gr$site=="Auba"]-0.15, y=gr$Dens_est[gr$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Bieb"]-0.05, y=gr$Dens_est[gr$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="O3"]+0.05, y=gr$Dens_est[gr$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="W1"]+0.15, y=gr$Dens_est[gr$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Auba"]-0.15, y=gr$Dens_est[gr$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="Bieb"]-0.05, y=gr$Dens_est[gr$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="O3"]+0.05, y=gr$Dens_est[gr$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=gr$yr[gr$site=="W1"]+0.15, y=gr$Dens_est[gr$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/Gammarus_roselii.svg",
  x = 2002,
  y = min(gr$Dens_est) * 30,
  width = 3.5,
  height = NULL
)
title(main="f. Gammarus roeselii",bty="n",cex.main=2)

#empty plot for where the labels will go
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

###############Hydropsyche siltalai 
hs[is.na(hs)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(hs$Dens_est),max(hs$Dens_est)), xlim=c(2000,2020))
#title(ylab="Density (XX)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=hs$yr[hs$site=="Auba"]-0.15, y=hs$Dens_est[hs$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Bieb"]-0.05, y=hs$Dens_est[hs$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="O3"]+0.05, y=hs$Dens_est[hs$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="W1"]+0.15, y=hs$Dens_est[hs$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Auba"]-0.15, y=hs$Dens_est[hs$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="Bieb"]-0.05, y=hs$Dens_est[hs$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="O3"]+0.05, y=hs$Dens_est[hs$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=hs$yr[hs$site=="W1"]+0.15, y=hs$Dens_est[hs$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/hydropsyche_siltalai.svg",
  x = 2015.5,
  y = min(hs$Dens_est) * 320,
  width = 5,
  height = NULL
)
title(main="g. Hydropsyche siltalai",bty="n",cex.main=2)
 
###############Orectochilus villosus
ov[is.na(ov)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(ov$Dens_est),max(ov$Dens_est)), xlim=c(2000,2020))
#title(ylab="Density (XX)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=ov$yr[ov$site=="Auba"]-0.15, y=ov$Dens_est[ov$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Bieb"]-0.05, y=ov$Dens_est[ov$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="O3"]+0.05, y=ov$Dens_est[ov$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="W1"]+0.15, y=ov$Dens_est[ov$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Auba"]-0.15, y=ov$Dens_est[ov$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="Bieb"]-0.05, y=ov$Dens_est[ov$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="O3"]+0.05, y=ov$Dens_est[ov$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=ov$yr[ov$site=="W1"]+0.15, y=ov$Dens_est[ov$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/Orectochilus_villosus.svg",
  x = 2017.5,
  y = min(ov$Dens_est) + 29,
  width = 2.25,
  height = NULL
)
title(main="h. Orectochilus villosus",bty="n",cex.main=2)
 
###############Prodiamesa olivacea  
po[is.na(po)] <- 0
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(po$Dens_est),max(po$Dens_est)), xlim=c(2000,2020))
#title(ylab="Density (XX)", line=2.7,cex.lab=1.5)
#title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)
points(x=po$yr[po$site=="Auba"]-0.15, y=po$Dens_est[po$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Bieb"]-0.05, y=po$Dens_est[po$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="O3"]+0.05, y=po$Dens_est[po$site=="O3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="W1"]+0.15, y=po$Dens_est[po$site=="W1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Auba"]-0.15, y=po$Dens_est[po$site=="Auba"], type="l", bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="Bieb"]-0.05, y=po$Dens_est[po$site=="Bieb"], type="l", bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="O3"]+0.05, y=po$Dens_est[po$site=="O3"], type="l", bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=po$yr[po$site=="W1"]+0.15, y=po$Dens_est[po$site=="W1"], type="l", bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
add_silhouette(
  upload_img = "Silhouettes/Prodiamesa_olivacea.svg",
  x = 2017.5,
  y = min(po$Dens_est) + 350,
  width = 4,
  height = NULL
)
title(main="i. Prodiamesa olivacea ",bty="n",cex.main=2)
title(xlab="Year", line=2.7,cex.lab=2)

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


