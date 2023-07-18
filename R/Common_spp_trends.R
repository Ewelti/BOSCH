##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

##load library
library(data.table)
library(nlme)
library(lme4)
library(scales)

# attach data
spp <- read.csv("RawData/RMO4sites.csv", header=T)
head(spp)

##scale everything
#function to add a new column onto the data with scaled vars (with s before their name)
scaleVars <- function(df){
  newd <- plyr::numcolwise(scale)(df)
  names(newd) <- sapply(names(newd),function(x)paste0("s",x))
  cbind(df, newd)
}
#apply function
spp <- scaleVars(spp)
head(spp)

spp$Lab <- log10(spp$Abundance_per_m2)

#subset by site
Au <- spp[which(spp$site=="Auba"), ]
Bi <- spp[which(spp$site=="Bieb"), ]
KiO3 <- spp[which(spp$site=="KiO3"), ]
KiW1 <- spp[which(spp$site=="KiW1"), ]

############Aubach taxa with unique yrs############################################
com <- NULL
for(i in unique(Au$ID_Art)){
    sub <- Au[Au$ID_Art == i, ]
    numY <- length(unique(sub$yr))
    dd <- data.frame(ID_Art = i, numY)
    com <- rbind(com, dd) ; rm(sub, numY, dd)
} ; rm(i)

unique(KiW1$yr)
##list spp present in at least 12 yrs
au_LT <- com[ which(com$numY > 11),]
nrow(au_LT)

##subset site data to only these more common spp
Au_com <-subset(Au, ID_Art %in% au_LT$ID_Art)
head(Au_com)

#calculate trends
trends <- NULL
for(i in unique(Au_com$ID_Art)){
  tryCatch({
    sub <- Au_com[Au_com$ID_Art == i, ]
    trend.i <- summary(gls(Lab ~ poly(sDOY,2) + syr,na.action=na.omit, data = sub))$tTable[4, c(1,2,4)]
    trend.i <- data.frame(ID_Art = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$ID_Art),conditionMessage(e), "\n")})
} ; rm(i)

site<-rep("Aubach",times=nrow(trends))
au_t<- cbind(site,trends)
head(au_t)

############Bieber taxa with unique yrs############################################
com <- NULL
for(i in unique(Bi$ID_Art)){
    sub <- Bi[Bi$ID_Art == i, ]
    numY <- length(unique(sub$yr))
    dd <- data.frame(ID_Art = i, numY)
    com <- rbind(com, dd) ; rm(sub, numY, dd)
} ; rm(i)

##list spp present in at least 12 yrs
Bi_LT <- com[ which(com$numY > 11),]
nrow(Bi_LT)

##subset site data to only these more common spp
Bi_com <-subset(Bi, ID_Art %in% Bi_LT$ID_Art)
head(Bi_com)

#calculate trends
trends <- NULL
for(i in unique(Bi_com$ID_Art)){
  tryCatch({
    sub <- Bi_com[Bi_com$ID_Art == i, ]
    trend.i <- summary(gls(Lab ~ poly(sDOY,2) + syr,na.action=na.omit, data = sub))$tTable[4, c(1,2,4)]
    trend.i <- data.frame(ID_Art = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$ID_Art),conditionMessage(e), "\n")})
} ; rm(i)

site<-rep("Bieber",times=nrow(trends))
Bi_t<- cbind(site,trends)
head(Bi_t)
TTRR <-rbind(au_t,Bi_t)

############KiO3 taxa with unique yrs############################################
com <- NULL
for(i in unique(KiO3$ID_Art)){
    sub <- KiO3[KiO3$ID_Art == i, ]
    numY <- length(unique(sub$yr))
    dd <- data.frame(ID_Art = i, numY)
    com <- rbind(com, dd) ; rm(sub, numY, dd)
} ; rm(i)

##list spp present in at least 12 yrs
KiO3_LT <- com[ which(com$numY > 11),]
nrow(KiO3_LT)

##subset site data to only these more common spp
KiO3_com <-subset(KiO3, ID_Art %in% KiO3_LT$ID_Art)
head(KiO3_com)

#calculate trends
trends <- NULL
for(i in unique(KiO3_com$ID_Art)){
  tryCatch({
    sub <- KiO3_com[KiO3_com$ID_Art == i, ]
    trend.i <- summary(gls(Lab ~ poly(sDOY,2) + syr,na.action=na.omit, data = sub))$tTable[4, c(1,2,4)]
    trend.i <- data.frame(ID_Art = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$ID_Art),conditionMessage(e), "\n")})
} ; rm(i)

site<-rep("KiO3",times=nrow(trends))
KiO3_t<- cbind(site,trends)
head(KiO3_t)
TTRR <-rbind(TTRR,KiO3_t)

############KiW1 taxa with unique yrs############################################
com <- NULL
for(i in unique(KiW1$ID_Art)){
    sub <- KiW1[KiW1$ID_Art == i, ]
    numY <- length(unique(sub$yr))
    dd <- data.frame(ID_Art = i, numY)
    com <- rbind(com, dd) ; rm(sub, numY, dd)
} ; rm(i)

##list spp present in at least 12 yrs
KiW1_LT <- com[ which(com$numY > 11),]
nrow(KiW1_LT)

##subset site data to only these more common spp
KiW1_com <-subset(KiW1, ID_Art %in% KiW1_LT$ID_Art)
head(KiW1_com)

#calculate trends
trends <- NULL
for(i in unique(KiW1_com$ID_Art)){
  tryCatch({
    sub <- KiW1_com[KiW1_com$ID_Art == i, ]
    trend.i <- summary(gls(Lab ~ poly(sDOY,2) + syr,na.action=na.omit, data = sub))$tTable[4, c(1,2,4)]
    trend.i <- data.frame(ID_Art = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$ID_Art),conditionMessage(e), "\n")})
} ; rm(i)

site<-rep("KiW1",times=nrow(trends))
KiW1_t<- cbind(site,trends)
head(KiW1_t)
TTRR <-rbind(TTRR,KiW1_t)
nrow(TTRR)
length(unique(TTRR$ID_Art))

########################################
################################################
#get sizes
sizes <- unique(spp[,c("ID_Art","spp_mm")])
si <- merge(TTRR,sizes,by="ID_Art")
#############################################
######################################

########################################
################################################
#save data
write.csv(si,"output_data/CommonSppTrends.csv")
#############################################
######################################

##all sites
sizemsite <- lm(si$Value ~ si$spp_mm + si$site)
summary(sizemsite)

sizem <- lmer(si$Value ~ si$spp_mm + (1|si$site))
##sizem <- lmer(si$Value ~ poly(si$spp_mm,2) + (1|si$site))
summary(sizem)
coefs <- data.frame(coef(summary(sizem)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

tiff(filename = "plots/CommonSppTrends_overSize.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
par(mar=c(4,4,0.2,0.2))

plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(-0.5,0.5), xlim=c(0,6))
title(ylab="Slope of species abundances", line=2.5,cex.lab=1.5)
title(xlab="Species size (mm)", line=2.5,cex.lab=1.5)
polygon(x=c(-4,-4,10,10),
        y=c(-100,0,0,-100), col = "grey80", border = "grey80")
box(lwd=3)
points(x=si$spp_mm[si$site=="Aubach"], y=si$Value[si$site=="Aubach"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=si$spp_mm[si$site=="Bieber"], y=si$Value[si$site=="Bieber"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=si$spp_mm[si$site=="KiO3"], y=si$Value[si$site=="KiO3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=si$spp_mm[si$site=="KiW1"], y=si$Value[si$site=="KiW1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
legend("bottomright", legend=c("Auba","Bieb","O3","W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2.5, cex=1.5)
abline(lm(si$Value ~ si$spp_mm), lwd=2)
abline(lm(si$Value[si$site=="Aubach"] ~ si$spp_mm[si$site=="Aubach"]),lwd=2,col=alpha(1,0.6),lty=2)
abline(lm(si$Value[si$site=="Bieber"] ~ si$spp_mm[si$site=="Bieber"]),lwd=2,col=alpha(2,0.6),lty=2)
abline(lm(si$Value[si$site=="KiO3"] ~ si$spp_mm[si$site=="KiO3"]),lwd=2,col=alpha(3,0.6),lty=2)
abline(lm(si$Value[si$site=="KiW1"] ~ si$spp_mm[si$site=="KiW1"]),lwd=2,col=alpha(4,0.6),lty=2)

dev.off()

#########################################

##site specific trends- not that interesting

#subset by site
Ausi <- si[which(si$site=="Aubach"), ]
Bisi <- si[which(si$site=="Bieber"), ]
KiO3si <- si[which(si$site=="KiO3"), ]
KiW1si <- si[which(si$site=="KiW1"), ]

##Aubach
sizem <- lm(Ausi$Value ~ Ausi$spp_mm)
summary(sizem)
plot(Ausi$Value ~ Ausi$spp_mm)

##Bieber
sizem <- lm(Bisi$Value ~ Bisi$spp_mm)
summary(sizem)
plot(Bisi$Value ~ Bisi$spp_mm)

##KiO3
sizem <- lm(KiO3si$Value ~ KiO3si$spp_mm)
summary(sizem)
plot(KiO3si$Value ~ KiO3si$spp_mm)

##KiW1
sizem <- lm(KiW1si$Value ~ KiW1si$spp_mm)
summary(sizem)
plot(KiW1si$Value ~ KiW1si$spp_mm)