##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

##load library
library(data.table)
library(nlme)
library(lme4)
library(scales)

# attach data
spp <- read.csv("RawData/RMO4sites_measuredTaxa_updatedR1.csv", header=T)
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

##test model- first spp
Af <- spp[which(spp$taxa=="Ancylus fluviatilis"), ]
trend.i <- summary(gls(Lab ~ poly(sDOY,2) + syr + site,na.action=na.omit, data = Af))#$tTable[4, c(1,2,4)]
trend.i

##################################
###################################
#calculate trends
trends <- NULL
for(i in unique(spp$taxa)){
  tryCatch({
    sub <- spp[spp$taxa == i, ]
    trend.i <- summary(gls(Lab ~ poly(sDOY,2) + syr + site,na.action=na.omit, data = sub))$tTable[4, c(1,2,4)]
    trend.i <- data.frame(taxa = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)

trends

########################################
################################################
#save data
write.csv(trends,"output_data/Measured9SppTrends_overall.csv")
#############################################
######################################

##test model- first spp
Af_A <- spp[which(spp$spp_site=="Ancylus fluviatilis Auba"), ]
trend.i <- summary(gls(Lab ~ poly(sDOY,2) + syr,na.action=na.omit, data = Af))$tTable[4, c(1,2,4)]
trend.i

##################################
###################################
#calculate trends- site specific
spp$spp_site <- paste(spp$taxa, spp$site)
head(spp)
trends <- NULL
for(i in unique(spp$spp_site)){
  tryCatch({
  sub <- spp[spp$spp_site == i, ]
	#sub<-sub[complete.cases(sub[, "Lab"]),]
    trend.i <- summary(gls(Lab ~ poly(sDOY,2) + syr,na.action=na.omit, data = sub))$tTable[4, c(1,2,4)]
    trend.i <- data.frame(spp_site = i, 
                        t(trend.i))
    trends <- rbind(trends, trend.i) ; rm(trend.i, sub)
    }, error=function(e){cat(unique(sub$taxa),conditionMessage(e), "\n")})
} ; rm(i)

trends

########################################
################################################
#save data
write.csv(trends,"output_data/Measured9SppTrends_bySite.csv")
#############################################
######################################
