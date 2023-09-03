#
#install.packages("effects")
library(effects)
library(stringr)
library(tidyr)
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")
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

int <- intra[complete.cases(intra$Yryly_Temp),]
head(int)

########calculate BL estimates for a given Temp (marginal effects) for each spp and site
int$spp_site <- paste(int$SiteShort, int$SppCode)
ests_s <- NULL
for(i in unique(int$spp_site)){
  tryCatch({
  sub <- int[int$spp_site == i, ]
	sub<-sub[complete.cases(sub[, "Ldens"]),]
	sub<-sub[complete.cases(sub[, "BL"]),]
	mod<-lm(BL~ Yryly_Temp + sDOY + Ldens, data=sub)
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
colnames(els)[3] ="temp"
colnames(els)[4] ="mean"
colnames(els)[5] ="lower95"
colnames(els)[6] ="upper95"
head(els)
##############

##plot
##tiff(filename = "plots/Temp_BL.tiff", width = 7, height = 6, units = 'in', res = 600, compression = 'lzw')



####dev.off()
##

##
