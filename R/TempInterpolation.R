##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

##load libraries
library("chillR")

# attach data
##KiO3
kio <- read.csv("RawData/Temp_data/KiO3_hourTemp.csv", header=T)
kio<-make_JDay(kio)
head(kio)

#Interpolate gaps
interp<-interpolate_gaps_hourly(hourtemps=kio,latitude=50.2)
kio3_interp_temp <- as.numeric(interp$weather$Temp)
n_kio<-cbind(kio,kio3_interp_temp)
head(n_kio)
plot(n_kio$kio3_interp_temp,type="l", ylab="Celsius", xlab="Hours") #check

##KiW1
kiw <- read.csv("RawData/Temp_data/KiW1_hourTemp.csv", header=T)
kiw<-make_JDay(kiw)
head(kiw)

#Interpolate gaps
interp<-interpolate_gaps_hourly(hourtemps=au,latitude=50.2)
kiw1_interp_temp <- as.numeric(interp$weather$Temp)
n_kiw<-cbind(kiw,kiw1_interp_temp)
head(n_kiw)
plot(n_kiw$kiw1_interp_temp,type="l", ylab="Celsius", xlab="Hours") #check

#merge
sub_b <- cbind("code"=n_kiw[,1],"kiw_interp_temp"=n_kiw[,8])
head(sub_b)
m_interp <- merge(n_kio, sub_b, by="code", all=T)
head(m_interp)

##Aubach
au <- read.csv("RawData/Temp_data/aubach_hourTemp.csv", header=T)
au<-make_JDay(au)
head(au)

#Interpolate gaps
interp<-interpolate_gaps_hourly(hourtemps=au,latitude=50.2)
aubach_interp_temp <- as.numeric(interp$weather$Temp)
n_au<-cbind(au,aubach_interp_temp)
head(n_au)
plot(n_au$aubach_interp_temp,type="l", ylab="Celsius", xlab="Hours") #check

#merge
sub_b <- cbind("code"=n_au[,1],"au_interp_temp"=n_au[,8])
head(sub_b)
m_interp <- merge(m_interp, sub_b, by="code", all=T)
head(m_interp)

##Bieber
bieb <- read.csv("RawData/Temp_data/bieb_hourTemp.csv", header=T)
bieb<-make_JDay(bieb)
head(bieb)

#Interpolate gaps
interp<-interpolate_gaps_hourly(hourtemps=bieb,latitude=50.2)
bieb_interp_temp <- as.numeric(interp$weather$Temp)
n_bieb<-cbind(bieb,bieb_interp_temp)
head(n_bieb)
plot(n_bieb$bieb_interp_temp,type="l", lwd=0.5,ylab="Celsius", xlab="Hours") #check

#merge
sub_b <- cbind("code"=n_bieb[,1],"bieb_interp_temp"=n_bieb[,8])
head(sub_b)
m_interp <- merge(m_interp, sub_b, by="code", all=T)
head(m_interp)

#write output
write.csv(m_interp, "output_data/InterpolatedTemp.csv")



