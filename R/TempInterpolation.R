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
plot(kio3_interp_temp,type="l", ylab="Celsius", xlab="Hours") #check

##gaps that are still problems
#1801-10706
#27937-29328
#72313-76071
length(n_kio$kio3_interp_temp) #92760
kio_gap <- c(n_kio$kio3_interp_temp[1:1800],rep(NA,(10706-1801+1)),
n_kio$kio3_interp_temp[10707:27936],rep(NA,(29328-27937+1)),
n_kio$kio3_interp_temp[29329:72312],rep(NA,(76071-72313+1)),
n_kio$kio3_interp_temp[76072:92760])
length(kio_gap)
plot(kio_gap~ c(1:92760), type="l", ylab="Celsius", xlab="Hours") #check

n_kio<-cbind(kio,kio_gap)
n_kio <- cbind(n_kio[,1:5],n_kio[,7:8])
head(n_kio)

##KiW1
kiw <- read.csv("RawData/Temp_data/KiW1_hourTemp.csv", header=T)
kiw<-make_JDay(kiw)
head(kiw)

#Interpolate gaps
interp<-interpolate_gaps_hourly(hourtemps=kiw,latitude=50.2)
kiw1_interp_temp <- as.numeric(interp$weather$Temp)
n_kiw<-cbind(kiw,kiw1_interp_temp)
head(n_kiw)
plot(n_kiw$kiw1_interp_temp,type="l", ylab="Celsius", xlab="Hours") #check

##gaps that are still problems
#27937-30168
#72305-76069
length(n_kiw$kiw1_interp_temp) #92760
kiw_gap <- c(n_kiw$kiw1_interp_temp[1:27936],rep(NA,(30168-27937+1)),
n_kiw$kiw1_interp_temp[30169:72304],rep(NA,(76069-72305+1)),
n_kiw$kiw1_interp_temp[76070:92760])
length(kiw_gap)
plot(kiw_gap~ c(1:92760), type="l", ylab="Celsius", xlab="Hours") #check

n_kiw<-cbind(kiw,kiw_gap)
head(n_kiw)

#merge
sub_b <- cbind("code"=n_kiw[,1],"kiw_gap"=n_kiw[,8])
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

##gaps that are still problems
#17292-19643
length(n_au$aubach_interp_temp) #80832
au_gap <- c(n_au$aubach_interp_temp[1:17291],rep(NA,(19643-17292+1)),
n_au$aubach_interp_temp[19644:80832])
length(au_gap)
plot(au_gap~ c(1:80832), type="l", ylab="Celsius", xlab="Hours") #check

n_au<-cbind(au,au_gap)
head(n_au)

#merge
sub_b <- cbind("code"=n_au[,1],"au_gap"=n_au[,8])
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

##gaps that are still problems
#15913-28346
length(n_bieb$bieb_interp_temp) #89520
bieb_gap <- c(n_bieb$bieb_interp_temp[1:15912],rep(NA,(28346-15913+1)),
n_bieb$bieb_interp_temp[28347:89520])
length(bieb_gap)
plot(bieb_gap~ c(1:89520), type="l", ylab="Celsius", xlab="Hours") #check

n_bieb<-cbind(bieb,bieb_gap)
head(n_bieb)

#merge
sub_b <- cbind("code"=n_bieb[,1],"bieb_gap"=n_bieb[,8])
head(sub_b)
m_interp <- merge(m_interp, sub_b, by="code", all=T)
head(m_interp)
m_interp <- m_interp[order(m_interp$Year, m_interp$JDay, m_interp$Hour),]

#write output
write.csv(m_interp, "output_data/InterpolatedTemp.csv")
###############################################################
###############################################
##################################
######################

##attach data
int <- read.csv("output_data/InterpolatedTemp.csv", header=T)

plot(int$kio_gap~c(1:dim(int)[1]), type='l')



