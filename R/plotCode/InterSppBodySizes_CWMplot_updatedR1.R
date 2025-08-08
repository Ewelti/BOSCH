##Set working directory
setwd("C:/Users/elwel/OneDrive/Desktop/aquatic_data/git/BOSCH/")

# load libraries
library(scales)

# attach data
inter <- read.csv("output_data/CWM.csv", header=T)
head(inter)

#######################################
####################################
#no sci notation
options(scipen = 999)
options(na.action = "na.omit")

tiff(filename = "plots/CWM_overYrandTemp.tiff", width = 10, height = 6, units = 'in', res = 600, compression = 'lzw')
par(mar=c(4,4,0.2,0.4),mfrow=c(1,2))

plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(inter$CWM),max(inter$CWM)), xlim=c(2000.7,2019.3))
legend(legend="a","topright", bty="n",cex=2)
title(ylab="CWM of body size (mg)", line=2.5,cex.lab=1.5)
title(xlab="Sampling year", line=2.5,cex.lab=1.5)
box(lwd=3)

mo_pmod <- lmer(CWM ~ year + (1|site), data = inter)
mo_e <- Effect("year", partial.residuals=T, mo_pmod)
mo_e <- data.frame(mo_e)
points(mo_e$fit ~ mo_e$year,type="l",col="gray20",lwd=4)
polygon(c(rev(mo_e$year), mo_e$year), c(rev(mo_e$upper), mo_e$lower), col=alpha(1,0.2), border = NA)

points(x=inter$year[inter$site=="Auba"], y=inter$CWM[inter$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=inter$year[inter$site=="Bieb"], y=inter$CWM[inter$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=inter$year[inter$site=="KiO3"], y=inter$CWM[inter$site=="KiO3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=inter$year[inter$site=="KiW1"], y=inter$CWM[inter$site=="KiW1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
#legend("topleft", legend=c("Aubach","Bieber","KiO3","KiW1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2.5, cex=1.5)
abline(lm(inter$CWM[inter$site=="Auba"] ~ inter$year[inter$site=="Auba"]),lwd=2,col=alpha(1,0.6),lty=1)
abline(lm(inter$CWM[inter$site=="Bieb"] ~ inter$year[inter$site=="Bieb"]),lwd=2,col=alpha(2,0.6),lty=1)
abline(lm(inter$CWM[inter$site=="KiO3"] ~ inter$year[inter$site=="KiO3"]),lwd=2,col=alpha(3,0.6),lty=1)
abline(lm(inter$CWM[inter$site=="KiW1"] ~ inter$year[inter$site=="KiW1"]),lwd=2,col=alpha(4,0.6),lty=1)

##########temp plot
plot(1, 1, type= "n",las=1,main="",cex.main=1.5,ylab="", xlab="", ylim=c(min(inter$CWM),max(inter$CWM)), xlim=c(min(inter$Yryly_Temp),max(inter$Yryly_Temp)))
legend(legend="b","topright", bty="n",cex=2)
#title(ylab="CWM of body size (mm)", line=2.5,cex.lab=1.5)
title(xlab="Mean temperature (\u00B0C) of prior 12 mo", line=2.5,cex.lab=1.5)
box(lwd=3)

mo_pmod <- lmer(CWM ~ Yryly_Temp + (1|site), data = inter)
mo_e <- Effect("Yryly_Temp", partial.residuals=T, mo_pmod)
mo_e <- data.frame(mo_e)
points(mo_e$fit ~ mo_e$Yryly_Temp,type="l",col="gray20",lwd=4)
polygon(c(rev(mo_e$Yryly_Temp), mo_e$Yryly_Temp), c(rev(mo_e$upper), mo_e$lower), col=alpha(1,0.2), border = NA)

points(x=inter$Yryly_Temp[inter$site=="Auba"], y=inter$CWM[inter$site=="Auba"], pch=21, bg=alpha(1,0.6),col=alpha(1,0.6),lwd=2,cex=2.5)
points(x=inter$Yryly_Temp[inter$site=="Bieb"], y=inter$CWM[inter$site=="Bieb"], pch=22, bg=alpha(2,0.6),col=alpha(2,0.6),lwd=2,cex=2.5)
points(x=inter$Yryly_Temp[inter$site=="KiO3"], y=inter$CWM[inter$site=="KiO3"], pch=23, bg=alpha(3,0.6),col=alpha(3,0.6),lwd=2,cex=2.5)
points(x=inter$Yryly_Temp[inter$site=="KiW1"], y=inter$CWM[inter$site=="KiW1"], pch=24, bg=alpha(4,0.6),col=alpha(4,0.6),lwd=2,cex=2.5)
legend("topleft", legend=c("Aubach","Bieber","Kinzig O3","Kinzig W1"),col=c(1,2,3,4),pt.bg=c(1,2,3,4),pt.lwd=1, pch=c(21,22,23,24),lty=0,lwd=2,bty="n",pt.cex=2.5, cex=1.5)

abline(lm(inter$CWM[inter$site=="Auba"] ~ inter$Yryly_Temp[inter$site=="Auba"]),lwd=2,col=alpha(1,0.6),lty=2)
abline(lm(inter$CWM[inter$site=="Bieb"] ~ inter$Yryly_Temp[inter$site=="Bieb"]),lwd=2,col=alpha(2,0.6),lty=1)
abline(lm(inter$CWM[inter$site=="KiO3"] ~ inter$Yryly_Temp[inter$site=="KiO3"]),lwd=2,col=alpha(3,0.6),lty=1)
abline(lm(inter$CWM[inter$site=="KiW1"] ~ inter$Yryly_Temp[inter$site=="KiW1"]),lwd=2,col=alpha(4,0.6),lty=2)

dev.off()

##
##