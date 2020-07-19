install.packages("readstata13")
install.packages("glmmfields")
install.packages("ape")
install.packages("spatialprobit")
install.packages("spatialreg")
install.packages("MASS")
install.packages("tidyverse")
install.packages("rgdal")
install.packages("spdep")
install.packages("spgwr")
install.packages("spsur")
install.packages("sf")
install.packages("dplyr")
library(spatialreg)
library(spatialprobit)
library(ape)
library(readstata13)
library(glmmfields)
library(MASS)
library(tidyverse)
library(rgdal)
library(spdep)
library(spgwr)
library(spsur)
library(sf)
library(dplyr)
spatialdataSAR<- read.dta13("/Users/Saiem 70/Dropbox/Research stuff/Alok bohara research/Dissertation/Chapter 1/All the data work/My data/With long form/SAR results/SAR Analysis data.dta")

############################
#Creating the weight matrix

require(Matrix)
sp_point_SAR <- cbind(spatialdataSAR$Longitude, spatialdataSAR$Latitude)
colnames(sp_point_SAR) <- c("Longitude","Latitude")
head(sp_point_SAR)
I_n <- sparseMatrix(i = 1:n, j = 1:n, x = 1)
nbSAR <- knn2nb(knearneigh(x=sp_point_SAR, k=8))
listwSAR <- nb2listw(nbSAR, style = "W",zero.policy = TRUE)


###############################################
#Simulatd Morans I and corelogram for recycling

M.recycle<-moran.test(spatialdataSAR$meanrecycle, listw=listwSAR, alternative="two.sided")
R<-M.recycle$estimate[1]

set.seed(123456)
MWTP_Recycle<-moran.mc(spatialdata$E,listw=listw,nsim=999)
hist(MWTP_Recycle$res, freq=TRUE, breaks=20, xlab="Simulated Moran's I")
abline(v=R,col="red")

cor.recycle<-sp.correlogram(nbSAR,spatialdataSAR$meanrecycle, order=8, method="I", style="W",zero.policy = TRUE)
print(cor.recycle); 
plot(cor.recycle, main="Spatial Correlogram", type="l", col="blue")

#Simulatd Morans I and corelogram for composting
moran.test(spatialdataSAR$meancomposting, listw=listwSAR, alternative="two.sided")
M.compost<-moran.test(spatialdataSAR$meancomposting, listw=listwSAR, alternative="two.sided")
C<-M.compost$estimate[1]

set.seed(123456)
MWTP_Composting<-moran.mc(spatialdataSAR$meancomposting,listw=listwSAR,nsim=999)
hist(MWTP_Composting$res, freq=TRUE, breaks=20, xlab="Simulated Moran's I")
abline(v=C,col="red")


cor.compost<-sp.correlogram(nbSAR,spatialdataSAR$meancomposting, order=8, method="I", style="W",zero.policy = TRUE)
print(cor.compost); 
plot(cor.compost, main="Spatial Correlogram", type="l", col="blue")

#Simulatd Morans I and corelogram for dumping
moran.test(spatialdataSAR$meandumping, listw=listwSAR, alternative="two.sided")
M.dumping<-moran.test(spatialdataSAR$meandumping, listw=listwSAR, alternative="two.sided")
D<-M.dumping$estimate[1]

set.seed(123456)
MWTP_Dumping<-moran.mc(spatialdataSAR$meandumping,listw=listwSAR,nsim=999)
hist(MWTP_Dumping$res, freq=TRUE, breaks=20, xlab="Simulated Moran's I")
abline(v=D,col="red")


cor.dumping<-sp.correlogram(nbSAR,spatialdataSAR$meandumping, order=8, method="I", style="W",zero.policy = TRUE)
print(cor.dumping); 
plot(cor.dumping, main="Spatial Correlogram", type="l", col="blue")

#Simulatd Morans I and corelogram for wastecollected55
moran.test(spatialdataSAR$meanwastecollected2, listw=listwSAR, alternative="two.sided")
M.wastecollected55<-moran.test(spatialdataSAR$meanwastecollected2, listw=listwSAR, alternative="two.sided")
W1<-M.wastecollected55$estimate[1]

set.seed(123456)
MWTP_wastecollected55<-moran.mc(spatialdataSAR$meanwastecollected2,listw=listwSAR,nsim=999)
hist(MWTP_wastecollected55$res, freq=TRUE, breaks=20, xlab="Simulated Moran's I")
abline(v=W1,col="red")


cor.wastecollected55<-sp.correlogram(nbSAR,spatialdataSAR$meanwastecollected2, order=8, method="I", style="W",zero.policy = TRUE)
print(cor.wastecollected55); 
plot(cor.wastecollected55, main="Spatial Correlogram", type="l", col="blue")

#Simulatd Morans I and corelogram for wastecollected80
moran.test(spatialdataSAR$meanwastecollected3, listw=listwSAR, alternative="two.sided")
M.wastecollected80<-moran.test(spatialdataSAR$meanwastecollected3, listw=listwSAR, alternative="two.sided")
W2<-M.wastecollected80$estimate[1]

set.seed(123456)
MWTP_wastecollected80<-moran.mc(spatialdataSAR$meanwastecollected3,listw=listwSAR,nsim=999)
hist(MWTP_wastecollected80$res, freq=TRUE, breaks=20, xlab="Simulated Moran's I")
abline(v=W2,col="red")


cor.wastecollected80<-sp.correlogram(nbSAR,spatialdataSAR$meanwastecollected3, order=8, method="I", style="W",zero.policy = TRUE)
print(cor.wastecollected80); 
plot(cor.wastecollected80, main="Spatial Correlogram", type="l", col="blue")


##################################################
#Separate SAR model for each of the attributes

K <- cbind(spatialdataSAR$male, spatialdataSAR$hhmember,  spatialdataSAR$Age, spatialdataSAR$High_Income, spatialdataSAR$High_Education)

sar1<- lagsarlm(spatialdataSAR$meanrecycle ~ Heard_Recycling_CampaignY+K, data = spatialdataSAR, listw = listwSAR, zero.policy = T)
summary(sar1)
AIC(sar1)

impacts.sarlm(sar1, listw=listwSAR)

sar2<- lagsarlm(spatialdataSAR$meancomposting ~ Heard_Recycling_CampaignY+Heard_Composting_CampaignY+Howsupportive_Landfill_SiteHI+Howsupportive_EfficiencyHI+K, data = spatialdataSAR, listw = listwSAR, zero.policy = T)
summary(sar2)
AIC(sar2)

sar3<- lagsarlm(spatialdataSAR$meandumping ~ Heard_Recycling_CampaignY+Heard_Composting_CampaignY+Howsupportive_Landfill_SiteHI+Howsupportive_EfficiencyHI+K, data = spatialdataSAR, listw = listwSAR, zero.policy = T)
summary(sar3)
AIC(sar3)


sar4<- lagsarlm(spatialdataSAR$meanwastecollected2 ~ Heard_Recycling_CampaignY+Heard_Composting_CampaignY+Howsupportive_Landfill_SiteHI+Howsupportive_EfficiencyHI+K, data = spatialdataSAR, listw = listwSAR, zero.policy = T)
summary(sar4)
AIC(sar4)

sar5<- lagsarlm(spatialdataSAR$meanwastecollected3 ~ Heard_Recycling_CampaignY+Heard_Composting_CampaignY+Howsupportive_Landfill_SiteHI+Howsupportive_EfficiencyHI+K, data = spatialdataSAR, listw = listwSAR, zero.policy = T)
summary(sar5)
AIC(sar5)

#################################
#SUR-SAR model for the attributes
formula_sur <- spatialdataSAR$meanrecycle |spatialdataSAR$meancomposting| spatialdataSAR$meandumping|spatialdataSAR$meanwastecollected2|spatialdataSAR$meanwastecollected3 ~ spatialdataSAR$Heard_Recycling_CampaignY+K|
               spatialdataSAR$Heard_Composting_CampaignY+ K|spatialdataSAR$Howsupportive_Landfill_SiteHI+K|spatialdataSAR$Howsupportive_EfficiencyHI+K|Howsupportive_EfficiencyHI+K

control <- list(fdHess = TRUE)
sur.slm1 <- spsurml(formula = formula_sur, listw = listwSAR, 
                    method = "LU", type = "slm", 
                    control = control,  
                   data = spatialdataSAR)
summary(sur.slm1)


sur.slm2 <- spsur3sls(formula = formula_sur, listw = listwSAR, 
                        maxlagW = 3, type = "slm", 
                        data = spatialdataSAR)
summary(sur.slm2)
