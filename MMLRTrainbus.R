# Title:                  Melbourne Train-bus MMLR
# Author details:         Laura Aston
# Affiliation:            Public Transport Research Group, Monash University
# Contact details:        laura.aston@monash.edu
# Script and data info:   This script performs pretreatment (factor analysis) followed by multivariate multiple linear regression of built environment and sociodemographic variables on transit ridership for a sample of train and bus locations in Melbourne.  
# Data:                   Ridership data includes average normal (school) weekday ridership 
#                         Train ridership, averaged for 2018, by Victorian Department of Transport. 
#                         Bus ridership, averaged for 4 months from August - November 2018, provided by Chris Loader from the bus planning team at the Victorian Department of Transport
#                         Built environment data is agrgegated to the transit walk catchment station area, and harvested from a variety of open data sources. For more information, refer to the publication related to this script. 
#                         Copyright statement: This script is the product of Laura Aston

# Multivariate Multiple Regression steps
# 1. Factor analysis if high collinearity between variables
# 2: Run correlation analysis (using VIF inspection) among the explanatory factors, and factors with high correlation are excluded
# 3.	Run simple regression (using one factor at a time) and keep only the factors that are significant at certain level for further analysis – unadjusted model
# 4.	Run multiple regression analysis using all the factors from step 2 – maximally adjusted model
# 5.	Remove statistically insignificant factors step by step from step 3 – parsimonious model
# 6.  Run diagnostic tests to check conformance of the solution with assumptions

#A useful resource for Multivariate Multiple Regression Analysis in R is: https://data.library.virginia.edu/getting-started-with-multivariate-multiple-regression/


#Senstitivty testing
#Sensitivity testing is needed to determine the impact of methodological choices on results (Aston et al 2019). The following variants will be tested for each step:
#Trainbus sample (different set of options explored for trainbus sample)
#a) Radius: 800m vs. 600m (factorisation, through to full solution)

#Set working directory
setwd("C:/Users/lkast1/Google Drive/PhD/2.Analysis/2. Empirical Analysis/BE-TR_Multi Country Samples/Melbourne/MMLR/BE-TU-Melbourne-MMLR")


#Install and load packages
install.packages("lm.beta")
install.packages("dplyr")
#install.packages("car")
install.packages("Hmisc")
install.packages("psych")
install.packages("car")
library(lm.beta)
library(dplyr)
#library(car)
library(Hmisc)
library(psych)
library(car) #needed for VIF

#turn off scientific notation
options(scipen = 999)

#read in data
MMLR_Data<-read.csv(file="Updated_MMLR_Data.csv")

#Transform patronage into natural log and add as a column to the dataframe
MMLR_Data<-mutate(MMLR_Data,
                  ln_Bus = log(Updated_Patronage_Bus),
                  ln_Tram = log(Patronage_Tram),
                  ln_Train = log(Patronage_Train),)

#optional - assign the Sample ID as the row names
row.names(MMLR_Data) <- MMLR_Data[,c(3)]

#First subset to test for trambus sample is all stops, with walk catchment radius of 600 metres
#subsetting for mode in the census of all eligible stops


Melb.Trainbus.600<- MMLR_Data[ which(MMLR_Data$Type=='Trainbus'
                                    & MMLR_Data$Radius=='600'
                                    & MMLR_Data$Set_Sample_ID == 'C'),]

Melb.Trainbus.800<- MMLR_Data[ which(MMLR_Data$Type=='Trainbus'
                                     & MMLR_Data$Radius=='800'
                                     & MMLR_Data$Set_Sample_ID == 'C'),]

#Assumption 1:  'Set ID' refers to the method of defining the 'walk catchment'. 'D' (distributed) means     catchments were separately calculated for the tram and bus points, then merged. 'C' (centroid) means the geographic centroid of the points was first found, and the catchment estimated from there. Data corresponding to the 'centroid' (C) catchments is used in this analysis. 

#Variables to include in the factor analysis (columns 19, 21,22, 24 - 43)
#X2_EmpDen+X4_PopDen+X5_Dwelling+X7_PropComm+X8_RetailEmp+X9_Balance+X10_Entropy+X11_HousingDiv+X12_Intersections+X13_PBN+X14_DestScore+X15_DestCount+X16_DistCBD+X17_ACCount+X18_ACNear+X19_FTZ+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X25_Rural+X26_Access

#Step 1.1: form a data frame that comprises only the variables to be included in the factor analysis (built environment)
#include variables measured on a count or continuous scale. This means "FTZ" and "AC_Count" should be excluded. Also "rural" is 0 for almost all responses, so exclude
fa.data.Melb.Trainbus.800<-Melb.Trainbus.800[,c(19, 21, 22, 24:33,35,41,43)]
fa.data.Melb.Trainbus.600<-Melb.Trainbus.600[,c(19, 21, 22, 24:33,35,41,43)]

#Step 1.2: Specify number of factors. Based on theory, will try four or five factors, consituting 1) density 2) diversity 3) design 4) regional accessibility and 5) local accessibility/walkability (Ewing and Cevero 2010, Voulgaris et al. 2017)

#alternatively, check scree plot
install.packages("nFactors")
library(nFactors)
ev_trainbus_600 <- eigen(cor(fa.data.Melb.Trainbus.600))
ev_trainbus_800 <- eigen(cor(fa.data.Melb.Trainbus.800))# get eigenvalues
ap_trainbus_600 <- parallel(subject=nrow(fa.data.Melb.Trainbus.600),var=ncol(fa.data.Melb.Trainbus.600),
               rep=100,cent=.05)
nS_trainbus_600 <- nScree(x=ev_trainbus_600$values, aparallel=ap_trainbus_600$eigen$qevpea)

ap_trainbus_800 <- parallel(subject=nrow(fa.data.Melb.Trainbus.800),var=ncol(fa.data.Melb.Trainbus.800),
                            rep=100,cent=.05)
nS_trainbus_800 <- nScree(x=ev_trainbus_800$values, aparallel=ap_trainbus_800$eigen$qevpea)

plotnScree(nS_trainbus_800) #3 eigenvalues
plotnScree(nS_trainbus_600) #3/4 eigenvalues


#Step 1.3 run factor analysis
#800
fa.Melb.Trainbus.800.3<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "none")
fa.Melb.Trainbus.800.3
#high uniqueness: proportion commerical, Entropy, DestScore, ACNear ,Urban

#600
fa.Melb.Trainbus.600.3<-factanal(fa.data.Melb.Trainbus.600, factors = 3, rotation = "none")
fa.Melb.Trainbus.600.3
#high uniqueness: proportion commerical, entropy, PBN, DistCBD, ACNear, Urban

fa.Melb.Trainbus.600.4<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "none")
fa.Melb.Trainbus.600.4
#high uniqueness: proportion commerical, entropy, ACNear

#Step 1.4 Remove variables with high uniqueness values: proportion commerical, ACNear, urban
fa.data.Melb.Trainbus.800<-Melb.Trainbus.800[,c(19, 21, 22, 25,26, 28:30, 32:33,43)]
fa.data.Melb.Trainbus.600<-Melb.Trainbus.600[,c(19, 21, 22, 25, 26, 28:33,41,43)]


fa.Melb.Trainbus.800.3<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "none")
fa.Melb.Trainbus.800.3

fa.Melb.Trainbus.800.3.promax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "promax")
fa.Melb.Trainbus.800.3.promax

fa.Melb.Trainbus.800.3.varimax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "varimax")
fa.Melb.Trainbus.800.3.varimax

#remove intersections, try promax again
fa.data.Melb.Trainbus.800<-Melb.Trainbus.800[,c(19, 21, 22, 25,26, 28,30, 32:33,43)]

fa.Melb.Trainbus.800.3.promax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "promax")
fa.Melb.Trainbus.800.3.promax

fa.Melb.Trainbus.800.3.varimax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "varimax")
fa.Melb.Trainbus.800.3.varimax

#remove access
fa.data.Melb.Trainbus.800<-Melb.Trainbus.800[,c(19, 21, 22, 25,26, 28,30, 32:33)]
fa.Melb.Trainbus.800.3.promax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "promax")
fa.Melb.Trainbus.800.3.promax

fa.Melb.Trainbus.800.3.varimax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "varimax")
fa.Melb.Trainbus.800.3.varimax
#remove access

fa.data.Melb.Trainbus.800<-Melb.Trainbus.800[,c(19, 21, 22, 25,26, 28,30, 33)]
fa.Melb.Trainbus.800.3.promax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "promax")
fa.Melb.Trainbus.800.3.promax

fa.Melb.Trainbus.800.3.varimax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "varimax")
fa.Melb.Trainbus.800.3.varimax

#remove housing div (cross-loading)
fa.data.Melb.Trainbus.800<-Melb.Trainbus.800[,c(19, 21, 22, 25,26, 30, 33)]
fa.Melb.Trainbus.800.3.promax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "promax")
fa.Melb.Trainbus.800.3.promax
#promax gives best solution, explaining 78.8% of variance (800m)

#600
fa.Melb.Trainbus.600.3<-factanal(fa.data.Melb.Trainbus.600, factors = 3, rotation = "none")
fa.Melb.Trainbus.600.3

fa.Melb.Trainbus.600.4<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "none")
fa.Melb.Trainbus.600.4

fa.Melb.Trainbus.600.4.promax<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "promax")
fa.Melb.Trainbus.600.4.promax

#remove urban
fa.data.Melb.Trainbus.600<-Melb.Trainbus.600[,c(19, 21, 22, 25, 26, 28:33,43)]
fa.Melb.Trainbus.600.4.promax<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "promax")
fa.Melb.Trainbus.600.4.promax

#remove intersections (high uniqueness)
fa.data.Melb.Trainbus.600<-Melb.Trainbus.600[,c(19, 21, 22, 25, 26, 28, 30:33,43)]

fa.Melb.Trainbus.600.4.promax<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "promax")
fa.Melb.Trainbus.600.4.promax

fa.Melb.Trainbus.600.4.varimax<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "varimax")
fa.Melb.Trainbus.600.4.varimax

#remove access (not loading on to any factors)

fa.data.Melb.Trainbus.600<-Melb.Trainbus.600[,c(19, 21, 22, 25:26, 28, 30:33)]
fa.Melb.Trainbus.600.4.promax<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "promax")
fa.Melb.Trainbus.600.4.promax
#promax gives best solution, explaining 74.7% of variance (600m)

#800m: promax, 3-factors (78.8% variance explained)
capture.output(fa.Melb.Trainbus.800.3.promax, file = "fa.Melb.Trainbus.800.3.promax.txt")
#600m: promax, 4-factors (74.7% of varance explained)
capture.output(fa.Melb.Trainbus.600.4.promax, file = "fa.Melb.Trainbus.600.4.promax.txt")

#rerun scree plots just to check with data subset
ev_trainbus_600_subset <- eigen(cor(fa.data.Melb.Trainbus.600))
ev_trainbus_800_subset <- eigen(cor(fa.data.Melb.Trainbus.800))# get eigenvalues
ap_trainbus_600_subset <- parallel(subject=nrow(fa.data.Melb.Trainbus.600),var=ncol(fa.data.Melb.Trainbus.600),
                            rep=100,cent=.05)
nS_trainbus_600_subset <- nScree(x=ev_trainbus_600_subset$values, aparallel=ap_trainbus_600$eigen_subset$qevpea)

ap_trainbus_800_subset <- parallel(subject=nrow(fa.data.Melb.Trainbus.800),var=ncol(fa.data.Melb.Trainbus.800),
                            rep=100,cent=.05)
nS_trainbus_800_subset <- nScree(x=ev_trainbus_800_subset$values, aparallel=ap_trainbus_800_subset$eigen$qevpea)

plotnScree(nS_trainbus_800_subset) #2/3 eigenvalues
plotnScree(nS_trainbus_600_subset) #3/4 eigenvalues

#Step 1.5 estimate factor scores and add to the master data frame
#800
Trainbus_fs_800 <- factor.scores(fa.data.Melb.Trainbus.800, fa.Melb.Trainbus.800.3.promax)     
Trainbus_fs_800 <- Trainbus_fs_800$scores            #get the columns of factor scores for each case
Melb.Trainbus.800<- cbind(Melb.Trainbus.800,Trainbus_fs_800)    #append factor scores to dataset (you can also #use merge()) or something comparable.

#600
Trainbus_fs_600 <- factor.scores(fa.data.Melb.Trainbus.600, fa.Melb.Trainbus.600.4.promax)     
Trainbus_fs_600 <- Trainbus_fs_600$scores         #get the columns of factor scores for each case
Melb.Trainbus.600<- cbind(Melb.Trainbus.600,Trainbus_fs_600)   #append factor scores to dataset (you can also #use merge()) or something comparable.


#Step 2. Estimate the covariance between dependent variables
#covariance of the logarithm
cov_Trainbus_ln<-cor.test(Melb.Trainbus.600$ln_Train, Melb.Trainbus.600$ln_Bus, method = "pearson", conf.level = 0.95)
cov_Trainbus_ln

capture.output(cov_Trainbus_ln,file="cov_Trainbus_ln.txt")

#covariance of the linear ridership
cov_Trainbus<-cor.test(Melb.Trainbus.600$Patronage_Train, Melb.Trainbus.600$Updated_Patronage_Bus, method = "pearson", conf.level = 0.95)
cov_Trainbus

capture.output(cov_Trainbus,file="cov_Trainbus.txt")

#Step 3. Check for colinearity among variables
#600m
library(car)
Melb.Trainbus.600.VIF<-vif(lm(ln_Train ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X25_Rural+X26_Access+Factor1+Factor2+Factor3+Factor4, data=Melb.Trainbus.600))
Melb.Trainbus.600.VIF
#no colinearity

Melb.Trainbus.800.VIF<-vif(lm(ln_Train ~ X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X15_DestCount+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X25_Rural+X26_Access+Factor1+Factor2+Factor3, data=Melb.Trainbus.800))
Melb.Trainbus.800.VIF
#Factor 1 vif = 5.116 --> remove

#Step 4. Get correlation matrix and check for low loading on variables
#600m
Corrdata.Trainbus.600<-Melb.Trainbus.600[,c(45, 47, 24, 27, 29, 34, 35, 37:41, 43, 48:51)]
#Option 1 for Correlation matrices with p-values
Corrdata.Trainbus.600<-rcorr(as.matrix(Corrdata.Trainbus.600))

#option 2 for flat correlation matrix
#Set up a custom function to flatten
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

options(max.print=1000000)

FlatCor.Trainbus.600<-flattenCorrMatrix(Corrdata.Trainbus.600$r,Corrdata.Trainbus.600$P)
capture.output(FlatCor.Trainbus.600,file="FlatCor.Trainbus.600.csv")

#not significant for ln_bus
#X23_MeanSize
#X26_Access
#Factor3


#ln_Train
#X23_MeanSize
#Factor2

#exclude mean size

#800m
Corrdata.Trainbus.800<-Melb.Trainbus.800[,c(45, 47, 24, 27:29, 31:32, 34, 35, 37:41, 49:50)]
#Option 1 for Correlation matrices with p-values
Corrdata.Trainbus.800<-rcorr(as.matrix(Corrdata.Trainbus.800))

#option 2 for flat correlation matrix
FlatCor.Trainbus.800<-flattenCorrMatrix(Corrdata.Trainbus.800$r,Corrdata.Trainbus.800$P)
capture.output(FlatCor.Trainbus.800,file="FlatCor.Trainbus.800.csv")

#not significant for ln_bus
#X23_MeanSize

#ln_Train
#X23_MeanSize

#do not include mean size

#step 4 Maximally adjusted models. 
#600m
Melb.Trainbus.600.MMLR.1<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor2+Factor3+Factor4, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.1)
Anova(Melb.Trainbus.600.MMLR.1)

#remove PropFTE
Melb.Trainbus.600.MMLR.2<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor2+Factor3+Factor4, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.2)
Anova(Melb.Trainbus.600.MMLR.2)

#remove factor2
Melb.Trainbus.600.MMLR.3<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor3+Factor4, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.3)
Anova(Melb.Trainbus.600.MMLR.3)

#remove access
Melb.Trainbus.600.MMLR.4<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X25_Rural+Factor1+Factor3+Factor4, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.4)
Anova(Melb.Trainbus.600.MMLR.4)

#remove factor4
Melb.Trainbus.600.MMLR.5<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X25_Rural+Factor1+Factor3, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.5)
Anova(Melb.Trainbus.600.MMLR.5)

#remove rural
Melb.Trainbus.600.MMLR.6<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor1+Factor3, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.6)
Anova(Melb.Trainbus.600.MMLR.6)

#remove factor 1
Melb.Trainbus.600.MMLR.7<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.7)
Anova(Melb.Trainbus.600.MMLR.7)

#remove intersectins
Melb.Trainbus.600.MMLR.8<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.8)
Anova(Melb.Trainbus.600.MMLR.8)

#remove entropy
Melb.Trainbus.600.MMLR.9<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600)
summary(Melb.Trainbus.600.MMLR.9)
Anova(Melb.Trainbus.600.MMLR.9)

#run diagnostics
par(mfrow=c(2,2))
plot(lm(ln_Train ~ X7_PropComm+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600))

#influential outliers: 1-600-C; 2-600-C; 3-600-C; 37-600-C
which(rownames(Melb.Trainbus.600) == "1-600-C") #135
which(rownames(Melb.Trainbus.600) == "2-600-C") #134
which(rownames(Melb.Trainbus.600) == "3-600-C") #133
which(rownames(Melb.Trainbus.600) == "37-600-C")#97

plot(lm(ln_Bus ~ X7_PropComm+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600))
#approaching influential outliers: 14-600-C;
which(rownames(Melb.Trainbus.600) == "14-600-C") #39
#remove influential outlier
Melb.Trainbus.600.rd2 <- Melb.Trainbus.600[-c(39, 97, 133, 134, 135),]

#rerun maximally adjusted model

Melb.Trainbus.600.MMLR.2.1<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor2+Factor3+Factor4, data=Melb.Trainbus.600.rd2)
summary(Melb.Trainbus.600.MMLR.2.1)
Anova(Melb.Trainbus.600.MMLR.2.1)

#remove factor1
Melb.Trainbus.600.MMLR.2.2<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor2+Factor3+Factor4, data=Melb.Trainbus.600.rd2)
summary(Melb.Trainbus.600.MMLR.2.2)
Anova(Melb.Trainbus.600.MMLR.2.2)

#remove factor 2
Melb.Trainbus.600.MMLR.2.3<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor3+Factor4, data=Melb.Trainbus.600.rd2)
summary(Melb.Trainbus.600.MMLR.2.3)
Anova(Melb.Trainbus.600.MMLR.2.3)

#remove AcCount
Melb.Trainbus.600.MMLR.2.4<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor3+Factor4, data=Melb.Trainbus.600.rd2)
summary(Melb.Trainbus.600.MMLR.2.4)
Anova(Melb.Trainbus.600.MMLR.2.4)

#remove access
Melb.Trainbus.600.MMLR.2.5<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+Factor3+Factor4, data=Melb.Trainbus.600.rd2)
summary(Melb.Trainbus.600.MMLR.2.5)
Anova(Melb.Trainbus.600.MMLR.2.5)

#remove factor 4
Melb.Trainbus.600.MMLR.2.6<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.600.rd2)
summary(Melb.Trainbus.600.MMLR.2.6)
Anova(Melb.Trainbus.600.MMLR.2.6)

#remove propFTE
Melb.Trainbus.600.MMLR.2.7<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.600.rd2)
summary(Melb.Trainbus.600.MMLR.2.7)
Anova(Melb.Trainbus.600.MMLR.2.7)

#remove MedInc
Melb.Trainbus.600.MMLR.2.8<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X18_ACNear+X20_LOS+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.600.rd2)
summary(Melb.Trainbus.600.MMLR.2.8)
Anova(Melb.Trainbus.600.MMLR.2.8)

plot(lm(ln_Train ~ X7_PropComm+X10_Entropy+X12_Intersections+X18_ACNear+X20_LOS+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.600.rd2))
#influential outliers: 163-600_C; 39-600-C

plot(lm(ln_Bus ~ X7_PropComm+X10_Entropy+X12_Intersections+X18_ACNear+X20_LOS+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.600.rd2))
#approaching influential outliers: 104-600-C; 39-600-C; 17-600-C

which(rownames(Melb.Trainbus.600.rd2) == "163-600-C") #130
which(rownames(Melb.Trainbus.600.rd2) == "39-600-C") #129
which(rownames(Melb.Trainbus.600.rd2) == "104-600-C") #8
which(rownames(Melb.Trainbus.600.rd2) == "17-600-C") #61

#remove influential outlier
Melb.Trainbus.600.rd3 <- Melb.Trainbus.600.rd2[-c(130, 129, 8, 61),]

#round 3 maximally adjusted model
#600m
Melb.Trainbus.600.MMLR.3.1<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor2+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.1)
Anova(Melb.Trainbus.600.MMLR.3.1)

#remove factor 2
Melb.Trainbus.600.MMLR.3.2<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.2)
Anova(Melb.Trainbus.600.MMLR.3.2)

#remove rural
Melb.Trainbus.600.MMLR.3.3<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X26_Access+Factor1+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.3)
Anova(Melb.Trainbus.600.MMLR.3.3)

#remove PropFTE
Melb.Trainbus.600.MMLR.3.4<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor1+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.4)
Anova(Melb.Trainbus.600.MMLR.3.4)

#remove factor 1
Melb.Trainbus.600.MMLR.3.5<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.5)
Anova(Melb.Trainbus.600.MMLR.3.5)

#remove access
Melb.Trainbus.600.MMLR.3.6<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.6)
Anova(Melb.Trainbus.600.MMLR.3.6)

#remove propcomm
Melb.Trainbus.600.MMLR.3.7<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.7)
Anova(Melb.Trainbus.600.MMLR.3.7)

#remove intersections
Melb.Trainbus.600.MMLR.3.8<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.8)
Anova(Melb.Trainbus.600.MMLR.3.8)

#remove ACNear
Melb.Trainbus.600.MMLR.3.9<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X17_ACCount+X20_LOS+X22_MedInc+X24_Urban+Factor3+Factor4, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.9)
Anova(Melb.Trainbus.600.MMLR.3.9)

#remove factor 4
Melb.Trainbus.600.MMLR.3.10<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X17_ACCount+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.10)
Anova(Melb.Trainbus.600.MMLR.3.10)

Melb.Trainbus.600.MMLR.3.10<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X17_ACCount+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600.rd3)
summary(Melb.Trainbus.600.MMLR.3.10)
Anova(Melb.Trainbus.600.MMLR.3.10)

#run diagnostics
plot(lm(ln_Train ~ X10_Entropy+X17_ACCount+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600.rd3))

plot(lm(ln_Bus ~ X10_Entropy+X17_ACCount+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600.rd3))

#approaching cook's distance: 15-600-C
which(rownames(Melb.Trainbus.600.rd3) == "15-600-C") #61

#remove influential outlier
Melb.Trainbus.600.rd4 <- Melb.Trainbus.600.rd3[-c(117),]

#rerun maximally adjusted model

Melb.Trainbus.600.MMLR.4.1<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor2+Factor3+Factor4, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.1)
Anova(Melb.Trainbus.600.MMLR.4.1)

#remove Factor 2

Melb.Trainbus.600.MMLR.4.2<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor3+Factor4, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.2)
Anova(Melb.Trainbus.600.MMLR.4.2)

#reove propFTE
Melb.Trainbus.600.MMLR.4.3<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor1+Factor3+Factor4, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.3)
Anova(Melb.Trainbus.600.MMLR.4.3)

#remove factor 1
Melb.Trainbus.600.MMLR.4.4<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor3+Factor4, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.4)
Anova(Melb.Trainbus.600.MMLR.4.4)

#remove rural
Melb.Trainbus.600.MMLR.4.5<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor3+Factor4, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.5)
Anova(Melb.Trainbus.600.MMLR.4.5)

#remove access
Melb.Trainbus.600.MMLR.4.6<-lm(cbind(ln_Train, ln_Bus) ~ X7_PropComm+X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3+Factor4, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.6)
Anova(Melb.Trainbus.600.MMLR.4.6)

#remove propr comm
Melb.Trainbus.600.MMLR.4.7<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+Factor3+Factor4, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.7)
Anova(Melb.Trainbus.600.MMLR.4.7)

#remove ACNear
Melb.Trainbus.600.MMLR.4.8<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X12_Intersections+X17_ACCount+X20_LOS+X22_MedInc+X24_Urban+Factor3+Factor4, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.8)
Anova(Melb.Trainbus.600.MMLR.4.8)

#remove factor 4
Melb.Trainbus.600.MMLR.4.9<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X12_Intersections+X17_ACCount+X20_LOS+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.9)
Anova(Melb.Trainbus.600.MMLR.4.9)

#remove medinc
Melb.Trainbus.600.MMLR.4.10<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X12_Intersections+X17_ACCount+X20_LOS+X24_Urban+Factor3, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.10)
Anova(Melb.Trainbus.600.MMLR.4.10)

#remove urban
Melb.Trainbus.600.MMLR.4.11<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X12_Intersections+X17_ACCount+X20_LOS+Factor3, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.11)
Anova(Melb.Trainbus.600.MMLR.4.11)

#remove intersections
Melb.Trainbus.600.MMLR.4.12<-lm(cbind(ln_Train, ln_Bus) ~ X10_Entropy+X17_ACCount+X20_LOS+Factor3, data=Melb.Trainbus.600.rd4)
summary(Melb.Trainbus.600.MMLR.4.12)
Anova(Melb.Trainbus.600.MMLR.4.12)

plot(lm(ln_Bus ~ X10_Entropy+X17_ACCount+X20_LOS+Factor3, data=Melb.Trainbus.600.rd4))

plot(lm(ln_Train ~ X10_Entropy+X17_ACCount+X20_LOS+Factor3, data=Melb.Trainbus.600.rd4))

capture.output(summary(Melb.Trainbus.600.MMLR.4.11),file = "Melb.trainbus.600.MMLR.4.11.csv")

#800m
Melb.Trainbus.800.MMLR.1 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X15_DestCount+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.1)
Anova(Melb.Trainbus.800.MMLR.1)

#removeDestCount
Melb.Trainbus.800.MMLR.2 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.2)
Anova(Melb.Trainbus.800.MMLR.2)

#remove MedInc
Melb.Trainbus.800.MMLR.3 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+X25_Rural+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.3)
Anova(Melb.Trainbus.800.MMLR.3)

#remove rural
Melb.Trainbus.800.MMLR.4<- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.4)
Anova(Melb.Trainbus.800.MMLR.4)

#remove entropy
Melb.Trainbus.800.MMLR.5<- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.5)
Anova(Melb.Trainbus.800.MMLR.5)

#remove destscore
Melb.Trainbus.800.MMLR.6<- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.6)
Anova(Melb.Trainbus.800.MMLR.6)

#remove access
Melb.Trainbus.800.MMLR.7<- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+Factor2+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.7)
Anova(Melb.Trainbus.800.MMLR.7)

#remove Factor 2
Melb.Trainbus.800.MMLR.8<- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.8)
Anova(Melb.Trainbus.800.MMLR.8)

#remove intersections
Melb.Trainbus.800.MMLR.9<- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.9)
Anova(Melb.Trainbus.800.MMLR.9)

#remove PropFTE
Melb.Trainbus.800.MMLR.10<- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X17_ACCount+X18_ACNear+X20_LOS+X24_Urban+Factor3, data=Melb.Trainbus.800)
summary(Melb.Trainbus.800.MMLR.10)
Anova(Melb.Trainbus.800.MMLR.10)

#run diagnostics
par(mfrow=c(2,2))
plot(lm(ln_Train ~ X7_PropComm+X11_HousingDiv+X17_ACCount+X18_ACNear+X20_LOS+X24_Urban+Factor3, data=Melb.Trainbus.800))

#influential outliers: 1-800-C; 2-800-C; 3-800-C

plot(lm(ln_Bus ~ X7_PropComm+X11_HousingDiv+X17_ACCount+X18_ACNear+X20_LOS+X24_Urban+Factor3, data=Melb.Trainbus.800))
#(Approaching influential outliers: 14-800-C, 39-800-C; 15-800-C

which(rownames(Melb.Trainbus.800) == "1-800-C") #135
which(rownames(Melb.Trainbus.800) == "2-800-C") #134
which(rownames(Melb.Trainbus.800) == "3-800-C") #133
which(rownames(Melb.Trainbus.800) == "14-800-C")#39
which(rownames(Melb.Trainbus.800) == "39-800-C")#131
which(rownames(Melb.Trainbus.800) == "15-800-C")#121

#remove influential outlier
Melb.Trainbus.800.rd2 <- Melb.Trainbus.800[-c(39, 131, 121, 133, 134, 135),]

#rerun maximally adjusted model
Melb.Trainbus.800.MMLR.2.1 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X15_DestCount+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800.rd2)
summary(Melb.Trainbus.800.MMLR.2.1)
Anova(Melb.Trainbus.800.MMLR.2.1)

#remove factor 2
Melb.Trainbus.800.MMLR.2.2 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X15_DestCount+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor3, data=Melb.Trainbus.800.rd2)
summary(Melb.Trainbus.800.MMLR.2.2)
Anova(Melb.Trainbus.800.MMLR.2.2)

#remove destcount
Melb.Trainbus.800.MMLR.2.3 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor3, data=Melb.Trainbus.800.rd2)
summary(Melb.Trainbus.800.MMLR.2.3)
Anova(Melb.Trainbus.800.MMLR.2.3)

#remove entropy
Melb.Trainbus.800.MMLR.2.4 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor3, data=Melb.Trainbus.800.rd2)
summary(Melb.Trainbus.800.MMLR.2.4)
Anova(Melb.Trainbus.800.MMLR.2.4)

#remove access
Melb.Trainbus.800.MMLR.2.5 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.800.rd2)
summary(Melb.Trainbus.800.MMLR.2.5)
Anova(Melb.Trainbus.800.MMLR.2.5)

#remove MedINc
Melb.Trainbus.800.MMLR.2.6 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.800.rd2)
summary(Melb.Trainbus.800.MMLR.2.6)
Anova(Melb.Trainbus.800.MMLR.2.6)

#remove ACNear
Melb.Trainbus.800.MMLR.2.7 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X20_LOS+X21_PropFTE+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.800.rd2)
summary(Melb.Trainbus.800.MMLR.2.7)
Anova(Melb.Trainbus.800.MMLR.2.7)

#run diagnostics
plot(lm(ln_Train ~ X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X20_LOS+X21_PropFTE+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.800.rd2))

#influential outlier: 163-800-C

plot(lm(ln_Bus ~ X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X20_LOS+X21_PropFTE+X24_Urban+X25_Rural+Factor3, data=Melb.Trainbus.800.rd2))
#influential outlier: 163-800-C

#Accept solution
which(rownames(Melb.Trainbus.800.rd2) == "163-800-C")#129

#remove influential outlier
Melb.Trainbus.800.rd3 <- Melb.Trainbus.800.rd2[-c(129),]

#rerun maximally adjusted model
Melb.Trainbus.800.MMLR.3.1 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X15_DestCount+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X25_Rural+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800.rd3)
summary(Melb.Trainbus.800.MMLR.3.1)
Anova(Melb.Trainbus.800.MMLR.3.1)

#remove rural
Melb.Trainbus.800.MMLR.3.2 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X15_DestCount+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X26_Access+Factor2+Factor3, data=Melb.Trainbus.800.rd3)
summary(Melb.Trainbus.800.MMLR.3.2)
Anova(Melb.Trainbus.800.MMLR.3.2)

#remove factor 2
Melb.Trainbus.800.MMLR.3.3 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X15_DestCount+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X26_Access+Factor3, data=Melb.Trainbus.800.rd3)
summary(Melb.Trainbus.800.MMLR.3.3)
Anova(Melb.Trainbus.800.MMLR.3.3)

#remove dest count
Melb.Trainbus.800.MMLR.3.4 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X10_Entropy+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X26_Access+Factor3, data=Melb.Trainbus.800.rd3)
summary(Melb.Trainbus.800.MMLR.3.4)
Anova(Melb.Trainbus.800.MMLR.3.4)

#remove entropy
Melb.Trainbus.800.MMLR.3.5 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+X26_Access+Factor3, data=Melb.Trainbus.800.rd3)
summary(Melb.Trainbus.800.MMLR.3.5)
Anova(Melb.Trainbus.800.MMLR.3.5)

#remove access
Melb.Trainbus.800.MMLR.3.6 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X24_Urban+Factor3, data=Melb.Trainbus.800.rd3)
summary(Melb.Trainbus.800.MMLR.3.6)
Anova(Melb.Trainbus.800.MMLR.3.6)

#remove MedInc
Melb.Trainbus.800.MMLR.3.7 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X24_Urban+Factor3, data=Melb.Trainbus.800.rd3)
summary(Melb.Trainbus.800.MMLR.3.7)
Anova(Melb.Trainbus.800.MMLR.3.7)

#remove ACNear
Melb.Trainbus.800.MMLR.3.8 <- lm(cbind(ln_Train, ln_Bus) ~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X20_LOS+X21_PropFTE+X24_Urban+Factor3, data=Melb.Trainbus.800.rd3)
summary(Melb.Trainbus.800.MMLR.3.8)
Anova(Melb.Trainbus.800.MMLR.3.8)

#run diagnostics
plot(lm(ln_Bus~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X20_LOS+X21_PropFTE+X24_Urban+Factor3, data=Melb.Trainbus.800.rd3))

plot(lm(ln_Train~X7_PropComm+X11_HousingDiv+X12_Intersections+X14_DestScore+X17_ACCount+X20_LOS+X21_PropFTE+X24_Urban+Factor3, data=Melb.Trainbus.800.rd3))

#accept solution 
capture.output(summary(Melb.Trainbus.800.MMLR.3.8), file = "Melb.Trainbus.800.MMLR.3.8.csv")
