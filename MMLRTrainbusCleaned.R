# Title:                  Melbourne Train-bus MMLR (cleaned version)
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
#1) Radius 800m and 600m catchments were tested; 800m catchment yielded higher explanatory power solution for both modes and will be used to proceed. 
#2) Outliers observed to affect results for both an 800m and 600 m catchment are removed: Sample ID: 1, 2, 3

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

#subsetting for mode in the census of all eligible stops

Melb.Trainbus.800<- MMLR_Data[ which(MMLR_Data$Type=='Trainbus'
                                     & MMLR_Data$Radius=='800'
                                     & MMLR_Data$Set_Sample_ID == 'C'),]

#Get row IDs for the outliers. 
which(rownames(Melb.Trainbus.800) == "1-800-C")#135
which(rownames(Melb.Trainbus.800) == "2-800-C")#134
which(rownames(Melb.Trainbus.800) == "3-800-C")#133

#remove outliers and create new file 'NO' ('N'o 'O'utliers)
Melb.Trainbus.800.NO <- Melb.Trainbus.800[-c(133, 134, 135),]

#Assumption 1:  'Set ID' refers to the method of defining the 'walk catchment'. 'D' (distributed) means     catchments were separately calculated for the tram and bus points, then merged. 'C' (centroid) means the geographic centroid of the points was first found, and the catchment estimated from there. Data corresponding to the 'centroid' (C) catchments is used in this analysis. 

#Variables to include in the factor analysis (columns 19, 21,22, 24 - 43)
#X2_EmpDen+X4_PopDen+X5_Dwelling+X7_PropComm+X8_RetailEmp+X9_Balance+X10_Entropy+X11_HousingDiv+X12_Intersections+X13_PBN+X14_DestScore+X15_DestCount+X16_DistCBD+X17_ACCount+X18_ACNear+X19_FTZ+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X25_Rural+X26_Access

#Step 1.1: form a data frame that comprises only the variables to be included in the factor analysis (built environment)
#include variables measured on a count or continuous scale. This means "FTZ" and "AC_Count" should be excluded. Also "rural" is 0 for almost all responses, so exclude
fa.data.Melb.Trainbus.800.NO<-Melb.Trainbus.800.NO[,c(19, 21, 22, 24:33,35,41,43)]

#Step 1.2: Specify number of factors. Based on theory, will try four or five factors, consituting 1) density 2) diversity 3) design 4) regional accessibility and 5) local accessibility/walkability (Ewing and Cevero 2010, Voulgaris et al. 2017)

#alternatively, check scree plot
install.packages("nFactors")
library(nFactors)
ev_trainbus_800.NO <- eigen(cor(fa.data.Melb.Trainbus.800.NO))# get eigenvalues
ap_trainbus_800.NO <- parallel(subject=nrow(fa.data.Melb.Trainbus.800.NO),var=ncol(fa.data.Melb.Trainbus.800.NO),
                            rep=100,cent=.05)
nS_trainbus_800.NO <- nScree(x=ev_trainbus_800.NO$values, aparallel=ap_trainbus_800.NO$eigen$qevpea)
plotnScree(nS_trainbus_800.NO) #3 eigenvalues

#Step 1.3 Run factor analysis
fa.Melb.Trainbus.800.NO.3<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 3, rotation = "none")
fa.Melb.Trainbus.800.NO.3

#remove high uniqueness variables: DestScore, ACNear ,Urban

fa.data.Melb.Trainbus.800.NO<-Melb.Trainbus.800.NO[,c(19, 21, 22, 24:30, 32:33,43)]

#rerun scree
ev_trainbus_800.NO <- eigen(cor(fa.data.Melb.Trainbus.800.NO))# get eigenvalues
ap_trainbus_800.NO <- parallel(subject=nrow(fa.data.Melb.Trainbus.800.NO),var=ncol(fa.data.Melb.Trainbus.800.NO),
                               rep=100,cent=.05)
nS_trainbus_800.NO <- nScree(x=ev_trainbus_800.NO$values, aparallel=ap_trainbus_800.NO$eigen$qevpea)
plotnScree(nS_trainbus_800.NO)

fa.Melb.Trainbus.800.NO.3<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 3, rotation = "none")
fa.Melb.Trainbus.800.NO.3

fa.Melb.Trainbus.800.NO.3.promax<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 3, rotation = "promax")
fa.Melb.Trainbus.800.NO.3.promax

fa.Melb.Trainbus.800.NO.3.varimax<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 3, rotation = "varimax")
fa.Melb.Trainbus.800.NO.3.varimax

#try 4-factor solution
fa.Melb.Trainbus.800.NO.4<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 4, rotation = "none")
fa.Melb.Trainbus.800.NO.4

fa.Melb.Trainbus.800.NO.4.promax<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 4, rotation = "promax")
fa.Melb.Trainbus.800.NO.4.promax

fa.Melb.Trainbus.800.NO.4.varimax<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 4, rotation = "varimax")
fa.Melb.Trainbus.800.NO.4.varimax

#four factor solution eliminates cross loading and is more interpretable than the three factor solutions
#try removing access
fa.data.Melb.Trainbus.800.NO<-Melb.Trainbus.800.NO[,c(19, 21, 22, 24:30, 32:33)]
fa.Melb.Trainbus.800.NO.4.promax<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 4, rotation = "promax")
fa.Melb.Trainbus.800.NO.4.promax

fa.Melb.Trainbus.800.NO.4.varimax<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 4, rotation = "varimax")
fa.Melb.Trainbus.800.NO.4.varimax

#remove PBN, specify 3 eigenvalues
fa.data.Melb.Trainbus.800.NO<-Melb.Trainbus.800.NO[,c(19, 21, 22, 24:29, 32:33)]
fa.Melb.Trainbus.800.NO.3.promax<-factanal(fa.data.Melb.Trainbus.800.NO, factors = 3, rotation = "promax")
fa.Melb.Trainbus.800.NO.3.promax

capture.output(fa.Melb.Trainbus.800.NO.3.promax,file="fa.Melb.Trainbus.800.NO.3.promax.txt")

#very intepretable, parsimonious solution (agrees with factors identified for clustering)

#chec scree plot again
ev_trainbus_800.NO <- eigen(cor(fa.data.Melb.Trainbus.800.NO))# get eigenvalues
ap_trainbus_800.NO <- parallel(subject=nrow(fa.data.Melb.Trainbus.800.NO),var=ncol(fa.data.Melb.Trainbus.800.NO),
                               rep=100,cent=.05)
nS_trainbus_800.NO <- nScree(x=ev_trainbus_800.NO$values, aparallel=ap_trainbus_800.NO$eigen$qevpea)
plotnScree(nS_trainbus_800.NO) #3 eigenvalues

#Step 1.5 estimate factor scores and add to the master data frame
Trainbus_fs_800.NO <- factor.scores(fa.data.Melb.Trainbus.800.NO, fa.Melb.Trainbus.800.NO.3.promax)
Trainbus_fs_800.NO <- Trainbus_fs_800.NO$scores           #get the columns of factor scores for each case
Melb.Trainbus.800.NO<- cbind(Melb.Trainbus.800.NO,Trainbus_fs_800.NO) #append factor scores to dataset (you can also #use merge()) or something comparable.

#Step 2 estimate covariance of the outcome variables 
cov_Trainbus_ln<-cor.test(Melb.Trainbus.800.NO$ln_Train, Melb.Trainbus.800.NO$ln_Bus, method = "pearson", conf.level = 0.95)
cov_Trainbus_ln

capture.output(cov_Trainbus_ln,file="cov_Trainbus_ln.NO.txt")

#covariance of the linear ridership
cov_Trainbus<-cor.test(Melb.Trainbus.800.NO$Patronage_Train, Melb.Trainbus.800.NO$Updated_Patronage_Bus, method = "pearson", conf.level = 0.95)
cov_Trainbus

capture.output(cov_Trainbus,file="cov_Trainbus.NO.txt")

#step 3 Check for multicolinearity
Melb.Trainbus.800.NO.VIF<-vif(lm(ln_Train ~ X13_PBN+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X26_Access+Factor1+Factor2+Factor3, data =Melb.Trainbus.800.NO))
Melb.Trainbus.800.NO.VIF

#no multicolinearity

#step 4 Simple correlations
Corrdata.Trainbus.800.NO<-Melb.Trainbus.800.NO[,c(45, 47, 30, 31, 34, 35, 37:41, 43, 48:50)]
#Option 1 for Correlation matrices with p-values
Corrdata.Trainbus.800.NO<-rcorr(as.matrix(Corrdata.Trainbus.800.NO))

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

FlatCor.Trainbus.800.NO<-flattenCorrMatrix(Corrdata.Trainbus.800.NO$r,Corrdata.Trainbus.800.NO$P)
capture.output(FlatCor.Trainbus.800.NO,file="FlatCor.Trainbus.800.NO.csv")

#not significant for ln_bus
#X21_PropFTE
#X23_MeanSize
#X26_Access
#Factor1

#ln_Train
#X21_PropFTE
#X23_MeanSize
#Factor2

#exclude mean size and prop FTE

#Step 4 maximally adjusted model
Melb.Trainbus.800.NO.MMLR.1<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor1+Factor2+Factor3, data =Melb.Trainbus.800.NO)
summary(Melb.Trainbus.800.NO.MMLR.1)
Anova(Melb.Trainbus.800.NO.MMLR.1)

#remove DestScore
Melb.Trainbus.800.NO.MMLR.2<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor1+Factor2+Factor3, data =Melb.Trainbus.800.NO)
summary(Melb.Trainbus.800.NO.MMLR.2)
Anova(Melb.Trainbus.800.NO.MMLR.2)

#remove factor 3
Melb.Trainbus.800.NO.MMLR.3<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor1+Factor2, data =Melb.Trainbus.800.NO)
summary(Melb.Trainbus.800.NO.MMLR.3)
Anova(Melb.Trainbus.800.NO.MMLR.3)

#remove MedInc
Melb.Trainbus.800.NO.MMLR.4<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X24_Urban+X26_Access+Factor1+Factor2, data =Melb.Trainbus.800.NO)
summary(Melb.Trainbus.800.NO.MMLR.4)
Anova(Melb.Trainbus.800.NO.MMLR.4)

#remove Access
Melb.Trainbus.800.NO.MMLR.5<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X24_Urban+Factor1+Factor2, data =Melb.Trainbus.800.NO)
summary(Melb.Trainbus.800.NO.MMLR.5)
Anova(Melb.Trainbus.800.NO.MMLR.5)

#remove factor 2
Melb.Trainbus.800.NO.MMLR.6<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X24_Urban+Factor1, data =Melb.Trainbus.800.NO)
summary(Melb.Trainbus.800.NO.MMLR.6)
Anova(Melb.Trainbus.800.NO.MMLR.6)

#remove ACNear
Melb.Trainbus.800.NO.MMLR.7<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X20_LOS+X24_Urban+Factor1, data =Melb.Trainbus.800.NO)
summary(Melb.Trainbus.800.NO.MMLR.7)
Anova(Melb.Trainbus.800.NO.MMLR.7)

#run diagnostics
par(mfrow=c(2,2))
plot(lm(ln_Train ~ X13_PBN+X17_ACCount+X20_LOS+X24_Urban+Factor1, data =Melb.Trainbus.800.NO))
#consider removing 163-800-C;39-800-C;104-800-C

plot(lm(ln_Bus ~ X13_PBN+X17_ACCount+X20_LOS+X24_Urban+Factor1, data =Melb.Trainbus.800.NO))
#consider removing (also) 104-800-C; 15-800-C and 14-800-C

which(rownames(Melb.Trainbus.800.NO) == "163-800-C") #132
which(rownames(Melb.Trainbus.800.NO) == "39-800-C") #131
which(rownames(Melb.Trainbus.800.NO) == "104-800-C") #8
which(rownames(Melb.Trainbus.800.NO) == "14-800-C") #39
which(rownames(Melb.Trainbus.800.NO) == "15-800-C") #121

#remove influential outlier
Melb.Trainbus.800.NO.rd2 <- Melb.Trainbus.800.NO[-c(8, 39, 121, 131, 132),]

#repeat maximally adjusted model
Melb.Trainbus.800.NO.MMLR.2.1<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor1+Factor2+Factor3, data =Melb.Trainbus.800.NO.rd2)
summary(Melb.Trainbus.800.NO.MMLR.2.1)
Anova(Melb.Trainbus.800.NO.MMLR.2.1)

#remove Factor 3
Melb.Trainbus.800.NO.MMLR.2.2<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor1+Factor2, data =Melb.Trainbus.800.NO.rd2)
summary(Melb.Trainbus.800.NO.MMLR.2.2)
Anova(Melb.Trainbus.800.NO.MMLR.2.2)

#remove Factor 1
Melb.Trainbus.800.NO.MMLR.2.3<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X14_DestScore+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor2, data =Melb.Trainbus.800.NO.rd2)
summary(Melb.Trainbus.800.NO.MMLR.2.3)
Anova(Melb.Trainbus.800.NO.MMLR.2.3)

#remove DestScore
Melb.Trainbus.800.NO.MMLR.2.4<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X24_Urban+X26_Access+Factor2, data =Melb.Trainbus.800.NO.rd2)
summary(Melb.Trainbus.800.NO.MMLR.2.4)
Anova(Melb.Trainbus.800.NO.MMLR.2.4)

#remove MedInc
#remove DestScore
Melb.Trainbus.800.NO.MMLR.2.5<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X24_Urban+X26_Access+Factor2, data =Melb.Trainbus.800.NO.rd2)
summary(Melb.Trainbus.800.NO.MMLR.2.5)
Anova(Melb.Trainbus.800.NO.MMLR.2.5)

#remove Urban
Melb.Trainbus.800.NO.MMLR.2.6<-lm(cbind(ln_Train, ln_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X26_Access+Factor2, data =Melb.Trainbus.800.NO.rd2)
summary(Melb.Trainbus.800.NO.MMLR.2.6)
Anova(Melb.Trainbus.800.NO.MMLR.2.6)

#run diagnostics
plot(lm(ln_Train ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X26_Access+Factor2, data =Melb.Trainbus.800.NO.rd2))
#no distinctive patterns

par(mfrow=c(2,2))

plot(lm(ln_Bus ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X26_Access+Factor2, data =Melb.Trainbus.800.NO.rd2))
#64=800-C approaching cook's distance. Inspect in Arc along with others, but probably leave as is

#Get Standardized regression coefficients
Trainbus_800_NO_train<-lm(ln_Train ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X26_Access+Factor2, data =Melb.Trainbus.800.NO.rd2)
Trainbus_800_NO_train<-lm.beta(Trainbus_800_NO_train)
summary(Trainbus_800_NO_train)

Trainbus_800_NO_bus<-lm(ln_Bus ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X26_Access+Factor2, data =Melb.Trainbus.800.NO.rd2)
Trainbus_800_NO_bus<-lm.beta(Trainbus_800_NO_bus)
summary(Trainbus_800_NO_bus)

capture.output(summary(Melb.Trainbus.800.NO.MMLR.2.6),file = "Melb.Trainbus.800.NO.MMLR.2.6.csv")
capture.output(summary(Trainbus_800_NO_bus),file = "Melb.Trainbus.800.bus.csv")
capture.output(summary(Trainbus_800_NO_train),file = "Melb.Trainbus.800.train.csv")

#intepreting diagnostic plots: https://data.library.virginia.edu/diagnostic-plots/