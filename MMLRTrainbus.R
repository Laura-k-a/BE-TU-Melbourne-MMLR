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
fa.data.Melb.Trainbus.800<-Melb.Trainbus.800[,c(19, 21, 22, 24:26, 28:33,35,41,43)]
fa.data.Melb.Trainbus.600<-Melb.Trainbus.600[,c(19, 21, 22, 24:26, 28:33,35,41,43)]

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
#high uniqueness: proportion commerical, PBN, DestScore, ACNear ,Urban

#600
fa.Melb.Trainbus.600.3<-factanal(fa.data.Melb.Trainbus.600, factors = 3, rotation = "none")
fa.Melb.Trainbus.600.3
#high uniqueness: proportion commerical, intersections, PBN, DistCBD, ACNear, Urban

fa.Melb.Trainbus.600.4<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "none")
fa.Melb.Trainbus.600.4
#high uniqueness: proportion commerical, ACNear, urban

#Step 1.4 Remove variables with high uniqueness values: proportion commerical, ACNear, urban
fa.data.Melb.Trainbus.800<-Melb.Trainbus.800[,c(19, 21, 22, 25:26, 28, 29, 32:33,43)]
fa.data.Melb.Trainbus.600<-Melb.Trainbus.600[,c(19, 21, 22, 25:26, 28:33,43)]


fa.Melb.Trainbus.800.3<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "none")
fa.Melb.Trainbus.800.3

fa.Melb.Trainbus.800.3.promax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "promax")
fa.Melb.Trainbus.800.3.promax

fa.Melb.Trainbus.800.3.varimax<-factanal(fa.data.Melb.Trainbus.800, factors = 3, rotation = "varimax")
fa.Melb.Trainbus.800.3.varimax

#promax solution is best (cross-loading with varimax solution)

fa.Melb.Trainbus.600.3<-factanal(fa.data.Melb.Trainbus.600, factors = 3, rotation = "none")
fa.Melb.Trainbus.600.3

fa.Melb.Trainbus.600.4<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "none")
fa.Melb.Trainbus.600.4

fa.Melb.Trainbus.600.4.promax<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "promax")
fa.Melb.Trainbus.600.4.promax

fa.Melb.Trainbus.600.4.varimax<-factanal(fa.data.Melb.Trainbus.600, factors = 4, rotation = "varimax")
fa.Melb.Trainbus.600.4.varimax

#promax gives best solution, explaining 71.6% of variance (600m)

#800m: promax, 3-factors (62.6% variance explained)
capture.output(fa.Melb.Trainbus.800.3.promax, file = "fa.Melb.Trainbus.800.3.promax.txt")
#600m: promax, 4-factors (71.65% of varance explained)
capture.output(fa.Melb.Trainbus.600.4.promax, file = "fa.Melb.Trainbus.600.4.promax.txt")

#Step 1.5 Evaluate the adequacy of the number of factors
#Methods for evaluating the appropriateness of the solution include (O'Hair 2014, pp. 106-109):
# important that their is a conceptual explanation for the identified factors
# a prior criterion: when number of factors is pre-specified
# cumulative variance >0.6

#sufficient variance explained and factors make sense. However, factor 3 does not have any variables loaded on it. Try rotation

fa.Melb.Trambus.600.noFTZ.promax<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "promax")
fa.Melb.Trambus.600.noFTZ.promax
#cross-loading of activity density

#try varimax
fa.Melb.Trambus.600.noFTZ.promax<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "varimax")
fa.Melb.Trambus.600.noFTZ.varimax

#activity density is still cross-loading, so try removing it. Substitute employment and population density, since they will no longer be singluar. 

fa.data.Melb.Trambus.600.noFTZ<-Melb.Trambus.600.noFTZ[,c(19, 21:22, 25, 28:29, 31:33,43)]
fa.data.Melb.Trambus.600.noFTZ<-as.matrix(fa.data.Melb.Trambus.600.noFTZ)
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "none")
fa.Melb.Trambus.600.noFTZ

fa.Melb.Trambus.600.noFTZ.promax<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "promax")
fa.Melb.Trambus.600.noFTZ.promax

fa.Melb.Trambus.600.noFTZ.varimax<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "varimax")
fa.Melb.Trambus.600.noFTZ.varimax

#varimax yields best solution - exlain 74.2% of variance. 

#note that the null hypothesis that 3 factors is sufficient is rejected. However, hypothesis testing has less significance in factor analysis than intepretability. 

capture.output(fa.Melb.Trambus.600.noFTZ.varimax,file ="fa.Melb.Trambus.600.noFTZ.varimax.csv")

#Also tried a 2-factor solutin (due to scree plot suggesting only 2 factors had egenvalues >1), however solution is less interpretable than the 3-factor solution, so will preference the 3-factor solution over 2-factor solution despite low loading on third factor. 

#Step 1.5 estimate factor scores and add to the master data frame
fa.Melb.Trambus.600.noFTZ.varimax<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "varimax", scores = "regression")
fa.Melb.Trambus.600.noFTZ.varimax
head(fa.Melb.Trambus.600.noFTZ.varimax$scores)

#obtain factor scores
#600m radius catchment
Trambus_fs_600 <- factor.scores(fa.data.Melb.Trambus.600.noFTZ, fa.Melb.Trambus.600.noFTZ.varimax)     
Trambus_fs_600 <- Trambus_fs_600$scores                 #get the columns of factor scores for each case
Melb.Trambus.600.noFTZ<- cbind(Melb.Trambus.600.noFTZ,Trambus_fs_600)    #append factor scores to dataset (you can also 
#use merge()) or something comparable.