# Title:                  Melbourne Tram-bus MMLR (cleaned version)
# Author details:         Laura Aston
# Affiliation:            Public Transport Research Group, Monash University
# Contact details:        laura.aston@monash.edu
# Script and data info:   This script performs pretreatment (factor analysis) followed by multivariate multiple linear regression of built environment and sociodemographic variables on transit ridership for a sample of tram and bus locations in Melbourne. 
# Data:                   Ridership data includes average normal (school) weekday ridership 
#                         Tram ridership, averaged for 2018, by Victorian Department of Transport. 
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
#1) Radius 600m and 400m catchments were tested; 600m catchment yielded higher explanatory power solution for both modes and will be used to proceed. 
#2) Stops in the free tram zone are not included since majority of nifluential outliers are from free tram zone which has a plausible reasno for distorting values (free fare zone)

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

Melb.Trambus.600.noFTZ<- MMLR_Data[ which(MMLR_Data$Type=='Trambus'
                                     & MMLR_Data$Radius=='600'
                                     & MMLR_Data$Set_Sample_ID == 'C'
                                     & MMLR_Data$X19_FTZ == '0'),]
#Assumption 1:  'Set ID' refers to the method of defining the 'walk catchment'. 'D' (distributed) means     catchments were separately calculated for the tram and bus points, then merged. 'C' (centroid) means the geographic centroid of the points was first found, and the catchment estimated from there. Data corresponding to the 'centroid' (C) catchments is used in this analysis. 

#Universal variable list
#X2_EmpDen+X4_PopDen+X5_Dwelling+X7_PropComm+X8_RetailEmp+X9_Balance+X10_Entropy+X11_HousingDiv+X12_Intersections+X13_PBN+X14_DestScore+X15_DestCount+X16_DistCBD+X17_ACCount+X18_ACNear+X19_FTZ+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X25_Rural+X26_Access

#Step 1.1: form a data frame that comprises only the variables to be included in the factor analysis (built environment)
#include variables measured on a count or continuous scale. This means "FTZ" and "AC_Count" should be excluded. Also "rural" is 0 for almost all responses, so exclude
fa.data.Melb.Trambus.600.NoFTZ<-Melb.Trambus.600.noFTZ[,c(19, 21, 22, 24:33,35,41,43)]

#Step 1.2: Specify number of factors. Based on theory, will try four or five factors, consituting 1) density 2) diversity 3) design 4) regional accessibility and 5) local accessibility/walkability (Ewing and Cevero 2010, Voulgaris et al. 2017)

#alternatively, check scree plot
install.packages("nFactors")
library(nFactors)

ev_trambus_600.noFTZ <- eigen(cor(fa.data.Melb.Trambus.600.NoFTZ))# get eigenvalues
ap_trambus_600.noFTZ <- parallel(subject=nrow(fa.data.Melb.Trambus.600.NoFTZ),var=ncol(fa.data.Melb.Trambus.600.NoFTZ),
                               rep=100,cent=.05)
nS_trambus_600.noFTZ <- nScree(x=ev_trambus_600.noFTZ$values, aparallel=ap_trambus_600.noFTZ$eigen$qevpea)
plotnScree(nS_trambus_600.noFTZ) #3 eigenvalues

#Step 1.3 Run factor analysis
fa.Melb.Trambus.600.No.FTZ.5<-factanal(fa.data.Melb.Trambus.600.NoFTZ, factors = 5, rotation = "none")
fa.Melb.Trambus.600.No.FTZ.5

#remove high uniqueness variables (on 4-factor solution): PBN, ACNear

fa.data.Melb.Trambus.600.NoFTZ<-Melb.Trambus.600.noFTZ[,c(19, 21, 22, 24:29, 31:33,41,43)]

fa.Melb.Trambus.600.No.FTZ.5<-factanal(fa.data.Melb.Trambus.600.NoFTZ, factors = 5, rotation = "none")
fa.Melb.Trambus.600.No.FTZ.5

fa.Melb.Trambus.600.No.FTZ.5.promax<-factanal(fa.data.Melb.Trambus.600.NoFTZ, factors = 5, rotation = "promax")
fa.Melb.Trambus.600.No.FTZ.5.promax

fa.Melb.Trambus.600.No.FTZ.5.varimax<-factanal(fa.data.Melb.Trambus.600.NoFTZ, factors = 5, rotation = "varimax")
fa.Melb.Trambus.600.No.FTZ.5.varimax

#doesn't make sense, try 4 factors
fa.Melb.Trambus.600.No.FTZ.4<-factanal(fa.data.Melb.Trambus.600.NoFTZ, factors = 4, rotation = "none")
fa.Melb.Trambus.600.No.FTZ.4

fa.Melb.Trambus.600.No.FTZ.4.promax<-factanal(fa.data.Melb.Trambus.600.NoFTZ, factors = 4, rotation = "promax")
fa.Melb.Trambus.600.No.FTZ.4.promax
#sound solution

fa.Melb.Trambus.600.No.FTZ.4.varimax<-factanal(fa.data.Melb.Trambus.600.NoFTZ, factors = 4, rotation = "varimax")
fa.Melb.Trambus.600.No.FTZ.4.varimax
#prop comm cross loads in the varimax solution so take the promax solution

#check scree plot again
ev_trambus_600.noFTZ <- eigen(cor(fa.data.Melb.Trambus.600.NoFTZ))# get eigenvalues
ap_trambus_600.noFTZ <- parallel(subject=nrow(fa.data.Melb.Trambus.600.NoFTZ),var=ncol(fa.data.Melb.Trambus.600.NoFTZ),
                                 rep=100,cent=.05)
nS_trambus_600.noFTZ <- nScree(x=ev_trambus_600.noFTZ$values, aparallel=ap_trambus_600.noFTZ$eigen$qevpea)
plotnScree(nS_trambus_600.noFTZ)#4 eigenvalues

capture.output(fa.Melb.Trambus.600.No.FTZ.4.promax, file= "fa.Melb.Trambus.600.No.FTZ.4.promax.txt")

fa.Melb.Trambus.600.No.FTZ.4.promax


#Step 1.5 estimate factor scores and add to the master data frame
Trambus_fs_600.noFTZ <- factor.scores(fa.data.Melb.Trambus.600.NoFTZ, fa.Melb.Trambus.600.No.FTZ.4.promax)
Trambus_fs_600.noFTZ <- Trambus_fs_600.noFTZ$scores      #get the columns of factor scores for each case
Melb.Trambus.600.noFTZ<- cbind(Melb.Trambus.600.noFTZ,Trambus_fs_600.noFTZ) #append factor scores to dataset


#Step 2 estimate covariance of the outcome variables 
