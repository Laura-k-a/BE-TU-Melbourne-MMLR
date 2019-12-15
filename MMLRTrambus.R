# Title:                  Melbourne Tram-bus MMLR
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
#Sensitivity testing is needed to determine the impact of methodological choices on results (Aston et al 2019). The following variants will be tested for each step:
#Trambus sample (different set of options explored for trainbus sample)
#A) With and without stops in the Free Tram Zone (factorisation only)
#B)Radius: 600m vs. 400m (factorisation, through to full solution)
#C) Factorisation method: minRes and maximum likelihood estimation (factorisation, through to full solution)


#Set working directory
setwd("C:/Users/lkast1/Google Drive/PhD/2.Analysis/2. Empirical Analysis/BE-TR_Multi Country Samples/Melbourne/MMLR/BE-TU-Melbourne-MMLR")


#Install and load packages
install.packages("lm.beta")
install.packages("dplyr")
#install.packages("car")
install.packages("Hmisc")
install.packages("psych")
library(lm.beta)
library(dplyr)
#library(car)
library(Hmisc)
library(psych)

#turn off scientific notation
options(scipen = 999)

#read in data
MMLR_Data<-read.csv(file="Updated_MMLR_Data.csv")

#Transform patronage into natural log and add as a column to the dataframe
MMLR_Data<-mutate(Trambus_MMLR_Data,
                           ln_Bus = log(Updated_Patronage_Bus),
                           ln_Tram = log(Patronage_Tram),
                           ln_Train = log(Patronage_Train),)

#optional - assign the Sample ID as the row names
row.names(MMLR_Data) <- MMLR_Data[,c(3)]

#First subset to test for trambus sample is all stops, with walk catchment radius of 600 metres
#subsetting for mode in the census of all eligible stops


Melb.Trambus.600<- MMLR_Data[ which(MMLR_Data$Type=='Trambus'
                                                  & MMLR_Data$Radius=='600'
                                                  & MMLR_Data$Set_Sample_ID == 'C'),]

#NB:  'Set ID' refers to the method of defining the 'walk catchment'. 'D' (distributed) means     catchments were separately calculated for the tram and bus points, then merged. 'C' (centroid) means the geographic centroid of the points was first found, and the catchment estimated from there. Data corresponding to the 'centroid' (C) catchments is used in this analysis. 

#Variables to include in the factor analysis (columns 18 - 43)
#X1_Emp+X2_EmpDen+X3_Pop+X4_PopDen+X5_Dwelling+X6_ActDen+X7_PropComm+X8_RetailEmp+X9_Balance+X10_Entropy+X11_HousingDiv+X12_Intersections+X13_PBN+X14_DestScore+X15_DestCount+X16_DistCBD+X17_ACCount+X18_ACNear+X19_FTZ+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X25_Rural+X26_Access

#form a data frame that comprises only the variables to be included in the factor analysis (built environment)
#include variables measured on a count or continuous scale. This means "FTZ" and "AC_Count" should be excluded. Also "rural" is 0 for almost all responses, so exclude
fa.Melb.Trambus.600<-Melb.Trambus.600[,c(18:33,35,41,43)]

#Step 1.1 Determine how many factors
Melb.Trambus.600.pca <- princomp(fa.Melb.Trambus.600)
summary(Melb.Trambus.600.pca)
plot(Melb.Trambus.600.pca)

#only one component. High correlation. How to fix? Could be a case of influential outliers. Try removing stops in the free tram fare zone.
#Could also be singularity - reove population density and employment density, the sum of which is equal to activity density. This is why 'population' and 'employment' in absolute values have also been included in the dataset


Melb.Trambus.600.noFTZ<- Melb.Trambus.600[ which(Melb.Trambus.600$X19_FTZ=='0'),]

fa.data.Melb.Trambus.600.noFTZ<-Melb.Trambus.600.noFTZ[,c(18,20,22:33,35,41,43)]
Melb.Trambus.600.pca.noFTZ <- princomp(fa.data.Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.pca.noFTZ)
plot(Melb.Trambus.600.pca.noFTZ)


#scree plot suggests there is only one principal component. Proceed to factor analysis and work backwards with unique factors removed. 

fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "none")
fa.Melb.Trambus.600.noFTZ

#eliminate high uniqueness variables (>0.6): Entropy, PBN, AC Near, Urban and rerun PCA
fa.data.Melb.Trambus.600.noFTZ<-Melb.Trambus.600.noFTZ[,c(18,20,22:26, 28:29, 31:33,43)]
Melb.Trambus.600.pca.noFTZ <- princomp(fa.data.Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.pca.noFTZ)
plot(Melb.Trambus.600.pca.noFTZ)
#still suggesting one component

#rerun factor analysis on 3 factors with high uniqueness variables removed
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "none")
fa.Melb.Trambus.600.noFTZ

#propr commercial has high uniqueness, try removing
fa.data.Melb.Trambus.600.noFTZ<-Melb.Trambus.600.noFTZ[,c(18,20,22:23, 25:26, 28:29, 31:33,43)]
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "none")
fa.Melb.Trambus.600.noFTZ

#balance has high uniqueness, try removing
fa.data.Melb.Trambus.600.noFTZ<-fa.data.Melb.Trambus.600.noFTZ[-c(6)]
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "none")
fa.Melb.Trambus.600.noFTZ

#remaingin factors seem to be appropriate for the factor solution. Rerun PCA
Melb.Trambus.600.pca.noFTZ <- princomp(fa.data.Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.pca.noFTZ)
plot(Melb.Trambus.600.pca.noFTZ)
#still suggesting only one component

#try two-factor solution

fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 2, rotation = "none")
fa.Melb.Trambus.600.noFTZ

#housing div has high uniqueness --> remove
fa.data.Melb.Trambus.600.noFTZ<-fa.data.Melb.Trambus.600.noFTZ[-c(6)]
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "none")
fa.Melb.Trambus.600.noFTZ

fa.Melb.Trambus.600.noFTZ.varimax<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "varimax")
fa.Melb.Trambus.600.noFTZ.varimax

fa.Melb.Trambus.600.noFTZ.promax<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "promax")
fa.Melb.Trambus.600.noFTZ.promax
#this solution seems to make sense. But what is the implication of the PCA only suggesting one principle component?

#try running with minres. 
