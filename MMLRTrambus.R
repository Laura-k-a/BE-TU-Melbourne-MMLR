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

#Assumption 1:  'Set ID' refers to the method of defining the 'walk catchment'. 'D' (distributed) means     catchments were separately calculated for the tram and bus points, then merged. 'C' (centroid) means the geographic centroid of the points was first found, and the catchment estimated from there. Data corresponding to the 'centroid' (C) catchments is used in this analysis. 

#Variables to include in the factor analysis (columns 18 - 43)
#X1_Emp+X2_EmpDen+X3_Pop+X4_PopDen+X5_Dwelling+X6_ActDen+X7_PropComm+X8_RetailEmp+X9_Balance+X10_Entropy+X11_HousingDiv+X12_Intersections+X13_PBN+X14_DestScore+X15_DestCount+X16_DistCBD+X17_ACCount+X18_ACNear+X19_FTZ+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X25_Rural+X26_Access

#form a data frame that comprises only the variables to be included in the factor analysis (built environment)
#include variables measured on a count or continuous scale. This means "FTZ" and "AC_Count" should be excluded. Also "rural" is 0 for almost all responses, so exclude


#Step 1.1 Specify number of factors. Based on theory, will try four or five factors, consituting 1) density 2) diversity 3) design 4) regional accessibility and 5) local accessibility/walkability (Ewing and Cevero 2010, Voulgaris et al. 2017)
fa.Melb.Trambus.600<-factanal(fa.data.Melb.Trambus.600, factors = 4, rotation = "none")


#system is computationally singular --> remove #Population and employment densities, which together constitute activity density (Keep absolute values of population and employment)
fa.data.Melb.Trambus.600<-fa.data.Melb.Trambus.600[-c(2,4)]
fa.Melb.Trambus.600<-factanal(fa.data.Melb.Trambus.600, factors = 4, rotation = "none")

#Unable to otopimise with 4 factors, try 5
fa.Melb.Trambus.600<-factanal(fa.data.Melb.Trambus.600, factors = 5, rotation = "none")
fa.Melb.Trambus.600
#high uniqueness for: housing diversity, PBN (i.e. bicycle connectivity), destination score, AC near. Remove these

fa.data.Melb.Trambus.600<-fa.data.Melb.Trambus.600[-c(9, 11, 12, 15)]
fa.Melb.Trambus.600<-factanal(fa.data.Melb.Trambus.600, factors = 5, rotation = "none")
fa.Melb.Trambus.600

#loadings on factor 4 and 5 are low, so try 3-factor solution
fa.Melb.Trambus.600<-factanal(fa.data.Melb.Trambus.600, factors = 3, rotation = "none")

#unable to optimise. Try 4
fa.Melb.Trambus.600<-factanal(fa.data.Melb.Trambus.600, factors = 4, rotation = "none")
fa.Melb.Trambus.600

#solution doesn't make a lot of sense and loadings on 3 and 4 still low. Try different rotation methods. 
fa.Melb.Trambus.600.promax<-factanal(fa.data.Melb.Trambus.600, factors = 4, rotation = "promax")
fa.Melb.Trambus.600.promax
#solution is intepretable, but some residual cross-loading. See if possible to get a solution with 3 factors

fa.Melb.Trambus.600.promax<-factanal(fa.data.Melb.Trambus.600, factors = 3, rotation = "promax") 
#unable to optimise. 

#2 factors?
fa.Melb.Trambus.600.promax<-factanal(fa.data.Melb.Trambus.600, factors = 2, rotation = "promax") 
fa.Melb.Trambus.600.promax

fa.Melb.Trambus.600.promax<-factanal(fa.data.Melb.Trambus.600, factors = 2, rotation = "varimax")
fa.Melb.Trambus.600.promax

#entropy and urban have high uniqueness --> would need to remove. Also some cross-loadig. Explore alternate rotations before removing these. 

fa.Melb.Trambus.600.varimax<-factanal(fa.data.Melb.Trambus.600, factors = 4, rotation = "varimax") 
fa.Melb.Trambus.600.varimax

#proportion commercial cross-loads. Try 3-factor solution
fa.Melb.Trambus.600.varimax<-factanal(fa.data.Melb.Trambus.600, factors = 3, rotation = "varimax")
fa.Melb.Trambus.600.varimax
#cannot optimise.


#try removing prop comm
fa.data.Melb.Trambus.600<-fa.data.Melb.Trambus.600[-c(5)]
fa.Melb.Trambus.600.varimax<-factanal(fa.data.Melb.Trambus.600, factors = 4, rotation = "varimax") 
fa.Melb.Trambus.600.varimax
#this introduces even more cross loading. Try promax rotation
fa.Melb.Trambus.600.promax<-factanal(fa.data.Melb.Trambus.600, factors = 4, rotation = "promax")
fa.Melb.Trambus.600.promax
#yields a logical solution with no cross-loading. Retain as best solution with FTZ in the sample

capture.output(fa.Melb.Trambus.600.promax, file = "fa.Melb.Trambus.600.promax.4.csv")

