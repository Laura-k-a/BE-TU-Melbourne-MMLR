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


Melb.Trambus.600<- MMLR_Data[ which(MMLR_Data$Type=='Trambus'
                                    & MMLR_Data$Radius=='600'
                                    & MMLR_Data$Set_Sample_ID == 'C'),]

#Assumption 1:  'Set ID' refers to the method of defining the 'walk catchment'. 'D' (distributed) means     catchments were separately calculated for the tram and bus points, then merged. 'C' (centroid) means the geographic centroid of the points was first found, and the catchment estimated from there. Data corresponding to the 'centroid' (C) catchments is used in this analysis. 

#Variables to include in the factor analysis (columns 18 - 43)
#X1_Emp+X2_EmpDen+X3_Pop+X4_PopDen+X5_Dwelling+X6_ActDen+X7_PropComm+X8_RetailEmp+X9_Balance+X10_Entropy+X11_HousingDiv+X12_Intersections+X13_PBN+X14_DestScore+X15_DestCount+X16_DistCBD+X17_ACCount+X18_ACNear+X19_FTZ+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+X25_Rural+X26_Access

#Step 1.1: form a data frame that comprises only the variables to be included in the factor analysis (built environment)
#include variables measured on a count or continuous scale. This means "FTZ" and "AC_Count" should be excluded. Also "rural" is 0 for almost all responses, so exclude
fa.data.Melb.Trambus.600<-Melb.Trambus.600[,c(18:26, 28:33,35,41,43)]

#Step 1.2: Specify number of factors. Based on theory, will try four or five factors, consituting 1) density 2) diversity 3) design 4) regional accessibility and 5) local accessibility/walkability (Ewing and Cevero 2010, Voulgaris et al. 2017)

#Step 1.3 run factor analysis
#system is computationally singular as activity density is the sum of population and employment density --> do not include these in the factor matrix. 
#No stable solution found
#Could be outliers. Evaluation of outliers, based on standardized scores (z-scores) revealed 17 "high" outliers, all of which were in Melbourne's free tram fare zone, suggesting systematic bias.
#Assumption 2: Remove all stops that lie within the free tram zone. 
# Could be singularity - remove population density and employment density, the sum of which is equal to activity density. This is why 'population' and 'employment' in absolute values have also been included in the dataset

Melb.Trambus.600.noFTZ<- Melb.Trambus.600[ which(Melb.Trambus.600$X19_FTZ=='0'),]

#Step 1.3 run factor analysis
fa.data.Melb.Trambus.600.noFTZ<-Melb.Trambus.600.noFTZ[,c(18, 20, 22:26, 28:33,35,41,43)]
fa.data.Melb.Trambus.600.noFTZ<-as.matrix(fa.data.Melb.Trambus.600.noFTZ)
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 5, rotation = "none")
#unable to optimize. Try 4-factor solution
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 4, rotation = "none")
#unable to optimize. Try 3-factor solution
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "none")
fa.Melb.Trambus.600.noFTZ

#Step 1.4 Remove variables with high uniqueness values: proportion commerical, balance, PBN, ACNear, urban
fa.data.Melb.Trambus.600.noFTZ<-Melb.Trambus.600.noFTZ[,c(18, 20, 22:23, 25, 28:29, 31:33,43)]
fa.data.Melb.Trambus.600.noFTZ<-as.matrix(fa.data.Melb.Trambus.600.noFTZ)
fa.Melb.Trambus.600.noFTZ<-factanal(fa.data.Melb.Trambus.600.noFTZ, factors = 3, rotation = "none")
fa.Melb.Trambus.600.noFTZ

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

#Step 2. Estimate the covariance between dependent variables
#covariance of the logarithm
cov_Trambus_ln<-cor.test(Melb.Trambus.600.noFTZ$ln_Tram, Melb.Trambus.600.noFTZ$ln_Bus, method = "pearson", conf.level = 0.95)
cov_Trambus_ln

#covariance of the linear ridership
cov_Trambus<-cor.test(Melb.Trambus.600.noFTZ$Patronage_Tram, Melb.Trambus.600.noFTZ$Updated_Patronage_Bus, method = "pearson", conf.level = 0.95)
cov_Trambus

capture.output(cov_trambus,file="cov_trambus.csv")

#Step 3. Check for colinearity among variables
#600m
Melb.Trambus.600.noFTZ.VIF<-vif(lm(ln_Tram ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data =Melb.Trambus.600.noFTZ))
Melb.Trambus.600.noFTZ.VIF
#prop FTE has VIF = 5.1, likely colinear with Medinc vif = 4.4. Remove PropFTE

Melb.Trambus.600.noFTZ.VIF<-vif(lm(ln_Tram ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data =Melb.Trambus.600.noFTZ))
Melb.Trambus.600.noFTZ.VIF

#Step 4. Get correlation matrix and check for low loading on variables
#600m
Corrdata.Trambus.600<-Melb.Trambus.600.noFTZ[,c(45, 46,24, 26, 27, 30, 34, 35, 37, 39:41, 48:50)]
#Option 1 for Correlation matrices with p-values
Corrdata.Trambus.600<-rcorr(as.matrix(Corrdata.Trambus.600))

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

FlatCor.Trambus.600<-flattenCorrMatrix(Corrdata.Trambus.600$r,Corrdata.Trambus.600$P)
capture.output(FlatCor.Trambus.600,file="FlatCor.Trambus.600.csv")

#not significant for ln_bus
#X10_Entropy
#X13_PBN
#X23_MeanSize
#Fac_3

#ln_Tram
#X9_Balance
#X22_MedInc

#retain all

#step 4 Maximally adjusted models. 
#600m
Melb.Trambus.600.MMLR1 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.MMLR1)
Anova(Melb.Trambus.600.MMLR1)

#remove ACCount
Melb.Trambus.600.MMLR2 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.MMLR2)
Anova(Melb.Trambus.600.MMLR2)

#remove MedInc
Melb.Trambus.600.MMLR3 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.MMLR3)
Anova(Melb.Trambus.600.MMLR3)

#remove MeanSize
Melb.Trambus.600.MMLR4 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.MMLR4)
Anova(Melb.Trambus.600.MMLR4)

#remove ACNEar
Melb.Trambus.600.MMLR5 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.MMLR5)
Anova(Melb.Trambus.600.MMLR5)

capture.output(summary(Melb.Trambus.600.MMLR5), file = "Melb.Trambus.600.MMLR5.csv")

#run diagnostic
par(mfrow=c(2,2))

#600m
plot(lm(ln_Tram ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ))
#influential outliers: 286-400-C; 454-400-6; 295-400-C
which(rownames(Melb.Trambus.600.noFTZ) == "286-600-C") #38
which(rownames(Melb.Trambus.600.noFTZ) == "454-600-C") #339
which(rownames(Melb.Trambus.600.noFTZ) == "295-600-C") #120
#remove influential outlier
Melb.Trambus.600.noFTZ.rd2 <- Melb.Trambus.600.noFTZ[-c(38, 339, 120),]

#Round 2, maximally adjusted models
Melb.Trambus.600.MMLR2.1 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.MMLR2.1)
Anova(Melb.Trambus.600.MMLR2.1)

#removeACCount
Melb.Trambus.600.MMLR2.2 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.MMLR2.2)
Anova(Melb.Trambus.600.MMLR2.2)

#removeMedInc
Melb.Trambus.600.MMLR2.3 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.MMLR2.3)
Anova(Melb.Trambus.600.MMLR2.3)

#removeMean_Size
Melb.Trambus.600.MMLR2.4 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.MMLR2.4)
Anova(Melb.Trambus.600.MMLR2.4)

#removeAcNear
Melb.Trambus.600.MMLR2.5 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.MMLR2.5)
Anova(Melb.Trambus.600.MMLR2.5)

#run diagnostics
plot(lm(ln_Tram ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd2))

plot(lm(ln_Bus ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd2))
#mostoutlying observation: 197-600-C
which(rownames(Melb.Trambus.600.noFTZ.rd2) == "197-600-C") #135

#remove influential outlier
Melb.Trambus.600.noFTZ.rd3 <- Melb.Trambus.600.noFTZ.rd2[-c(135),]

#maximally adjusted model
Melb.Trambus.600.MMLR3.1 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.MMLR3.1)
Anova(Melb.Trambus.600.MMLR3.1)

#remove ACCount
Melb.Trambus.600.MMLR3.2 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.MMLR3.2)
Anova(Melb.Trambus.600.MMLR3.2)

#RemoveMedInc
Melb.Trambus.600.MMLR3.3 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.MMLR3.3)
Anova(Melb.Trambus.600.MMLR3.3)

#removeMeanSize
Melb.Trambus.600.MMLR3.4 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.MMLR3.4)
Anova(Melb.Trambus.600.MMLR3.4)

#removeAcNear
Melb.Trambus.600.MMLR3.5 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.MMLR3.5)
Anova(Melb.Trambus.600.MMLR3.5)

#run diagnostics
plot(lm(ln_Tram ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd3))

plot(lm(ln_Bus ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd3))

#not outside cook's distance, but see what happens if remove: 525-600-C, 325-600-C and 345-600-C
#mostoutlying observation: 197-600-C
which(rownames(Melb.Trambus.600.noFTZ.rd3) == "525-600-C") #175
which(rownames(Melb.Trambus.600.noFTZ.rd3) == "325-600-C") #5
which(rownames(Melb.Trambus.600.noFTZ.rd3) == "345-600-C") #2
#remove influential outlier
Melb.Trambus.600.noFTZ.rd4 <- Melb.Trambus.600.noFTZ.rd3[-c(175, 5, 2),]

#maximally adjusted model
Melb.Trambus.600.MMLR4.1 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd4)
summary(Melb.Trambus.600.MMLR4.1)
Anova(Melb.Trambus.600.MMLR4.1)

#remove AcCount
Melb.Trambus.600.MMLR4.2 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X22_MedInc+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd4)
summary(Melb.Trambus.600.MMLR4.2)
Anova(Melb.Trambus.600.MMLR4.2)

#removeMedInc
Melb.Trambus.600.MMLR4.3 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X23_MeanSize+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd4)
summary(Melb.Trambus.600.MMLR4.3)
Anova(Melb.Trambus.600.MMLR4.3)

#remove meansize
Melb.Trambus.600.MMLR4.4 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X18_ACNear+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd4)
summary(Melb.Trambus.600.MMLR4.4)
Anova(Melb.Trambus.600.MMLR4.4)

#remove ACNear
Melb.Trambus.600.MMLR4.5 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd4)
summary(Melb.Trambus.600.MMLR4.5)
Anova(Melb.Trambus.600.MMLR4.5)

#run diagnostics
plot(lm(ln_Tram ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd4))

plot(lm(ln_Bus ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd4))

#not outside cook's distance, but see what happens if remove: 413-600-C, 352-600-C and 343-600-C, 328-600-C,
#mostoutlying observation: 197-600-C
which(rownames(Melb.Trambus.600.noFTZ.rd4) == "413-600-C") #307
which(rownames(Melb.Trambus.600.noFTZ.rd4) == "352-600-C") #5
which(rownames(Melb.Trambus.600.noFTZ.rd4) == "343-600-C") #167
which(rownames(Melb.Trambus.600.noFTZ.rd4) == "328-600-C") #7

Melb.Trambus.600.noFTZ.rd5 <- Melb.Trambus.600.noFTZ.rd4[-c(307, 5, 167,7),]

Melb.Trambus.600.MMLR.rd5 <- lm(cbind(ln_Tram, ln_Bus) ~X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd5)
summary(Melb.Trambus.600.MMLR.rd5)
Anova(Melb.Trambus.600.MMLR.rd5)

#run diagnostics
plot(lm(ln_Tram ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd5))

plot(lm(ln_Bus ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd5))
#leave as is
#get the covariance matrix of the remaining observations
#covariance of the logarithm
cov_Trambus_ln.rd5<-cor.test(Melb.Trambus.600.noFTZ.rd5$ln_Tram, Melb.Trambus.600.noFTZ.rd5$ln_Bus, method = "pearson", conf.level = 0.95)
cov_Trambus_ln.rd5
#Higher; at 0.36 --> this probably does warrant keeping the MMLR

#the direction of factors 2 and 3 changes depending on the catchment buffer being 400 or 600m. 
#slightly higher explanatory power for both modes with the 600m catchment solution, so accept this. 

capture.output(summary(Melb.Trambus.600.MMLR.rd5), file = "Melb.Trambus.600.MMLR.rd5.csv")

Melb.Trambus.600.MMLR.tram<-lm(ln_Tram ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd5)

Melb.Trambus.600.MMLR.tram.beta<-lm.beta(Melb.Trambus.600.MMLR.tram)
summary(Melb.Trambus.600.MMLR.tram.beta)
capture.output(summary(Melb.Trambus.600.MMLR.tram.beta),file="Melb.Trambus.600.MMLR.tram.beta.csv")
Melb.Trambus.600.MMLR.bus<-lm(ln_Bus ~ X7_PropComm+X9_Balance+X10_Entropy+X13_PBN+X20_LOS+X24_Urban+Fac_1+Fac_2+Fac_3, data = Melb.Trambus.600.noFTZ.rd5)

Melb.Trambus.600.MMLR.bus.beta<-lm.beta(Melb.Trambus.600.MMLR.bus)
summary(Melb.Trambus.600.MMLR.bus.beta)
capture.output(summary(Melb.Trambus.600.MMLR.bus.beta),file="Melb.Trambus.600.MMLR.bus.beta.csv")
