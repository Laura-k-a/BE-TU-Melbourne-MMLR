# Title:                  Melbourne All Stops (cleaned version)
# Author details:         Laura Aston
# Affiliation:            Public Transport Research Group, Monash University
# Contact details:        laura.aston@monash.edu
# Script and data info:   This script performs pretreatment (factor analysis) of built environment variables measured for transit stations in Melbourne
# Data:                   Refer to [insert doi for figshare] for ontology and reference for built environment.
# Copyright statement: This script is the product of Laura Aston

# Sample attributes
#modes: all
#variables: BE and sociodemographic
#unit of analysis: 800m
#outlier: all in

#Set working directory
setwd("C:/Users/lkast1/Google Drive/PhD/2.Analysis/2. Empirical Analysis/BE-TR_Multi Country Samples/Melbourne/Melb.All.Stops/Melb.AllStops.Repo")

install.packages("nFactors")
install.packages("REdaS")
library(nFactors)
library(psych)
library(REdaS) #package for Bartlett's sphericity test
citation("REdaS")
citation("psych")

#master variable list [cols 19:53]

#20_Emp	21_Empden	22_Pop	23_Popden	24_Dwell	25_ActDen	26_PropComm	27_RetailEmp	28_Balance	29_LUEntropy	30_HousingDiv	31_PedConnect	32_PBN	33_DestScore	34_DestCount	35_Parkiteer	36_CDBDist	37_ACDist	38_ACCount	39 _FTZ	40_Parking	41_PropUrban	42_PropRural	43_EmpAccess	44_C_LOS	45_O_Bus_LOS	46_O_Tram_LOS	47_O_Train_LOS	48_O_LOS	49_PropFTE	50_MeanSize	51_MedInc	52_PropOS	53_PropBach

#Step 1:Test the variables are adequately correlated to use factor analysis (significant Bartlett's sphericity test)
Allmodes.Melb<- read.csv("Melb.AllStops.Data.18Dec19.csv", header=TRUE, sep=",")

#subset for 800m catchment

row.names(Allmodes.Melb) <- Allmodes.Melb[,c(2)]
Allmodes.Melb.800<-Allmodes.Melb[ which(Allmodes.Melb$Target.break=='800'),]

cormat.BE.Allmodes.800<-Allmodes.Melb.800[,c(21, 23, 24, 26:34, 36:37, 41:43)]
cormat.SD.Allmodes.800<-Allmodes.Melb.800[,c(49:53)]

bartletts_BE.Allmodes.Melb.800<-bart_spher(cormat.BE.Allmodes.800, use = c("complete.obs"))
bartletts_BE.Allmodes.Melb.800

bartletts_SD.Allmodes.Melb.800<-bart_spher(cormat.SD.Allmodes.800, use = c("complete.obs"))
bartletts_SD.Allmodes.Melb.800

#scree plot
install.packages("nFactors")
library(nFactors)


ev_BE <- eigen(cor(cormat.BE.Allmodes.800))# get eigenvalues
ap_BE <- parallel(subject=nrow(cormat.BE.Allmodes.800),var=ncol(cormat.BE.Allmodes.800),
                  rep=100,cent=.05)
nS_BE<- nScree(x=ev_BE$values, aparallel=ap_BE$eigen$qevpea)
plotnScree(nS_BE) #4 eigenvalues

#BE_Unrotated factor matrix
fa.BE.800.4 <- factanal(cormat.BE.Allmodes.800, factors = 4, rotation = "none", lower =0.01) #unable to optimise
fa.BE.800.3 <- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "none", lower =0.01) #unable to optimise
fa.BE.800.5 <- factanal(cormat.BE.Allmodes.800, factors = 5, rotation = "none", lower =0.01) #unable to optimise
fa.BE.800.2 <- factanal(cormat.BE.Allmodes.800, factors = 2, rotation = "none", lower =0.01)
fa.BE.800.2
#can't optimise
#remove highest uniqueness variables (>0.8)
#X29 LU Entropy
#X32_PBN
#X37_ACDist
#X42_PropRural


#20_Emp	21_Empden	22_Pop	23_Popden	24_Dwell	25_ActDen	26_PropComm	27_RetailEmp	28_Balance	29_LUEntropy	30_HousingDiv	31_PedConnect	32_PBN	33_DestScore	34_DestCount	35_Parkiteer	36_CDBDist	37_ACDist	38_ACCount	39 _FTZ	40_Parking	41_PropUrban	42_PropRural	43_EmpAccess	44_C_LOS	45_O_Bus_LOS	46_O_Tram_LOS	47_O_Train_LOS	48_O_LOS	49_PropFTE	50_MeanSize	51_MedInc	52_PropOS	53_PropBach

cormat.BE.Allmodes.800<-Allmodes.Melb.800[,c(21, 23, 24, 26, 27, 28, 30, 31, 33, 34, 36, 41:43)]
#could also be due to outliers
#update scree plot
ev_BE <- eigen(cor(cormat.BE.Allmodes.800))# get eigenvalues
ap_BE <- parallel(subject=nrow(cormat.BE.Allmodes.800),var=ncol(cormat.BE.Allmodes.800),
                  rep=100,cent=.05)
nS_BE<- nScree(x=ev_BE$values, aparallel=ap_BE$eigen$qevpea)
plotnScree(nS_BE) #4 eigenvalues

fa.BE.800.4 <- factanal(cormat.BE.Allmodes.800, factors = 4, rotation = "none", lower =0.01)
fa.BE.800.4

#remove high uniqueness variables (>0.6)
#X24 Prop Rural
cormat.BE.Allmodes.800<-Allmodes.Melb.800[,c(21, 23, 24, 26, 27, 28, 30, 31, 33, 34, 36, 41,43)]
ev_BE <- eigen(cor(cormat.BE.Allmodes.800))# get eigenvalues
ap_BE <- parallel(subject=nrow(cormat.BE.Allmodes.800),var=ncol(cormat.BE.Allmodes.800),
                  rep=100,cent=.05)
nS_BE<- nScree(x=ev_BE$values, aparallel=ap_BE$eigen$qevpea)
plotnScree(nS_BE) #3 eigenvalues

fa.BE.800.3 <- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "none", lower =0.01)
fa.BE.800.3

fa.BE.800.4 <- factanal(cormat.BE.Allmodes.800, factors = 4, rotation = "none", lower =0.01)
fa.BE.800.4

fa.BE.800.4.promax <- factanal(cormat.BE.Allmodes.800, factors = 4, rotation = "promax", lower =0.01)
fa.BE.800.4.promax

fa.BE.800.4.varimax <- factanal(cormat.BE.Allmodes.800, factors = 4, rotation = "varimax", lower =0.01)
fa.BE.800.4.varimax
#Access crossload in varimax solution; balance cross loads in promax solution. Try 3 factors

fa.BE.800.3.promax <- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "promax", lower =0.01)
fa.BE.800.3.promax

fa.BE.800.3.varimax <- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "varimax", lower =0.01)
fa.BE.800.3.varimax
#no clear solution; remove propcomm and propurban which have high uniqueness under 3-factor solution

cormat.BE.Allmodes.800<-Allmodes.Melb.800[,c(21, 23, 24, 27, 28, 30, 31, 33, 34, 36, 43)]
ev_BE <- eigen(cor(cormat.BE.Allmodes.800))# get eigenvalues
ap_BE <- parallel(subject=nrow(cormat.BE.Allmodes.800),var=ncol(cormat.BE.Allmodes.800),
                  rep=100,cent=.05)
nS_BE<- nScree(x=ev_BE$values, aparallel=ap_BE$eigen$qevpea)
plotnScree(nS_BE) #2 eigenvalues

fa.BE.800.3 <- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "none", lower =0.01)
fa.BE.800.3

fa.BE.800.3.promax <- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "promax", lower =0.01)
fa.BE.800.3.promax

fa.BE.800.3.varimax <- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "varimax", lower =0.01)
fa.BE.800.3.varimax

#remove balance and try 2-factor solution
cormat.BE.Allmodes.800<-Allmodes.Melb.800[,c(21, 23, 24, 27, 30, 31, 33, 34, 36, 43)]
fa.BE.800.2<- factanal(cormat.BE.Allmodes.800, factors = 2, rotation = "none", lower =0.01)
fa.BE.800.2

fa.BE.800.2.promax<- factanal(cormat.BE.Allmodes.800, factors = 2, rotation = "promax", lower =0.01)
fa.BE.800.2.promax

fa.BE.800.2.varimax<- factanal(cormat.BE.Allmodes.800, factors = 2, rotation = "varimax", lower =0.01)
fa.BE.800.2.varimax
#try 3-factor solution with balance removed
fa.BE.800.3<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "none", lower =0.01)
fa.BE.800.3

fa.BE.800.3.promax<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "promax", lower =0.01)
fa.BE.800.3.promax

fa.BE.800.3.varimax<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "varimax", lower =0.01)
fa.BE.800.3.varimax

#employment access has persistent cross-loading. Remove
cormat.BE.Allmodes.800<-Allmodes.Melb.800[,c(21, 23, 24, 27, 30, 31, 33, 34, 36)]

fa.BE.800.2<- factanal(cormat.BE.Allmodes.800, factors = 2, rotation = "none", lower =0.01)
fa.BE.800.2

fa.BE.800.2.promax<- factanal(cormat.BE.Allmodes.800, factors = 2, rotation = "promax", lower =0.01)
fa.BE.800.2.promax #preferred

fa.BE.800.2.varimax<- factanal(cormat.BE.Allmodes.800, factors = 2, rotation = "varimax", lower =0.01)
fa.BE.800.2.varimax

#try 3-factors
fa.BE.800.3<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "none", lower =0.01)
fa.BE.800.3

fa.BE.800.3.promax<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "promax", lower =0.01)
fa.BE.800.3.promax #ped connect cross-loading

fa.BE.800.3.varimax<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "varimax", lower =0.01)
fa.BE.800.3.varimax

#pedconnect has persistent cross-loading, try removing
cormat.BE.Allmodes.800<-Allmodes.Melb.800[,c(21, 23, 24, 27, 30, 33, 34, 36)]
fa.BE.800.3<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "none", lower =0.01)
fa.BE.800.3

fa.BE.800.3.promax<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "promax", lower =0.01)
fa.BE.800.3.promax 
#makes sense

fa.BE.800.3.varimax<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "varimax", lower =0.01)
fa.BE.800.3.varimax
#crossloading of Dest count --> accept promax solution

#check scree
ev_BE <- eigen(cor(cormat.BE.Allmodes.800))# get eigenvalues
ap_BE <- parallel(subject=nrow(cormat.BE.Allmodes.800),var=ncol(cormat.BE.Allmodes.800),
                  rep=100,cent=.05)
nS_BE<- nScree(x=ev_BE$values, aparallel=ap_BE$eigen$qevpea)
plotnScree(nS_BE) #2 eigenvalues; levels off at 2 and again at 3. Accept 3-factor solution. 

fa.BE.800.3.promax<- factanal(cormat.BE.Allmodes.800, factors = 3, rotation = "promax", lower =0.01, scores = "regression")
capture.output(fa.BE.800.3.promax , file = "fa.BE.800.3.promax.csv")
capture.output(fa.BE.800.3, file = "fa.BE.800.3.evs.csv")
fa.BE.800.3

#Step 2 - append factor scores to dataset

Melb_All_fs <- factor.scores(cormat.BE.Allmodes.800, fa.BE.800.3.promax)
Melb_All_fs <-Melb_All_fs$scores    #get the columns of factor scores for each case
Allmodes.Melb.800<- cbind(Allmodes.Melb.800,Melb_All_fs) #append factor scores to dataset

Melb_All_fs
options(max.print=10000000)
capture.output(Allmodes.Melb.800, file = "Allmodes.Melb.800.csv")
capture.output(Melb_All_fs, file = "Melb_All_fsb.800.csv")

#factor1 = Local access
#factor2 = Employment density
#factor 3 = Residential density

##Sociodemographic factors
ev_SD <- eigen(cor(cormat.SD.Allmodes.800))# get eigenvalues
ap_SD <- parallel(subject=nrow(cormat.SD.Allmodes.800),var=ncol(cormat.SD.Allmodes.800),
rep=100,cent=.05)
nS_SD<- nScree(x=ev_SD$values, aparallel=ap_SD$eigen$qevpea)
plotnScree(nS_SD) #e = 2

fa.SD.800.2<- factanal(cormat.SD.Allmodes.800, factors = 2, rotation = "none", lower =0.01)
fa.SD.800.2

#Prop FTE, Mean Size ,Prop OS have high uniqueness. Remove
#Can't run factor analysis with fewer than 3 variables. No factors for demography but may need to choose one of median income or bachelor's degree to resolve colineaity. 