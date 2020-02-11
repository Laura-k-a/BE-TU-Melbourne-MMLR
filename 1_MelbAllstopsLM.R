# Title:                  Melbourne BE-TR All stops linear Model
# Author details:         Laura Aston
# Affiliation:            Public Transport Research Group, Monash University
# Contact details:        laura.aston@monash.edu
# Script and data info:   This script performs pretreatment (factor analysis) followed by multivariate multiple linear regression of built environment and sociodemographic variables on transit ridership for a sample of tram and bus locations in Melbourne. 
# Data:                   Ridership data includes average normal (school) weekday ridership 
#                         Train, Tram ridership, averaged for 2018, by Victorian Department of Transport. 
#                         Bus ridership, averaged for 4 months from August - November 2018, provided by Chris Loader from the bus planning team at the Victorian Department of Transport
#                         Refer to [insert doi for figshare] for ontology and reference for built environment.
#                         Copyright statement: This script is the product of Laura Aston


#This is the script for the first run through of factor and cluster analysis, involving
#modes: train/tram/bus
#variables: all
#unit of analysis: 800m (train)/600 (tram)/ 400 (bus)
#outlier: 
#analysis method: Linear regression

options(scipen = 999)
#Set working directory
setwd("C:\Users\lkast1\Google Drive\PhD\2.Analysis\2. Empirical Analysis\BE-TR_Multi Country Samples\Melbourne\Melb.All.Stops\Melb.AllStops.Repo\Data")

library(dplyr)# mutate function
library(car)#VIF function
library(Hmisc)#rcorr function
library(lm.beta)


#Refer to sampling script
MMLR_Data<- read.csv("Melb.AllStops.Data.1Jan20.csv", header=TRUE, sep=",")
row.names(MMLR_Data) <- MMLR_Data[,c(2)]

BusSample.400<-MMLR_Data[which (MMLR_Data$Mode=='bus'
                                & MMLR_Data$Target.break=='400'),]
TramSample.600<-MMLR_Data[which (MMLR_Data$Mode=='tram'
                                 & MMLR_Data$Target.break=='600'),]
TrainSample.800<-MMLR_Data[which (MMLR_Data$Mode=='train'
                                  & MMLR_Data$Target.break=='800'),]


#Col headers: ln_Patronage + PropComm +	Balance +	LUEntropy	+ PedConnect+	PBN	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ FTZ	+ Parking	+ PropUrban	+ PropRural	+ EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Tram_LOS	+ O_Train_LOS	+ O_LOS	MedInc +	PropOS +	PropBach	+ censored_PropFTE	+ censored_MeanSize	+ 	ln_Emp +	ln_Pop


setwd("C:\Users\lkast1\Google Drive\PhD\2.Analysis\2. Empirical Analysis\BE-TR_Multi Country Samples\Melbourne\Melb.All.Stops\Melb.AllStops.Repo\Regression outputs\Census")

#Bus
#step 3 Check for multicolinearity
Bus.LM.VIF<-vif(lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect+	PBN	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban	+ PropRural + EmpAccess	+ C_LOS	+ O_Tram_LOS	+ O_Train_LOS	+ MedInc +	PropOS +	PropBach	+ censored_PropFTE	+ censored_MeanSize	+	ln_Pop + ln_Emp, data =BusSample.400))
#removed  overlapping (total) level of service, FTZ to get rid of singularity
Bus.LM.VIF

#step 4 Simple correlations
Corrdata.bus<-BusSample.400[,c(19:27, 29:33, 35, 36, 38:42, 49:51)]

#Option 1 for Correlation matrices with p-values
Corrdata.bus<-rcorr(as.matrix(Corrdata.bus))

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

Corrdata.bus<-flattenCorrMatrix(Corrdata.bus$r,Corrdata.bus$P)
capture.output(Corrdata.bus,file="Corrdata.bus.csv")

#not significant for ln_bus
#censored_PropFTE

#Step 4 maximally adjusted model
Melb.bus.LM.1<-lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect+	PBN	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban	+ PropRural + EmpAccess	+ C_LOS	+ O_Tram_LOS	+ O_Train_LOS	+ MedInc +	PropOS +	PropBach	+ censored_MeanSize	+	ln_Pop + ln_Emp, data =BusSample.400)
summary(Melb.bus.LM.1)

Melb.bus.LM.1<-lm.beta(Melb.bus.LM.1)
capture.output(summary(Melb.bus.LM.1), file = "bus.MA.txt")

#function for backward (starts with maximal model)
Melb.bus.LM.1.regb<-step(Melb.bus.LM.1,
                                direction = "backward",
                                trace = 0) #don't print steps
summary(Melb.bus.LM.1.regb)

Melb.bus.LM.1.regb<-lm.beta(Melb.bus.LM.1.regb)
capture.output(summary(Melb.bus.LM.1.regb), file = "bus.PM.txt")

#diagnostics
par(mfrow=c(2,2))
plot(Melb.bus.LM.1.regb)
#1498 is influential outlier. It is more pronounced an outlier than it was in the other sample

which(rownames(BusSample.400) == "1498-bus") #212
BusSample.400.rd2 <- BusSample.400[-c(212),]

#bus round 2 maximally adjusted
#Step 4 maximally adjusted model
Melb.bus.LM.2.1<-lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect+	PBN	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban	+ PropRural + EmpAccess	+ C_LOS	+ O_Tram_LOS	+ O_Train_LOS	+ MedInc +	PropOS +	PropBach	+ censored_MeanSize	+	ln_Pop + ln_Emp, data =BusSample.400.rd2)
summary(Melb.bus.LM.2.1)

Melb.bus.LM.2.1<-lm.beta(Melb.bus.LM.2.1)
capture.output(summary(Melb.bus.LM.2.1), file = "bus.MA.rd.2.txt")

#function for backward (starts with maximal model)
Melb.bus.LM.2.regb<-step(Melb.bus.LM.2.1,
                         direction = "backward",
                         trace = 0) #don't print steps
summary(Melb.bus.LM.2.regb)
Melb.bus.LM.2.regb<-lm.beta(Melb.bus.LM.2.regb)

capture.output(summary(Melb.bus.LM.2.regb), file = "bus.PM.rd2.txt")

plot(Melb.bus.LM.2.regb)
#looks like some naturally outlying values but majority of points occur in a linear region of all graphs. Investigate: 1416-bus, 1144-bus and 1131-bus

#tram
#step 3 Check for multicolinearity
Tram.LM.VIF<-vif(lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect+	PBN	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban+ FTZ+ EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Train_LOS+	PropOS +	PropBach	+ censored_PropFTE	+ censored_MeanSize	+	ln_Pop, data =TramSample.600))
#removed  overlapping (total) level of service, propRural to get rid of singularity
#removed medInc and ln_Emp, high VIF
Tram.LM.VIF

#step 4 Simple correlations
Corrdata.tram<-TramSample.600[,c(19:30, 32:34, 36, 39:42, 49, 51)]

#Option 1 for Correlation matrices with p-values
Corrdata.tram<-rcorr(as.matrix(Corrdata.tram))

options(max.print=1000000)

Corrdata.tram<-flattenCorrMatrix(Corrdata.tram$r,Corrdata.tram$P)
capture.output(Corrdata.tram,file="Corrdata.tram.csv")

#not significant for ln_tram
#Parking
#PropBach

#Step 4 maximally adjusted model
Melb.tram.LM.1<-lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect+	PBN	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ PropUrban+ FTZ+ EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Train_LOS+	PropOS	+ censored_PropFTE	+ censored_MeanSize	+	ln_Pop, data =TramSample.600)
summary(Melb.tram.LM.1)

Melb.tram.LM.1<-lm.beta(Melb.tram.LM.1)
capture.output(summary(Melb.tram.LM.1), file = "tram.MA.txt")

#function for backward (starts with maximal model)
Melb.tram.LM.1.regb<-step(Melb.tram.LM.1,
                         direction = "backward",
                         trace = 0) #don't print steps
summary(Melb.tram.LM.1.regb)

Melb.tram.LM.1.regb<-lm.beta(Melb.tram.LM.1.regb)
capture.output(summary(Melb.tram.LM.1.regb), file = "tram.PM.txt")

#diagnostics
par(mfrow=c(2,2))
plot(Melb.tram.LM.1.regb)
#Some naturally outlying values. Most seem to fit with assumptions. Most outlying are: 170, 204 and 99. 

#train
#step 3 check for multicolinearity
Train.LM.VIF<-vif(lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect+	PBN	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban+ PropRural +  EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Tram_LOS + PropOS +	PropBach	+ censored_PropFTE	+ censored_MeanSize	+	ln_Pop, data =TrainSample.800))
#removed  overlapping (total) level of service to get rid of singularity
#removed medInc and ln_Emp, high VIF
Train.LM.VIF

#step 4 Simple correlations
Corrdata.train<-TrainSample.800[,c(19:35, 39:42, 49, 51)]

#Option 1 for Correlation matrices with p-values
Corrdata.train<-rcorr(as.matrix(Corrdata.train))

options(max.print=1000000)

Corrdata.train<-flattenCorrMatrix(Corrdata.train$r,Corrdata.train$P)
capture.output(Corrdata.train,file="Corrdata.train.csv")

#not significant for ln_tram
#PBN
#PropRural
#censored_PropFTE

#Step 4 maximally adjusted model
Melb.train.LM.1<-lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban+  EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Tram_LOS + PropOS +	PropBach	+ censored_MeanSize	+	ln_Pop, data =TrainSample.800)

summary(Melb.train.LM.1)

Melb.train.LM.1<- lm.beta(Melb.train.LM.1)
capture.output(summary(Melb.train.LM.1), file = "train.MA.txt")

#function for backward (starts with maximal model)
Melb.train.LM.1.regb<-step(Melb.train.LM.1,
                          direction = "backward",
                          trace = 0) #don't print steps
summary(Melb.train.LM.1.regb)

Melb.train.LM.1.regb<- lm.beta(Melb.train.LM.1.regb)
capture.output(summary(Melb.train.LM.1.regb), file = "train.PM.txt")

#diagnostics
plot(Melb.train.LM.1.regb)
#Samllest sample - shows its non-linearity and heteroskedasticity. In particular, 597, 621 and 628 do not fit the curves at all. 4 is close to coko's distance (investigate)

which(rownames(TrainSample.800) == "597-train") #93
which(rownames(TrainSample.800) == "621-train") #214
which(rownames(TrainSample.800) == "628-train") #210
TrainSample.800.rd2 <- TrainSample.800[-c(93, 214, 210),]

#rd 2
#Step 4 maximally adjusted model
Melb.train.LM.2.1<-lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban+  EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Tram_LOS + PropOS +	PropBach	+ censored_MeanSize	+	ln_Pop, data =TrainSample.800.rd2)

summary(Melb.train.LM.2.1)

Melb.train.LM.2.1<-lm.beta(Melb.train.LM.2.1)
capture.output(summary(Melb.train.LM.2.1), file = "train.MA.rd2.txt")

#function for backward (starts with maximal model)
Melb.train.LM.2.regb<-step(Melb.train.LM.2.1,
                           direction = "backward",
                           trace = 0) #don't print steps
summary(Melb.train.LM.2.regb)

Melb.train.LM.2.regb<-lm.beta(Melb.train.LM.2.regb)
capture.output(summary(Melb.train.LM.2.regb), file = "train.PM.rd2.txt")

#diagnostics
plot(Melb.train.LM.2.regb)
#more clear misfits: 627, 609 636
which(rownames(TrainSample.800.rd2) == "627-train") #73
which(rownames(TrainSample.800.rd2) == "609-train") #72
which(rownames(TrainSample.800.rd2) == "636-train") #79
TrainSample.800.rd3 <- TrainSample.800.rd2[-c(73, 72, 79),]

#rd 3
#Step 4 maximally adjusted model
Melb.train.LM.3.1<-lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban+  EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Tram_LOS + PropOS +	PropBach	+ censored_MeanSize	+	ln_Pop, data =TrainSample.800.rd3)

summary(Melb.train.LM.3.1)

Melb.train.LM.3.1<-lm.beta(Melb.train.LM.3.1)
capture.output(summary(Melb.train.LM.3.1), file = "train.MA.rd3.txt")

#function for backward (starts with maximal model)
Melb.train.LM.3.regb<-step(Melb.train.LM.3.1,
                           direction = "backward",
                           trace = 0) #don't print steps
summary(Melb.train.LM.3.regb)

Melb.train.LM.3.regb<-lm.beta(Melb.train.LM.3.regb)
capture.output(summary(Melb.train.LM.3.regb), file = "train.PM.rd3.txt")

plot(Melb.train.LM.3.regb)
#look like last remaining outliers (since at tail of residuals and scale-location plots): 648, 642, 632

which(rownames(TrainSample.800.rd3) == "648-train") #208
which(rownames(TrainSample.800.rd3) == "642-train") #66
which(rownames(TrainSample.800.rd3) == "632-train") #72
TrainSample.800.rd4 <- TrainSample.800.rd3[-c(208, 66, 72),]

#rd 4
#Step 4 maximally adjusted model
Melb.train.LM.4.1<-lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban+  EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Tram_LOS + PropOS +	PropBach	+ censored_MeanSize	+	ln_Pop, data =TrainSample.800.rd4)

summary(Melb.train.LM.4.1)
#explanatory power went significantly down

Melb.train.LM.4.1<-lm.beta(Melb.train.LM.4.1)
capture.output(summary(Melb.train.LM.4.1), file = "train.MA.rd4.txt")

#function for backward (starts with maximal model)
Melb.train.LM.4.regb<-step(Melb.train.LM.4.1,
                           direction = "backward",
                           trace = 0) #don't print steps
summary(Melb.train.LM.4.regb)

Melb.train.LM.4.regb<-lm.beta(Melb.train.LM.4.regb)
capture.output(summary(Melb.train.LM.4.regb), file = "train.PM.rd4.txt")

plot(Melb.train.LM.4.regb)

#615 influential outlier (but otherwise assumptions are met)
which(rownames(TrainSample.800.rd4) == "615-train") #70
TrainSample.800.rd5 <- TrainSample.800.rd4[-c(70),]

#rd 5
#Step 4 maximally adjusted model
Melb.train.LM.5.1<-lm(ln_Patronage ~ PropComm +	Balance +	LUEntropy	+ PedConnect	+ DestScore	+ Parkiteer	+ ACDist	+ ACCount	+ Parking	+ PropUrban+  EmpAccess	+ C_LOS	+ O_Bus_LOS	+ O_Tram_LOS + PropOS +	PropBach	+ censored_MeanSize	+	ln_Pop, data =TrainSample.800.rd5)

summary(Melb.train.LM.5.1)
#explanatory power back up (second highest)

Melb.train.LM.5.1<-lm.beta(Melb.train.LM.5.1)
capture.output(summary(Melb.train.LM.5.1), file = "train.MA.rd5.txt")

#function for backward (starts with maximal model)
Melb.train.LM.5.regb<-step(Melb.train.LM.5.1,
                           direction = "backward",
                           trace = 0) #don't print steps
summary(Melb.train.LM.5.regb)


Melb.train.LM.5.regb<-lm.beta(Melb.train.LM.5.regb)
capture.output(summary(Melb.train.LM.5.regb), file = "train.PM.rd5.txt")

plot(Melb.train.LM.5.regb)
#looks good; investigate most outlying values still in the model: 4, 559, 467, 2
