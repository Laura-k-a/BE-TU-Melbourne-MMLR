#not significant for bus
#X13_PBN
#X23_MeanSize
#Factor3
#Factor4

#Tram
#All significant

Melb.Trambus.600.noFTZ.MMLR.L1<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.noFTZ.MMLR.L1)
Anova(Melb.Trambus.600.noFTZ.MMLR.L1)

#removeACCount
Melb.Trambus.600.noFTZ.MMLR.L2<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.noFTZ.MMLR.L2)
Anova(Melb.Trambus.600.noFTZ.MMLR.L2)

#remove mean size
Melb.Trambus.600.noFTZ.MMLR.L3<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.noFTZ.MMLR.L3)
Anova(Melb.Trambus.600.noFTZ.MMLR.L3)

#remove medinc
Melb.Trambus.600.noFTZ.MMLR.L4<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.noFTZ.MMLR.L4)
Anova(Melb.Trambus.600.noFTZ.MMLR.L4)

#remove ACNEar
Melb.Trambus.600.noFTZ.MMLR.L5<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X20_LOS+X21_PropFTE+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ)
summary(Melb.Trambus.600.noFTZ.MMLR.L5)
Anova(Melb.Trambus.600.noFTZ.MMLR.L5)

#diagnostic plots
par(mfrow=c(2,2))
plot(lm(Patronage_Tram ~ X13_PBN+X20_LOS+X21_PropFTE+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ))
#influential outlier: 295-600-C, 347-600-C, 345-600-C

plot(lm(Updated_Patronage_Bus ~ X13_PBN+X20_LOS+X21_PropFTE+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ))
#influential: 344-600-C, 346-600-C, 295-600-C

which(rownames(Melb.Trambus.600.noFTZ) == "295-600-C") #120
which(rownames(Melb.Trambus.600.noFTZ) == "347-600-C") #1
which(rownames(Melb.Trambus.600.noFTZ) == "345-600-C") #2
which(rownames(Melb.Trambus.600.noFTZ) == "344-600-C") #4
which(rownames(Melb.Trambus.600.noFTZ) == "346-600-C") #3
#remove influential outlier
Melb.Trambus.600.noFTZ.rd2 <- Melb.Trambus.600.noFTZ[-c(1, 2, 3, 4, 120),]

#rerun maximally adjusted model
Melb.Trambus.600.noFTZ.MMLR.L2.1<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.noFTZ.MMLR.L2.1)
Anova(Melb.Trambus.600.noFTZ.MMLR.L2.1)

#remove AC Count
Melb.Trambus.600.noFTZ.MMLR.L2.2<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.noFTZ.MMLR.L2.2)
Anova(Melb.Trambus.600.noFTZ.MMLR.L2.2)

#remove mean size
Melb.Trambus.600.noFTZ.MMLR.L2.3<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.noFTZ.MMLR.L2.3)
Anova(Melb.Trambus.600.noFTZ.MMLR.L2.3)

#remove Factor4
Melb.Trambus.600.noFTZ.MMLR.L2.4<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3, data =Melb.Trambus.600.noFTZ.rd2)
summary(Melb.Trambus.600.noFTZ.MMLR.L2.4)
Anova(Melb.Trambus.600.noFTZ.MMLR.L2.4)

#solution makes a lot of sense although less exlanatory power
#run diagnostics

plot(lm(Patronage_Tram ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3, data=Melb.Trambus.600.noFTZ.rd2))
#residual plots suggest non-adherence to assumption of linearity --> use logarithms

plot(lm(Updated_Patronage_Bus ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3, data=Melb.Trambus.600.noFTZ.rd2))
#influential outlier: 348-600-C

which(rownames(Melb.Trambus.600.noFTZ.rd2) == "348-600-C") #4
#remove influential outlier
Melb.Trambus.600.noFTZ.rd3 <- Melb.Trambus.600.noFTZ.rd2[-c(4),]

#rerun maximally adjusted model
Melb.Trambus.600.noFTZ.MMLR.L3.1<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X17_ACCount+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.noFTZ.MMLR.L3.1) #explanatory power has decreased even further
Anova(Melb.Trambus.600.noFTZ.MMLR.L3.1)

#remove ACCount
Melb.Trambus.600.noFTZ.MMLR.L3.2<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+X23_MeanSize+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.noFTZ.MMLR.L3.2)
Anova(Melb.Trambus.600.noFTZ.MMLR.L3.2)

#remove MeanSize
Melb.Trambus.600.noFTZ.MMLR.L3.3<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3+Factor4, data =Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.noFTZ.MMLR.L3.3)
Anova(Melb.Trambus.600.noFTZ.MMLR.L3.3)

#remove factor 4
Melb.Trambus.600.noFTZ.MMLR.L3.4<-lm(cbind(Patronage_Tram, Updated_Patronage_Bus) ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3, data =Melb.Trambus.600.noFTZ.rd3)
summary(Melb.Trambus.600.noFTZ.MMLR.L3.4)
Anova(Melb.Trambus.600.noFTZ.MMLR.L3.4)

#run diagnostics

plot(lm(Patronage_Tram ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3, data =Melb.Trambus.600.noFTZ.rd3))
#residual plots suggest non-adherence to assumption of linearity --> use logarithms

plot(lm(Updated_Patronage_Bus ~ X13_PBN+X18_ACNear+X20_LOS+X21_PropFTE+X22_MedInc+Factor1+Factor2+Factor3, data =Melb.Trambus.600.noFTZ.rd3))

#as outliers are being removed, explanatory power is decreasing and residual plots are getting worse, suggesting the relationship is non-linear