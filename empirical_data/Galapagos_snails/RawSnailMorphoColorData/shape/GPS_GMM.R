# Written by Andrew Kraemer - last edited 12 April 2019
# modified by A Rominger on 24 April 2019 


library(geomorph)
library(abind)

gps_gmm_owd <- setwd('~/Dropbox/Research/continuousity/MESS/empirical_data/galapagos_snails/RawSnailMorphoColorData/shape')

# Read in data (shells were digitized twice, and the code here will get shape variables and an average for each shell before analyses)
scansA<-readland.tps("scansA2019.tps",specID="ID",readcurves=T)
scansB<-readland.tps("scansB2019.tps",specID="ID",readcurves=T)
AandB<-abind(scansA,scansB)
IDs<-read.csv('names.csv')
curves<-as.matrix(read.csv('curveslide.csv'))

# calculate Procrustes 
GPS.gpa<-gpagen(AandB,curves=curves)
A.gpa<-GPS.gpa$coords[,,1:1789]
B.gpa<-GPS.gpa$coords[,,1790:3578]
mean.gpa<-(A.gpa+B.gpa)/2
mean.Csize<-(GPS.gpa$Csize[1:1789]+GPS.gpa$Csize[1790:3578])/2

# project into tangent space
snail.tangent.scores<-plotTangentSpace(mean.gpa,axis1=1,axis2=2)$pc.scores
size<-mean.Csize

# Extract Groups
species<-IDs$species
island<-IDs$island

snail.data.frame<-geomorph.data.frame(shape=mean.gpa,species=species,island=island)
procD.lm(shape~species,data= snail.data.frame,iter=999,RRPP=F)

#mean tangent scores for each species

modelspshape<-lm(snail.tangent.scores~species)
summary(manova(modelspshape))
yhat.full<-predict(manova(modelspshape))
LSmeans.spshape<-NULL
  for (i in 1:ncol(yhat.full)){
    temp<-tapply(yhat.full[,i],species,mean)
    LSmeans.spshape<-cbind(LSmeans.spshape,temp)
  }
SD.spshape<-NULL
  for (i in 1:ncol(snail.tangent.scores)){
    temp<-tapply(snail.tangent.scores[,i],species,sd)
    SD.spshape<-cbind(SD.spshape,temp)
  }
  
modelspsize<-lm(size~species)
yhat.full<-predict(aov(modelspsize))
LSmeans.spCsize<-tapply(yhat.full,species,mean)
SD.spCsize<-tapply(size,species,sd)
write.csv(LSmeans.spCsize,'LSmeans.spCsize.csv')
write.csv(SD.spCsize,'SD.spCsize.csv')
write.csv(LSmeans.spshape,'LSmeans.spShape.csv')
shape.var<-morphol.disparity(snail.tangent.scores~species,groups=~species)
write.csv(shape.var$Procrustes.var,'shapevar.csv')


plot(LSmeans.spshape[,1],LSmeans.spshape[,2],asp=1,pch=21,bg='magenta')
text(LSmeans.spshape[,1],LSmeans.spshape[,2],labels= rownames(LSmeans.spshape),cex = 0.5,pos=2)

setwd(gps_gmm_owd)
