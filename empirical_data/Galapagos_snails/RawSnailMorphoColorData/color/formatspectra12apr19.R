# Written by Andrew Kraemer - last edited 12 April 2019

#basic loads
rm(list=ls())

#read in data
snails<-read.csv("masterspectra14to18.csv")

#format color data
brightness<-snails$brightness
spectra<-snails[,-(1:11)]
spectra.std<-spectra/brightness #if you want to evaluate spectra, first take a pca of 'spectra.std' - this has 'brightness' information stripped away, leaving color