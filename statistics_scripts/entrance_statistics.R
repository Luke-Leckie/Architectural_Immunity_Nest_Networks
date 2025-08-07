rm(list=ls())
gc()
library(gam)
library(interp)
library(cowplot)
library(ggplot2)
library(glmm)
library(lme4)
library(lmerTest)
library(dplyr)
library(report)
library(rsq)
library(car)

### Chamber analysis --------------
direc<-''
data_dir<-paste0(direc,'network_and_video_data/')

##Entrance-entrance distance
NE_NE <- read.csv(paste0(data_dir,'entrance_entrance.csv'), header = TRUE)
NE_NE<-subset(NE_NE, time>=24)

NE_NE<-subset(NE_NE, name!='R1C1SK_FRI')
NE_NE144<-subset(NE_NE, day=='MON')

model <- lmer(variable~ time+treatment+ (1|colony/subset)+(1|week), data = NE_NE)
car::Anova(model, type='II');summary(model)

model144 <- lmer(variable ~ treatment+ (1|colony/subset)+(1|week), data = NE_NE144)
car::Anova(model144, type='II')
summary(model144)
