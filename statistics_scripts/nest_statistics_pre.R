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
data_dir<-'/home/ll16598/Documents/GAUDI/Data/publication/'

## Test of global nest properties ------------------
DF <- read.csv(paste0(data_dir,'architectural_immunity_nest_analysis.csv'), header = TRUE)
DF0<-subset(DF, time==0)
DF<-subset(DF, DF$name!='R1C1SK_FRI')
DF$volume<-DF$volume/1000

vol_0<-lmer(sqrt(volume) ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(vol_0)
summary(vol_0)

sa_0<-lmer(sqrt(surface.area+1) ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(sa_0)
summary(sa_0)

j_0<-lmer(sqrt(junctions+1) ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(j_0)
summary(j_0)

tn_0<-lmer(sqrt(edge_num+1) ~ treatment+ (1|colony)+(1|week), data = DF_0)##
car::Anova(tn_0);summary(tn_0)

en_0<-lmer(sqrt(new_ne_count) ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(en_0);summary(en_0)

cn_0<-lmer(log(new_cham_count+1) ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(cn_0, type=2);summary(cn_0)

eff_0<-lmer(log(efficiency) ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(eff_0);summary(eff_0)

den_0<-lmer(log(density) ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(den_0);summary(den_0)

mod_0<-lmer(log(modularity+1) ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(mod_0);summary(mod_0)
hist(resid(mod_0))

dh_0<-glmer(degree_hetero ~ treatment+ (1|colony)+(1|week),family='Gamma', data = DF_0)
car::Anova(dh_0);summary(dh_0)

cl_0<-lmer(mean_clustering ~ treatment+ (1|colony)+(1|week), data = DF_0)
car::Anova(cl_0);summary(cl_0)

wd_0<-glmer(log(weighted_diameter) ~ treatment+ (1|colony)+(1|week),family='Gamma', data = DF_0)
car::Anova(wd_0);summary(wd_0)

uwd_0<-glmer(log(unweighted_diameter) ~ treatment+ (1|colony)+(1|week),family='Gamma', data = DF_0)
car::Anova(uwd_0);summary(uwd_0)

