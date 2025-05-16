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
data_dir<-''

## Test of global nest properties ------------------
DF <- read.csv(paste0(data_dir,'architectural_immunity_nest_analysis.csv'), header = TRUE)
DF0<-subset(DF, time==0)
DF<-subset(DF, DF$name!='R1C1SK_FRI')
DF$volume<-DF$volume/1000

DF_DATA<-subset(DF, time>=24)

## mixed-effects models with time × treatment interaction ------------------
vol_dynamics  <- lmer(sqrt(volume+1)                   ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
sa_dynamics   <- lmer(surface.area    ~ treatment * log(time) + (1|colony/subset) + (1|week), data = DF_DATA)
j_dynamics    <- lmer(sqrt(junctions+1)                 ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
tn_dynamics   <- lmer(sqrt(edge_num + 1)        ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
en_dynamics   <- lmer(sqrt(new_ne_count + 0.1)  ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
cn_dynamics   <- lmer(log(new_cham_count + 0.1) ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
eff_dynamics  <- lmer(log(efficiency+1)                ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
den_dynamics  <- lmer(log(density+1)                   ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
mod_dynamics  <- lmer(log(modularity+1)                ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
dh_dynamics   <- lmer(log(degree_hetero + 1)    ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
cl_dynamics   <- lmer(log(mean_clustering + 1) ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)
wd_dynamics   <- lmer(log(weighted_diameter + 1)   ~ treatment + time + (1|colony/subset) + (1|week), data = DF_DATA)
uwd_dynamics  <- lmer(unweighted_diameter ~ treatment * time + (1|colony/subset) + (1|week), data = DF_DATA)

car::Anova(vol_dynamics,  type = 3);  summary(vol_dynamics)
car::Anova(sa_dynamics,   type = 3);  summary(sa_dynamics)
car::Anova(j_dynamics,    type = 3);  summary(j_dynamics)
car::Anova(tn_dynamics,   type = 3);  summary(tn_dynamics)
car::Anova(en_dynamics,   type = 3);  summary(en_dynamics)
car::Anova(cn_dynamics,   type = 3);  summary(cn_dynamics)
car::Anova(eff_dynamics,  type = 3);  summary(eff_dynamics)
car::Anova(den_dynamics,  type = 3);  summary(den_dynamics)
car::Anova(mod_dynamics,  type = 3);  summary(mod_dynamics)
car::Anova(dh_dynamics,   type = 3);  summary(dh_dynamics)
car::Anova(cl_dynamics,   type = 3);  summary(cl_dynamics)
car::Anova(wd_dynamics,   type = 3);  summary(wd_dynamics)
car::Anova(uwd_dynamics,  type = 3);  summary(uwd_dynamics)


## mixed-effects models for test at 6-days ------------------
DF_144<-subset(DF, time==144)
vol_144<-lmer(volume ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(vol_144, type=2)
summary(vol_144)

sa_144<-lmer(sqrt(surface.area+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(sa_144, type=2); summary(sa_144)

j_144<-lmer(junctions ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(j_144, type=2);summary(j_144)

tn_144<-lmer(sqrt(edge_num+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(tn_144, type=2);summary(tn_144)

en_144<-lmer(sqrt(new_ne_count+0.1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(en_144, type=2);summary(en_144)

cn_144<-lmer(log(new_cham_count+0.1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(cn_144, type=2);summary(cn_144)

eff_144<-lmer(efficiency ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(eff_144, type=2);summary(eff_144)

den_144<-lmer(density ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(den_144, type=2);summary(den_144)

mod_144<-lmer(modularity ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(mod_144, type=2);summary(mod_144)
hist(resid(mod_144))

dh_144<-lmer(log(degree_hetero+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(dh_144, type=2);summary(dh_144)

cl_144<-lmer(sqrt(mean_clustering+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(cl_144, type=2);summary(cl_144)

wd_144<-lmer(sqrt(weighted_diameter+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(wd_144, type=2);summary(wd_144)

uwd_144<-lmer(sqrt(unweighted_diameter+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(uwd_144, type=2);summary(uwd_144)


### Chamber analysis --------------
##inter-chamber
chamber_dat <- read.csv(paste0(data_dir,'inter_chamber.csv'), header = TRUE)
chamber_dat<-subset(chamber_dat, name!='R1C1SK_FRI')
chamber_dat$variable <- chamber_dat$chamber_chamber_degs
chamber_dat$log_variable <- log(chamber_dat$variable+1)
chamber_dat$sqrt_variable <- sqrt(chamber_dat$variable+1)
DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
model <- lmer(log_variable ~ treatment+ (1|colony/week),
              data = DF_144)
car::Anova(model, type=2);summary(model)

chamber_dat_POST<-subset(chamber_dat, day!='TUE')
model <-lmer(log_variable  ~ time*treatment+ (1|colony/subset)+(1|week)+(1|chamber_id), data = chamber_dat_POST)
car::Anova(model, type='III');summary(model)
  
##degree
chamber_dat <- read.csv(paste0(data_dir,'degree.csv'), header = TRUE)
chamber_dat<-subset(chamber_dat, name!='R1C1SK_FRI')
chamber_dat$variable <- chamber_dat$chamber_degree
DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
model <- glmer(variable ~ treatment+ (1|colony)+(1|week),family='poisson',
              data = DF_144)
car::Anova(model, type=2);summary(model)

chamber_dat_POST<-subset(chamber_dat, day!='TUE')
model <-glmer(variable  ~ time+treatment+ (1|colony/subset)+(1|week)+(1|chamber_id),family='poisson', data = chamber_dat_POST)
car::Anova(model, type='II');summary(model)

##
chamber_dat <- read.csv(paste0(data_dir,'betweenness.csv'), header = TRUE)
chamber_dat<-subset(chamber_dat, name!='R1C1SK_FRI')
chamber_dat$variable <- chamber_dat$cham_flow_centrality_list
chamber_dat$log_variable <- log(chamber_dat$variable+0.1)

DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
model <- lmer(log_variable ~ treatment+ (1|colony)+(1|week),
              data = DF_144)
car::Anova(model, type=2);summary(model)

chamber_dat_POST<-subset(chamber_dat, day!='TUE')
model <-lmer(log_variable  ~ time*treatment+ (1|colony/subset)+(1|week)+(1|chamber_id), data = chamber_dat_POST)
car::Anova(model, type='III');summary(model)

chamber_dat <- read.csv(paste0(data_dir,'betweenness.csv'), header = TRUE)
chamber_dat<-subset(chamber_dat, name!='R1C1SK_FRI')
chamber_dat$variable <- chamber_dat$cham_flow_centrality_list
chamber_dat$log_variable <- log(chamber_dat$variable+0.1)

DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
model <- lmer(log_variable ~ treatment+ (1|colony)+(1|week),
              data = DF_144)
car::Anova(model, type=2);summary(model)

chamber_dat_POST<-subset(chamber_dat, day!='TUE')
model <-lmer(log_variable  ~ time*treatment+ (1|colony/subset)+(1|week)+(1|chamber_id), data = chamber_dat_POST)
car::Anova(model, type='III');summary(model)

##Efficiency centrality
chamber_dat <- read.csv(paste0(data_dir,'efficiency.csv'), header = TRUE)
chamber_dat<-subset(chamber_dat, name!='R1C1SK_FRI')
chamber_dat$variable <- chamber_dat$chamber_info_centrality
chamber_dat$log_variable <- log(chamber_dat$variable+1)

DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
model <- lmer(log_variable ~ treatment+ (1|colony)+(1|week),
              data = DF_144)

car::Anova(model, type=2);summary(model)
chamber_dat_POST<-subset(chamber_dat, day!='TUE')
model <- glmer(variable+0.1 ~ time*treatment+ (1|colony/subset)+(1|week),family='Gamma'(link='log'),control = glmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 1e5)), data = chamber_dat_POST)
car::Anova(model, type='III');summary(model)

##Entrance-entrance distance
NE_NE <- read.csv(paste0(data_dir,'entrance_entrance.csv'), header = TRUE)
NE_NE$variable <- NE_NE$any_ne_any_ne
NE_NE<-subset(NE_NE, time>=24)
NE_NE<-subset(NE_NE, name!='R1C1SK_FRI')
NE_NE144<-subset(NE_NE, day=='MON')

model <- lmer(variable~ time+treatment+ (1|colony/subset), data = NE_NE)
car::Anova(model, type='III');summary(model)

model144 <- lmer(variable ~ treatment+ (1|colony/subset)+(1|week), data = NE_NE144)
resid<-residuals(model144)
shapiro.test(resid)
hist(resid)
plot(model144)
car::Anova(model144, type='II')
summary(model144)
