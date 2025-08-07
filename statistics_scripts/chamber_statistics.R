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
direc<''
### Chamber analysis --------------
##inter-chamber
data_dir<-paste0(direc,'network_and_video_data/')
chamber_dat <- read.csv(paste0(data_dir,'inter_chamber.csv'), header = TRUE)
chamber_dat<-subset(chamber_dat, name!='R1C1SK_FRI')
chamber_dat$variable <- chamber_dat$chamber_chamber_degs
chamber_dat$log_variable <- log(chamber_dat$variable+1)

chamber_dat_POST<-subset(chamber_dat, day!='TUE')
model <-lmer(log_variable  ~ time*treatment+ (1|colony/subset)+(1|week)+(1|chamber_id), data = chamber_dat_POST)
car::Anova(model, type='III');summary(model)

DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
model <- lmer(log_variable ~ treatment+ (1|colony)+(1|week),
              data = DF_144)
car::Anova(model, type=2);summary(model)

  
##degree
chamber_dat <- read.csv(paste0(data_dir,'degree.csv'), header = TRUE)
chamber_dat<-subset(chamber_dat, name!='R1C1SK_FRI')
chamber_dat$variable <- chamber_dat$chamber_degree
DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
chamber_dat_POST<-subset(chamber_dat, day!='TUE')
model <-glmer(variable  ~ time+treatment+ (1|colony/subset)+(1|week)+(1|chamber_id),family='poisson', data = chamber_dat_POST)
car::Anova(model, type='II');summary(model)
model <- glmer(variable ~ treatment+ (1|colony)+(1|week),family='poisson',
              data = DF_144)
car::Anova(model, type=2);summary(model)


##betweenness
chamber_dat <- read.csv(paste0(data_dir,'betweenness.csv'), header = TRUE)
chamber_dat<-subset(chamber_dat, name!='R1C1SK_FRI')
chamber_dat$variable <- chamber_dat$cham_flow_centrality_list
chamber_dat$log_variable <- log(chamber_dat$variable+0.1)
chamber_dat_POST<-subset(chamber_dat, day!='TUE')
model <-lmer(log_variable  ~ time*treatment+ (1|colony/subset)+(1|week)+(1|chamber_id), data = chamber_dat_POST)
car::Anova(model, type='III');summary(model)

DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
model <- lmer(log_variable ~ treatment+ (1|colony)+(1|week),
              data = DF_144)
car::Anova(model, type=2);summary(model)

##Efficiency centrality
chamber_dat <- read.csv(paste0(data_dir,'efficiency.csv'), header = TRUE)
chamber_dat <- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/ANALYSIS/CENTRALITY_ANALYSIS/CHAM_INFO.csv', header = TRUE)

chamber_dat$variable <- chamber_dat$chamber_info_centrality
chamber_dat$log_variable <- log(chamber_dat$variable+1)
chamber_dat_POST<-subset(chamber_dat, time>=24)

model <- glmer(variable+0.1 ~ time*treatment+ (1|colony/subset)+(1|week),family='Gamma'(link='log'),control = glmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 1e5)), data = chamber_dat_POST)
car::Anova(model, type='III');summary(model)

DF_144<-subset(chamber_dat, time==144)
nrow(DF_144)
model <- lmer(log_variable ~ treatment+ (1|colony)+(1|week),
              data = DF_144)

car::Anova(model, type=2);summary(model)

