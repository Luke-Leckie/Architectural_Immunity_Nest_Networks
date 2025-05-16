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

data_unmarked=read.csv(paste0(data_dir, 'unmarked_exits.csv'))
data_treated=read.csv(paste0(data_dir, 'marked_exits.csv'))

data_treated$treatment <- as.factor(data_treated$treatment)
data_treated$treatment <- relevel(data_treated$treatment, ref = "sham")

data_treated$log_marked_exits=log(data_treated$marked_exits+1)

model_marked_lm=lmer(log_marked_exits~t2*treatment+(1|colony/subset)+(1|week), data=data_treated)
car::Anova(model_marked_lm, type=3)
summary(model_marked_lm)
hist(resid(model_marked_lm))


model_unmarked_lm=lmer(log(unmarked_exits+1)~t2*treatment+(1|colony/subset)+(1|week),data=data_unmarked)
car::Anova(model_unmarked_lm, type=3)
summary(model_unmarked_lm)
