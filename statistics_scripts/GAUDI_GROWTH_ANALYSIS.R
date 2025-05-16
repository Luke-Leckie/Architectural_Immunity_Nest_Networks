rm(list=ls())
gc()
library(gam)
library(interp)
library(cowplot)
library(ggplot2)
library(glmm)
library(lme4)
library(lmerTest)
library(multcomp)
library(dplyr)
library(report)
library(car)


qq_plot <- function(model, main) {
  qqnorm(resid(model), main = main)
  qqline(resid(model))
}
THEME='light'
pres_factor=5
if (THEME=='light'){
  pointsize=1
  alpha=1
  alpha2=0.1
  labsize=7
  labsize2=2.5
  titlesize=8.5
  linalpha=0.2
  linwidth=0.8
  tick_length=0.1
  linewidth2=1
  linewidth3=0.4
  figwidth= 7.0866
  panels=4
  figheight=(figwidth*panels)/3.333
  text_col='black'
  
  col_path='#FFAE32'
  col_sham='#39568CFF'}

if (THEME=='presentation'){
  figwidth= 7.0866
  panels=4
  figheight=(figwidth*panels)/3.333
  pointsize=1*pres_factor
  alpha=1
  alpha2=0.1
  labsize=7*pres_factor
  labsize2=2.5*pres_factor
  titlesize=8.5*pres_factor
  linalpha=0.2
  linwidth=0.8*pres_factor
  tick_length=0.1
  linewidth2=1*pres_factor
  linewidth3=0.4
  
  text_col='white'
  myColors <- c( "PATHOGEN" = "orange","SHAM" = "cyan")
  col_path='orange'
  col_sham='cyan'}

DF_V<-read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/VOL_DATA.csv')
DF <- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/CURRENT_Gs/01-05-23_2/_STD_GRAPH_ANALYSIS_WIDTH_NOV_20.csv', header = TRUE)
#DF<-read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/Publication_data/architectural_immunity_nest_analysis.csv', header=TRUE)
DF<-subset(DF, DF$name!='R1C1SK_FRI')
DF_V<-subset(DF_V, DF_V$name!='R1C1SK_FRI')
DF
#DF <- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/CURRENT_Gs/01-05-23_2/_STD_GRAPH_ANALYSIS_WIDTH_DEG_HET.csv', header = TRUE)
#DF <- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/CURRENT_Gs/01-05-23_2/_STD_GRAPH_ANALYSIS_WIDTH_DEG.csv', header = TRUE)
#DF_SW <- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/CURRENT_Gs/01-05-23_2/GRAPH_ANALYSIS_SIGMA.csv', header = TRUE)
nrow(DF_V)
DF_TUE<-subset(DF_V, day=='TUE')
nrow(DF_TUE)
model_sqrt_v <- lmer(sqrt(volume/1000) ~ treatment+ (1|colony)+(1|week),
                     data = DF_TUE)

DF_144<-subset(DF2, day=='MON')
model_sqrt_v <- lmer(sqrt(volume/1000) ~ treatment+ (1|colony)+(1|week),
                     data = DF_144)
model_sqrt_e <- lmer(sqrt(efficiency+1) ~ treatment+ (1|colony)+(1|week),
                     data = DF_144)

summary(model_sqrt_e)
car::Anova(model_sqrt_e)
hist(resid(model_sqrt_e))
DF_144<-subset(DF, time==144)
nrow(DF_144)
DF_144$modularity
model <- lmer(log(modularity+1) ~ treatment+ (1|colony)+(1|week),
                     data = DF_144)
anova(model)

DF_144$junctions
DF_144$chamber_ids

shapiro.test(resid(model_sqrt_e))
hist(resid(model_sqrt_e))
car::Anova(model_sqrt_v, type=2)
car::Anova(model_sqrt_e, type=2)

summary(model_sqrt_e)
DF_TUE<-subset(DF, day=='TUE')
DF<-subset(DF, time>=24)
DF_V<-subset(DF_V, time>=24)
DF$volume<-DF_V$volume/1000
result <- DF %>%
  filter(treatment == "SHAM") %>%
  group_by(time) %>%
  summarise(
    mean_volume = mean(volume),
    sd_volume = sd(volume)
  )
model_sqrt_v <- lmer(sqrt(surface.area) ~ time*treatment+ (1|colony/subset)+(1|week),
                     data = DF_V)
car::Anova(model_sqrt_v, type=3)
summary(model_sqrt_v)
# View the result
print(result)
DF$sqrt_vol<-sqrt(DF$c)
head(DF)
DF$mean_clustering
model_c <- lmer(log(mean_clustering+1) ~ treatment + time*treatment+ (1|colony/subset)+(1|order),
                     data = DF)
car::Anova(model_c, type=3)
hist(resid(model_c))

model_sqrt_v <- lmer(sqrt_vol ~ treatment + time*treatment+ (1|colony/subset)+(1|order),
                     data = DF)
car::Anova(model_sqrt_v, type=3)

summary(model_sqrt_v)
DF$volfit<-predict(model_sqrt_v)
DF$volume<-DF$sqrt_vol

DF_TUE<-(subset(DF, day=='TUE'))
DF_TUE$log_density<-log(DF_TUE$density)
model_den <- lmer(log_density ~ treatment+ (1|colony)+(1|week), data = DF_TUE)
car::Anova(model_den, type=3)
summary(model_den)
model_tue <- lmer(log(efficiency) ~ treatment+ (1|colony)+(1|week), data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)
model_tue <- lmer(modularities ~ treatment+ (1|colony)+(1|week), data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)
model_tue <- glmer(weighted_degree_hetero ~ treatment+ (1|colony)+(1|week),family = 'Gamma', data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)
model_tue <- glmer(degree_hetero ~ treatment+ (1|colony)+(1|week),family = 'Gamma', data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)
model_mog <- lmer(log(mean_clustering+1) ~ time*treatment+ (1|colony/subset)+(1|week), data = DF_DATA)
car::Anova(model_mog,type =3)
summary(model_mog)

DF_DATA<-(subset(DF, day!='TUE'))
DF_DATA$log_density<-log(DF_DATA$density)
model_den <- lmer(log_density ~ time*treatment+ (1|colony/subset)+(1|order), data = DF_DATA)
car::Anova(model_den, type='III')
summary(model_den)
# Create breakpoint variables
DF$segment1 <- ifelse(DF$time <= 48, DF$time, 48)
DF$segment2 <- ifelse(DF$time >= 48 & DF$time <= 144, DF$time - 48, 0)

# Now fit the model using these segmented variables

DF$log_weighted_degree_hetero<-log(DF$weighted_degree_hetero)
model <- lmer(log_weighted_degree_hetero ~ segment1 * treatment + segment2 * treatment+ 
                       (1 | colony/subset) + (1 | order),
                     data = DF)
car::Anova(model, type='III')
resid<-resid(model)
plot(resid)
hist(resid)
summary(model)

NEW_CH_DATA<-DF
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
length(NEW_CH_DATA$time)
newdat$segment1 <- ifelse(newdat$time <= 48, newdat$time, 48)
newdat$segment2 <- ifelse(newdat$time >= 48, newdat$time - 48, 0)

newdat$log_weighted_degree_hetero <- predict(model,newdat,re.form=NA)
mm <- model.matrix(terms(model),newdat)


pvar1 <- diag(mm %*% tcrossprod(vcov(model),mm))
cmult <- 1.96 ## could use 1.96
dim(mm)
dim(vcov(model))
terms(model)
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$log_weighted_degree_hetero-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$log_weighted_degree_hetero+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$log_weighted_degree_hetero-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$log_weighted_degree_hetero+sqrt(pvar1)
)

newdat[c("weighted_degree_hetero","ci_lo","ci_hi","se_lo","se_hi")] <- exp(newdat[c("log_weighted_degree_hetero","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_weight_deg_het <- ggplot(data = NEW_CH_DATA, aes(x = time, y = weighted_degree_hetero)) +
  #geom_ribbon(data = filter(newdat_sham, time<=48), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = filter(newdat_pathogen, time<=48), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham, time<=48), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen, time<=48), size=linwidth, color = col_path) +
  #geom_ribbon(data = filter(newdat_sham, time>=48), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = filter(newdat_pathogen, time>=48), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham, time>=48), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen, time>=48), size=linwidth, color = col_path) +
  
   geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  #scale_y_continuous(breaks = c(-2, -1.5,-1,-0.5,0)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Weighted degree heterogeneity")))+
  annotate("text", x = 24, y = 2, label = "24h,\np=0.016", size = labsize2, fontface='bold')+
  annotate("text", x = 48, y = 2, label = "24h-48h,\np=0.059", size = labsize2)+
  annotate("text", x = 110, y = 2, label = "48h-144h, p=0.842", size = labsize2)
  # annotate("text", x = 60, y = 2, label = "24h - 48h ***", size = labsize2)
#annotate("text", x = 110, y = 2, label = "24h - 144h *", size = labsize2)
p_weight_deg_het  


###Efficiency
{
DF_DATA$log_efficiency<-log(DF_DATA$efficiency)
model_eff <- lmer(log_efficiency ~ time*treatment+ (1|colony/subset)+(1|week), data = DF_DATA)
car::Anova(model_eff, type=3)
summary(model_eff)
NEW_CH_DATA<-DF_DATA

time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
length(NEW_CH_DATA$time)

newdat$log_efficiency <- predict(model_eff,newdat,re.form=NA)
mm <- model.matrix(terms(model_eff),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_eff),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$log_efficiency-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$log_efficiency+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$log_efficiency-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$log_efficiency+sqrt(pvar1)
)
newdat[c("efficiency","ci_lo","ci_hi","se_lo","se_hi")] <- exp(newdat[c("log_efficiency","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_efficiency <- ggplot(data = NEW_CH_DATA, aes(x = time, y = efficiency)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = filter(newdat_pathogen, treatment == "PATHOGEN"), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  #geom_point(data=NEW_CH_DATA,aes(x = jitter_time), pch = 16, size = 3, color = factor(NEW_CH_DATA$treatment, levels = c("PATHOGEN", "SHAM"))) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Efficiency")))#+
  
  #annotate("text", x = 36.1, y = 0.39, label = "24h,\np=0.802", size = labsize2)+
  #annotate("text", x = 110, y = 0.39, label = "Treatment X time, p=0.032 *", size = labsize2, fontface='bold')

p_efficiency}


{
DF_DATA<-subset(DF, time>=24)
DF_DATA$log_density<-log(DF_DATA$density)
model_den <- lmer(log_density ~ time*treatment+ (1|colony/subset)+(1|order), data = DF_DATA)
car::Anova(model_den, type=3)
NEW_CH_DATA<-DF_DATA
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
length(NEW_CH_DATA$time)

newdat$log_density <- predict(model_den,newdat,re.form=NA)
mm <- model.matrix(terms(model_den),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_den),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$log_density-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$log_density+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$log_density-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$log_density+sqrt(pvar1)
)
newdat[c("density","ci_lo","ci_hi","se_lo","se_hi")] <- exp(newdat[c("log_density","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_density <- ggplot(data = NEW_CH_DATA, aes(x = time, y = density)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = filter(newdat_pathogen, treatment == "PATHOGEN"), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Density")))#+
#  annotate("text", x = 36, y = 0.139, label = "24h,\np=0.805", size = labsize2, color=text_col)+
  #annotate("text", x = 110, y = 0.139, label = "Treatment X time, p=0.024 *", size = labsize2, color=text_col, fontface='bold')
p_density
}
##diameter
model_tue <- glmer(log(max_edge_diameter_list) ~ treatment+ (1|colony)+(1|week),family = 'Gamma', data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)

model_tue <- glmer(log(max_edge_spatial_path_list) ~ treatment+ (1|colony)+(1|week),family = 'Gamma', data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)
DF<-DF_DATA
DF$log_variable<-log(DF$max_edge_spatial_path_list)
DF$variable<-DF$max_edge_spatial_path_list
model_log <- lmer(log_variable ~ time+treatment+ (1|colony/subset)+(1|week)+(1|order), data = DF)
car::Anova(model_log, type=3)
summary(model_log)
DF$fit_diam<-exp(predict(model_log))
NEW_CH_DATA<-DF
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
length(NEW_CH_DATA$time)

newdat$log_variable <- predict(model_log,newdat,re.form=NA)
mm <- model.matrix(terms(model_log),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_log),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$log_variable-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$log_variable+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$log_variable-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$log_variable+sqrt(pvar1)
)
newdat[c("variable","ci_lo","ci_hi","se_lo","se_hi")] <- exp(newdat[c("log_variable","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_diam <- ggplot(data = NEW_CH_DATA, aes(x = time, y = variable)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = filter(newdat_pathogen, treatment == "PATHOGEN"), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Weighted diameter")))#+
  #annotate("text", x = 110, y = 430, label = "Treatment, p=0.022 *", size = labsize2, color=text_col, fontface='bold')
p_diam





####mod
DF$variable<-DF$modularities
DF$log_variable<-log(DF$modularities+1)
model_log <- lmer(log_variable ~ time*treatment+ (1|colony/subset)+(1|week), data = DF)
#model_log <- lmer(variable ~ time*treatment+ (1|colony/subset)+(1|week)+(1|order), data = DF)

hist(resid(model_log))
car::Anova(model_log, type=3)
summary(model_log)
DF$fit_diam<-exp(predict(model_log))
NEW_CH_DATA<-DF
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
length(NEW_CH_DATA$time)

newdat$log_variable <- predict(model_log,newdat,re.form=NA)
mm <- model.matrix(terms(model_log),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_log),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$log_variable-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$log_variable+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$log_variable-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$log_variable+sqrt(pvar1)
)
newdat[c("variable","ci_lo","ci_hi","se_lo","se_hi")] <- exp(newdat[c("log_variable","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])-1
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_mod <- ggplot(data = NEW_CH_DATA, aes(x = time, y = variable)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = filter(newdat_pathogen, treatment == "PATHOGEN"), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Modularity")))#+
#  annotate("text", x = 110, y = 0.7, label = "Treatment X time, p=0.22", size = labsize2, color=text_col)
p_mod
#myColors <- c( "SHAM" = "#39568CFF","PATHOGEN" = "#FFAE32") #DCE319FF - light green
#colScale <- scale_colour_manual(name = "Treatment",values = myColors, labels=c("control", "pathogen"))
#olScale2<- scale_fill_manual(name = "Treatment",values = myColors, labels=c("control", "pathogen"))
# plots=list(p_vol, p_eff, p_den, p_mod, p_dim)
# plots=list(p_vol, p_eff, p_den, p_mod, p_dim, p_deghet)
#plots=list(p_eff, p_den, p_dim, p_deghet)
plots=list(p_diam,p_efficiency, p_density, p_mod)
if (THEME=='presentation'){
  plots=list(p_diam, p_density)
  
}else{
  myColors <- c( "SHAM" = "#39568CFF","PATHOGEN" = "#FFAE32")
}
colScale <- scale_colour_manual(name = "Treatment",values = myColors, labels=c("pathogen", "control"))
colScale2<- scale_fill_manual(name = "Treatment",values = myColors, labels=c("pathogen", "control"))
labelss=c("A)", "B)", "C)", "D)", "E)", "F)")
if (THEME=='light'){
  pointsize=0.5
  alpha=0.5
  alpha2=0.1
  labsize=6.5
  labsize2=2.25
  titlesize=7.5
  linalpha=0.2
  linwidth=0.8
  tick_length=0.1
  linewidth2=1
  linewidth3=0.4
  figwidth= 7.0866
  panels=1
  figheight=(figwidth*panels)/3.333
  figwidth=figwidth/2}
for (i in seq_along(plots)){
  if (THEME=="presentation"){plots[[i]] = plots[[i]] +
    theme(
      plot.background = element_rect(fill = "black"), # set overall plot background to black
      panel.background = element_rect(fill = "black"), # set the panel background to black
      axis.text.x = element_text(size = labsize, colour = "white"),
      plot.title=element_text(size=titlesize),
      axis.text.y = element_text(size = labsize, colour = "white"),
      axis.ticks = element_line(colour = "white"), # set the axis ticks to white
      axis.title.x = element_text(size = titlesize, colour = "white"),
      axis.title.y = element_text(size = titlesize, colour = "white"),
      #axis.title.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(l = 20, r = 10, t = 20, b = 35),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(),
      axis.line = element_line(size = 1, colour = "white"))+
    colScale+colScale2  }
  if (THEME=="light"){
    plots[[i]]=plots[[i]]#+labs(title=(labelss[i]))
    if (i %in% c(1, 2)) {
      plot_margin = margin(l = 1, r = 1, t = 5, b = 1)
    } else if (i %in% c(3, 4)) {
      plot_margin = margin(l = 1, r = 1, t = 1, b = 10)
    } else {
      plot_margin = margin(l = 1, r = 1, t = 1, b = 1)}
    plot_margin = margin(l = 1, r = 1, t = 1, b = 1)
    #    if (i==1){
    #   plots[[i]] = plots[[i]] +
    #     theme_bw()+
    #     theme(axis.text.x = element_text(size = labsize),
    #           #axis.text.x = element_blank(),# Change the text size of x-axis tick labels
    #           axis.text.y = element_text(size = labsize), # Change the text size of y-axis tick labels
    #           #axis.title.x = element_text(size = 18), # Change the text size of x-axis label
    #           axis.title.x = element_blank(),
    #           plot.title=element_text(size=titlesize),
    #           legend.text = element_text(size=labsize),
    #           axis.ticks = element_line(colour = "black"),
    #           axis.ticks.length = unit(tick_length, "cm"),
    #           legend.title = element_text(size=labsize),
    #           legend.key.size = unit(0.3, "cm"),
    #           legend.position = c(0.65, 0.98), # Adjust these values to fine-tune the positio
    #           axis.title.y = element_text(size = titlesize),
    #           #axis.title.y = element_blank(),
    #           plot.margin = plot_margin, 
    #           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
    #           axis.line = element_line(size = linewidth3))+
    #     colScale+colScale2  
    #   next
    # }
    plots[[i]] = plots[[i]] +
      theme_bw()+
      theme(axis.text.x = element_text(size = labsize),
            #axis.text.x = element_blank(),# Change the text size of x-axis tick labels
            axis.text.y = element_text(size = labsize), # Change the text size of y-axis tick labels
            #axis.title.x = element_text(size = 18), # Change the text size of x-axis label
            axis.title.x = element_blank(),
            plot.title=element_text(size=titlesize),
            axis.ticks = element_line(colour = "black"),
            axis.ticks.length = unit(tick_length, "cm"),
            
            axis.title.y = element_text(size = titlesize),
            #axis.title.y = element_blank(),
            legend.position = "none",
            plot.margin = plot_margin, 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
            axis.line = element_line(size = 0.5))+
      colScale+colScale2  
  }}
dir_save='/home/ll16598/Documents/GAUDI/plots'
single_figwidth=2.33
single_figheight=2
name_list=c('diam', 'eff', 'den','mod')
for (i in seq_along(plots)){
  
  dir_save2=paste(dir_save, name_list[[i]],sep='/')
  ggsave(filename = paste(dir_save2, ".png", sep=''), plot = plots[[i]], width = single_figwidth, height = single_figheight, dpi = 800, limitsize = F)}
  
  
p_vol=plots[[1]]; p_eff=plots[[2]]; p_den=plots[[3]]; p_mod=plots[[4]]#; p_dim=plots[[5]]; p6=plots[[6]]
p6 <- ggplot(data.frame(), aes(0, 0)) + 
  geom_blank() + 
  theme_void()

if (THEME=='presentation'){
  figwidth= 30.0866
  combined_plot<-cowplot::plot_grid(p_vol,p_eff, ncol = 2,nrow=1,
                                    rel_widths = c(1,1),
                                    align = "hv",
                                    axis = "l",
                                    hspace = 0,  # Adjust horizontal spacing
                                    vspace = 0   # Adjust vertical spacing
  )
  combined_plot <- combined_plot + 
    theme(plot.background = element_rect(fill = "white"))
  final_plot <- ggdraw(combined_plot)
  #   draw_label("Time since treatment (hours)", x =0.5, y = 0.016, angle = 0, size = titlesize)
  ggsave(filename = paste(dir_save, "nest_props_pres.png", sep='/'), plot = final_plot, width = figwidth, height = figheight*1.5, dpi = 1000, limitsize = F)
}
#combined_plot<-plot_grid(p_vol,p_eff,p_den, p_mod, ncol = 2,nrow=2, rel_widths = c(1, 1), align = "hv", axis = "l")
#combined_plot<-plot_grid(p_vol, p_eff,p_den, p_mod, p_dim, ncol = 2,nrow=2, rel_widths = c(1, 1), align = "hv", axis = "l")
combined_plot<-plot_grid(p_vol, p_eff,p_den, p_mod,   ncol = 2,
                         nrow = 2,
                         rel_widths = c(1,1),
                         align = "hv",
                         axis = "l",
                         hspace = 0,  # Adjust horizontal spacing
                         vspace = 0   # Adjust vertical spacing
)
combined_plot <- combined_plot + 
  theme(plot.background = element_rect(fill = "white"))
final_plot <- ggdraw(combined_plot) + 
draw_label("Time since treatment (hours)", x =0.52, y = 0.015, angle = 0, size = titlesize)
final_plot
ggsave(filename = paste(dir_save, "nest_network.png", sep='/'), plot = final_plot, width = figwidth, height = figheight*0.4, dpi = 1000, limitsize = F)






########components

DF<-subset(DF, time>=24)
DF2<-subset(DF, time>=24)
model_tue <- glmer(sqrt(edge_num) ~ treatment+ (1|colony)+(1|week),family = 'Gamma', data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)
model_tue <- glmer(sqrt(junction_num) ~ treatment+ (1|colony)+(1|week),family = 'Gamma', data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)

model_ju <- lmer(sqrt(junctions) ~ time*treatment+ (1|colony/subset)+(1|week), data = DF)
car::Anova(model_ju, type='III')
summary(model_ju)
DF$edges<-DF$edge_num
DF$log_edges <- log(DF$edge_num)
DF$sqrt_edges <- sqrt(DF$edge_num)

model_sqrt_ed <- lmer(sqrt_edges ~ treatment + time*treatment+ (1|colony/subset)+(1|week),
                      data = DF)
summary(model_sqrt_ed)
car::Anova(model_sqrt_ed, type='III')
THEME='light'

DF_144<-subset(DF, time==144)
nrow(DF_144)
model <- lmer(log_edges ~ treatment+ (1|colony)+(1|week),
              data = DF_144)
plot(model)
car::Anova(model, type=2)
summary(model)
hist(resid(model))
shapiro.test(resid(model))
qqnorm(resid(model))



NEW_CH_DATA<-DF
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
length(NEW_CH_DATA$time)

newdat$sqrt_edges <- predict(model_sqrt_ed,newdat,re.form=NA)
mm <- model.matrix(terms(model_sqrt_ed),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_sqrt_ed),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$sqrt_edges-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$sqrt_edges+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$sqrt_edges-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$sqrt_edges+sqrt(pvar1)
)
newdat[c("edges","ci_lo","ci_hi","se_lo","se_hi")] <- (newdat[c("sqrt_edges","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])^2
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_tunnels <- ggplot(data = NEW_CH_DATA, aes(x = time, y = edges)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = newdat_pathogen, aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  #geom_point(data=NEW_CH_DATA,aes(x = jitter_time), pch = 16, size = 3, color = factor(NEW_CH_DATA$treatment, levels = c("PATHOGEN", "SHAM"))) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Number of tunnels")))+
  labs(title="B)")
  #annotate("text", x = 110, y = 200, label = "24h - 144h *", size = labsize2)
  #annotate("text", x = 36.1, y = 200, label = "24h,\np=0.755 ", size = labsize2, color=text_col)+
  #annotate("text", x = 110, y = 200, label = "24h-144h, p=0.026", size = labsize2, color=text_col, fontface='bold')

p_tunnels      

DF2<- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/CURRENT_Gs/01-05-23_2/_STD_GRAPH_ANALYSIS_WIDTH_NOV_20.csv', header = TRUE)
DF_TUE<-subset(DF2, day=='TUE')
DF2<-subset(DF2, DF2$name!='R1C1SK_FRI')

model_tue <- glmer(sqrt(nest_entrance_num) ~ treatment+ (1|colony)+(1|week), data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)

DF2<-subset(DF2, time>=24)
library(lme4)
DF2$variable<-DF2$new_ne_count
DF2$sqrt_variable<-sqrt(DF2$variable+0.1)
model_sqrt <- lmer(sqrt_variable ~ time*treatment+ (1|colony/subset)+(1|week), data = DF2)
##model_log <- lmer(log(variable) ~ log(time)*treatment+ (1|colony/subset)+(1|week), data = DF2)
model1<-model_sqrt
resid<-residuals(model1)
shapiro.test(resid)
hist(resid)
DF2$fit<-predict(model_sqrt)^2
car::Anova(model1, type=3)
summary(model1)


NEW_CH_DATA<-DF2
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
newdat$sqrt_variable <- predict(model_sqrt,newdat,re.form=NA)
mm <- model.matrix(terms(model_sqrt),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_sqrt_ed),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$sqrt_variable-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$sqrt_variable+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$sqrt_variable-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$sqrt_variable+sqrt(pvar1)
)
newdat[c("variable","ci_lo","ci_hi","se_lo","se_hi")] <- (newdat[c("sqrt_variable","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])^2
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_entrances <- ggplot(data = NEW_CH_DATA, aes(x = time, y = variable)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = newdat_pathogen, aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  #geom_point(data=NEW_CH_DATA,aes(x = jitter_time), pch = 16, size = 3, color = factor(NEW_CH_DATA$treatment, levels = c("PATHOGEN", "SHAM"))) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Number of nest entrances")))+
 # annotate("text", x = 36.1, y = 12, label = "24h,\np=0.494 ", size = labsize2, color=text_col)+
#  annotate("text", x = 110, y = 12, label = "24h-144h, p=0.177", size = labsize2, color=text_col)+
  labs(title = "C)")
p_entrances


NE_VOL <- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/ANALYSIS/PHYSICAL_ANALYSIS/ENTRANCE_DIMS.csv', header = TRUE)
NE_VOL <- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/ANALYSIS/PHYSICAL_ANALYSIS/ENTRANCE_DIMS2.csv', header = TRUE)

DF_TUE<-subset(NE_VOL, day=='TUE')

model_tue <- lmer(log(nest_entrance_dims_list) ~ treatment+ (1|colony)+(1|week), data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)
NE_VOL$variable <- NE_VOL$nest_entrance_dims_list
NE_VOL$log_variable <- log(NE_VOL$variable)
NE_VOL$sqrt_variable <- sqrt(NE_VOL$variable)
NE_VOL_TUE<-subset(NE_VOL, day=='TUE')
#model_tue <- lmer(sqrt_variable ~ treatment+ (1|colony/subset)+(1|week), data = NE_VOL_TUE)
car::Anova(model_tue)
NE_VOL<-subset(NE_VOL, day!='TUE')
nrow(NE_VOL)
model <- lmer(variable ~ time*treatment+ (1|colony/subset)+(1|week), data = NE_VOL)

model_log <- lmer(log_variable ~ time*treatment+ (1|colony/subset)+(1|week), data = NE_VOL)
model_rank <- lmer(rank(nest_entrance_dims_list) ~ time*treatment+ (1|colony/subset)+(1|week), data = NE_VOL)

anova(model_log, model_rank)
model1<-model_log
resid<-residuals(model1)
shapiro.test(resid)
hist(resid)
plot(model1)
car::Anova(model1, type='III')
summary(model1)
NE_VOL$fit<-exp(predict(model_log))

NEW_CH_DATA<-NE_VOL
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
newdat$log_variable <- predict(model_log,newdat,re.form=NA)
mm <- model.matrix(terms(model_log),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_sqrt_ed),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$log_variable-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$log_variable+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$log_variable-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$log_variable+sqrt(pvar1)
)
newdat[c("variable","ci_lo","ci_hi","se_lo","se_hi")] <- exp(newdat[c("log_variable","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_entrance_dim <- ggplot(data = NEW_CH_DATA, aes(x = time, y = variable)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = newdat_pathogen, aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha2, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  #geom_point(data=NEW_CH_DATA,aes(x = jitter_time), pch = 16, size = 3, color = factor(NEW_CH_DATA$treatment, levels = c("PATHOGEN", "SHAM"))) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Nest entrance diameter (cm)")))+
  labs(title="D)")

p_entrance_dim


# 
# CH_VOL<-read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/VOL_ANALYSIS_NOV.csv', header = TRUE)
# 
# 
# CH_VOL<-subset(CH_VOL, time>=24)
# CH_VOL$variable<-CH_VOL$cham_diam_list
# CH_VOL$time<-as.numeric(CH_VOL$time)
# CH_VOL$variable<-CH_VOL$junc_vol_list
# CH_VOL$variable<-CH_VOL$cham_vol_list
# CH_VOL$variable<-CH_VOL$junc_diam_list
# 
# process_row <- function(row) {
#   # Extract numeric values from cham_diam_list
#   numbers_str <- unlist(strsplit(gsub("\\[|\\]", "", row$variable), ", "))
#   numbers <- as.numeric(numbers_str)
#   
#   # Replicate the other columns to match the length of numbers
#   replicated_row <- lapply(row[-which(names(row) == "variable")], function(column_value) {
#     rep(column_value, length(numbers))
#   })
#   
#   # Combine the replicated columns with the numbers
#   new_row <- cbind(data.frame(variable = numbers), data.frame(replicated_row))
#   return(new_row)
# }
# 
# 
# # Apply process_row to each row of CH_VOL
# list_of_dfs <- apply(CH_VOL, 1, function(row) process_row(as.data.frame(t(row))))
# # Combine all the data frames into one
# new_df <- do.call(rbind, list_of_dfs)
# new_df$time<-as.numeric(new_df$time)
# aggregated_df <- new_df %>%
#   group_by(colony, subset, time, week, treatment) %>%
#   summarise(
#     mean_variable = mean(variable, na.rm = TRUE),
#     .groups = 'drop'  # This drops the grouping structure from the result
#   )
# 
# 
# DF_TUE<-subset(new_df, day=='TUE')
# CH_VOL
# model_tue <- lmer(sqrt(variable) ~ treatment+ (1|colony/subset)+(1|week), data = DF_TUE)
# car::Anova(model_tue, type=3)
# summary(model_tue)
# 
# model_sqrt <- lmer(sqrt(variable) ~ time*treatment+ (1|colony/subset)+(1|week), data = new_df)
# summary(model_sqrt)
# model1<-model_sqrt
# car::Anova(model1, type='III')
# mean_values <- tapply(new_df$variable, new_df$time, mean)
# mean_data <- aggregate(variable ~ time, data = new_df, FUN = function(x) mean(x, na.rm = TRUE))
# print(mean_data)
# 
# 
# mean(big_list)
# length(big_list)
CH_VOL<-read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/ANALYSIS/PHYSICAL_ANALYSIS/ch_vol.csv', header = TRUE)
CH_VOL$variable <- CH_VOL$cham_vol_list
#mean_data <- aggregate(variable ~ day, data = CH_VOL, FUN = mean)
CH_VOL<-subset(CH_VOL, CH_VOL$name!='R1C1SK_FRI')

CH_VOL$log_variable <- log(CH_VOL$variable)
CH_VOL$sqrt_variable <- sqrt(CH_VOL$variable)
CH_VOL_TUE<-subset(CH_VOL, day=='TUE')
model_tue <- lmer(sqrt_variable ~ treatment+ (1|colony/subset)+(1|week), data = CH_VOL_TUE)
car::Anova(model_tue)
CH_VOL<-subset(CH_VOL, day!='TUE')

model_sqrt <- lmer(sqrt_variable ~ time*treatment+ (1|colony/subset)+(1|week), data = CH_VOL)
model1<-model_sqrt
resid<-residuals(model1)
shapiro.test(resid)
hist(resid)
plot(model1)
car::Anova(model1, type='III')
summary(model1)
CH_VOL$fit<-predict(model_sqrt)^2
nrow(CH_VOL)
NEW_CH_DATA<-CH_VOL
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
newdat$sqrt_variable <- predict(model_sqrt,newdat,re.form=NA)
mm <- model.matrix(terms(model_sqrt),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_sqrt_ed),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$sqrt_variable-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$sqrt_variable+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$sqrt_variable-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$sqrt_variable+sqrt(pvar1)
)
newdat[c("variable","ci_lo","ci_hi","se_lo","se_hi")] <- (newdat[c("sqrt_variable","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])^2
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_ch_vol <- ggplot(data = NEW_CH_DATA, aes(x = time, y = variable)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = newdat_pathogen, aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha2, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  #geom_point(data=NEW_CH_DATA,aes(x = jitter_time), pch = 16, size = 3, color = factor(NEW_CH_DATA$treatment, levels = c("PATHOGEN", "SHAM"))) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  scale_y_continuous(breaks=c(0.0, 1.0,2.0,3.0,4.0,5.0,6.0,7.0),
                     labels = sprintf("%.1f", c(0.0, 1.0,2.0,3.0,4.0,5.0,6.0,7.0))) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Chamber volume (cm"^"3",")")))#+
  #a#nnotate("text", x = 36.1, y = 5.5, label = "24h,\np=0.844 ", size = labsize2, color=text_col)+
  #annotate("text", x = 110, y = 5.5, label = "24h-144h, p=0.440", size = labsize2, color=text_col)
p_ch_vol

##Volume

DF_V<-subset(DF_V, time>=24)
DF$volume<-DF_V$volume/1000
DF$sqrt_vol<-sqrt(DF$volume)
model_sqrt <- lmer(sqrt_vol ~ treatment + time*treatment+ (1|colony/subset)+(1|week),
                     data = DF)

NEW_CH_DATA<-DF
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
newdat$sqrt_vol <- predict(model_sqrt,newdat,re.form=NA)
mm <- model.matrix(terms(model_sqrt),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_sqrt_ed),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$sqrt_vol-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$sqrt_vol+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$sqrt_vol-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$sqrt_vol+sqrt(pvar1)
)
newdat[c("volume","ci_lo","ci_hi","se_lo","se_hi")] <- (newdat[c("sqrt_vol","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])^2
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_vol <- ggplot(data = NEW_CH_DATA, aes(x = time, y = volume)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = newdat_pathogen, aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  #geom_point(data=NEW_CH_DATA,aes(x = jitter_time), pch = 16, size = 3, color = factor(NEW_CH_DATA$treatment, levels = c("PATHOGEN", "SHAM"))) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Chamber volume (cm"^"3",")")))+
  ylab(expression(paste("Nest volume (cm"^"3",")")))#+
 # annotate("text", x = 36.1, y = 25, label = "24h,\np=0.703 ", size = labsize2, color=text_col)+
#  annotate("text", x = 110, y = 25, label = "24h-144h, p=0.015", size = labsize2, color=text_col, fontface='bold')
  #annotate("text", x = 110, y = 23, label = "24h - 144h *", size = labsize2)



DF2 <- read.csv('/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/CURRENT_Gs/01-05-23_2/_STD_GRAPH_ANALYSIS_WIDTH_NOV_20.csv', header = TRUE)
DF2<-subset(DF2, DF2$name!='R1C1SK_FRI')

DF2$chamber_num<-DF2$new_cham_count



DF_TUE<-subset(DF2, day=='TUE')
model_tue <- lmer(log(chamber_num+0.1) ~ treatment+ (1|colony)+(1|week), data = DF_TUE)
car::Anova(model_tue, type=3)
summary(model_tue)
DF2<-subset(DF2, time>=24)
sum(DF2$chamber_num)
DF2$sqrt_chamber_number<-sqrt(DF2$chamber_num)
DF2$log_chamber_number<-log(DF2$chamber_num+0.1)

##144
DF_144<-subset(DF2, time==144)
nrow(DF_144)
model <- lmer(log_chamber_number ~ treatment+ (1|colony)+(1|week),
              data = DF_144)
plot(model)
car::Anova(model, type=2)
summary(model)
hist(resid(model))
shapiro.test(resid(model))
qqnorm(resid(model))
###

model_log <- lmer(log_chamber_number ~ time*treatment+ (1|colony/subset)+(1|week), data = DF2)
model_sqrt <- lmer(sqrt_chamber_number ~ time*treatment+ (1|colony/subset)+(1|week), data = DF2)
summary(model_log)
Anova(model_log,type=3)
hist(resid(model_log))
NEW_CH_DATA<-DF2
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
newdat$log_chamber_number <- predict(model_log,newdat,re.form=NA)
mm <- model.matrix(terms(model_log),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_log),mm))
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , sqrt_ci_lo = newdat$log_chamber_number-cmult*sqrt(pvar1)
  , sqrt_ci_hi = newdat$log_chamber_number+cmult*sqrt(pvar1)
  ,sqrt_se_lo  = newdat$log_chamber_number-sqrt(pvar1)
  ,sqrt_se_hi  = newdat$log_chamber_number+sqrt(pvar1)
)
newdat[c("chamber_num","ci_lo","ci_hi","se_lo","se_hi")] <- exp(newdat[c("log_chamber_number","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")])-0.1
newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))

p_chamber_num <- ggplot(data = NEW_CH_DATA, aes(x = time, y = chamber_num)) +
  #geom_ribbon(data = filter(newdat_sham), aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_sham, linalpha), color = NA) +
  #geom_ribbon(data = newdat_pathogen, aes(ymin = se_lo, ymax = se_hi), fill = alpha(col_path, linalpha), size=linwidth, color = NA) +
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  #geom_point(data=NEW_CH_DATA,aes(x = jitter_time), pch = 16, size = 3, color = factor(NEW_CH_DATA$treatment, levels = c("PATHOGEN", "SHAM"))) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Number of chambers")))#+
 # annotate("text", x = 36.1, y = 16, label = "24h,\np=0.308 ", size = labsize2, color=text_col)+
#  annotate("text", x = 110, y = 16, label = "24h-144h, p=0.182", size = labsize2, color=text_col)
car::Anova(model_log, type=3)
summary(model_log)
p_chamber_num
myColors <- c( "SHAM" = "#39568CFF", "PATHOGEN" = "#FFAE32") #DCE319FF - light green
colScale <- scale_colour_manual(name = "Treatment",values = myColors, labels=c( "control", "pathogen"))
colScale2<- scale_fill_manual(name = "Treatment",values = myColors, labels=c( "control", "pathogen"))
# plots=list(p_vol, p_eff, p_den, p_mod, p_dim, p_deghet)
#plots=list(p_eff, p_den, p_dim, p_deghet)
plots=list(p_vol, p_tunnels, p_chamber_num, p_ch_vol, p_entrances)
if (THEME=='presentation'){
  plots=list(p_interchamber_density, p_entrancetochamber)}else{
    myColors <- c( "SHAM" = "#39568CFF","PATHOGEN" = "#FFAE32")}
myColors <- c( "SHAM" = "#39568CFF","PATHOGEN" = "#FFAE32")
colScale <- scale_colour_manual(name = "Treatment",values = myColors, labels=c("pathogen", "control"))
colScale2<- scale_fill_manual(name = "Treatment",values = myColors, labels=c("pathogen", "control"))
shapescale<-scale_shape_manual(name = "Treatment",values=c("SHAM"=10, "PATHOGEN"=21), labels=c("pathogen", "control"))

labelss=c("A)", "B)", "C)", "D)", "E)", "F)")
for (i in seq_along(plots)){
  if (THEME=="presentation"){plots[[i]] = plots[[i]] +
    theme(
      plot.background = element_rect(fill = "black"), # set overall plot background to black
      panel.background = element_rect(fill = "black"), # set the panel background to black
      axis.text.x = element_text(size = labsize, colour = "white"),
      plot.title=element_text(size=titlesize),
      axis.text.y = element_text(size = labsize, colour = "white"),
      axis.ticks = element_line(colour = "white"), # set the axis ticks to white
      axis.title.x = element_text(size = titlesize, colour = "white"),
      axis.title.y = element_text(size = titlesize, colour = "white"),
      #axis.title.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(l = 10, r = 10, t = 20, b = 10),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(),
      axis.line = element_line(size = 1, colour = "white"))+
    colScale+colScale2  }
  if (THEME=="light"){
    plots[[i]]=plots[[i]]+labs(title=(labelss[i]))
    if (i %in% c(1, 2)) {
      plot_margin = margin(l = 1, r = 1, t = 5, b = 1)
    } else if (i %in% c(5, 6)) {
      plot_margin = margin(l = 1, r = 1, t = 1, b = 10)
    } else {
      plot_margin = margin(l = 1, r = 1, t = 1, b = 1)
    }
    
    # if (i==1){
    #   plots[[i]] = plots[[i]] +
    #     theme_bw()+
    #     theme(axis.text.x = element_text(size = labsize),
    #           #axis.text.x = element_blank(),# Change the text size of x-axis tick labels
    #           axis.text.y = element_text(size = labsize), # Change the text size of y-axis tick labels
    #           #axis.title.x = element_text(size = 18), # Change the text size of x-axis label
    #           axis.title.x = element_blank(),
    #           plot.title=element_text(size=titlesize),
    #           legend.text = element_text(size=labsize),
    #           axis.ticks = element_line(colour = "black"),
    #           axis.ticks.length = unit(tick_length, "cm"),
    #           legend.title = element_text(size=labsize),
    #           legend.key.size = unit(0.3, "cm"),
    #           legend.position = c(0.65, 0.3), # Adjust these values to fine-tune the positio
    #           axis.title.y = element_text(size = titlesize),
    #           #axis.title.y = element_blank(),
    #           legend.background = element_rect(fill = NA, colour = NA), # Make legend box transparent
    #           plot.margin = margin(l = 5, r = 5, t = 5, b = 5),
    #           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
    #           axis.line = element_line(size = linewidth3))+
    #     colScale+colScale2  + shapescale
    #   next
    # }
    plots[[i]] = plots[[i]] +
      theme_bw()+
      theme(axis.text.x = element_text(size = labsize),
            #axis.text.x = element_blank(),# Change the text size of x-axis tick labels
            axis.text.y = element_text(size = labsize), # Change the text size of y-axis tick labels
            #axis.title.x = element_text(size = 18), # Change the text size of x-axis label
            axis.title.x = element_blank(),
            plot.title=element_text(size=titlesize),
            axis.ticks = element_line(colour = "black"),
            axis.ticks.length = unit(tick_length, "cm"),
            plot.margin = plot_margin, 
            axis.title.y = element_text(size = titlesize),
            #axis.title.y = element_blank(),
            legend.position = "none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
            axis.line = element_line(size = linewidth3))+
      colScale+colScale2+ shapescale
  }}
  dir_save='/home/ll16598/Documents/ARCHITECTURAL_IMMUNITY/PLOTS/FINAL'
  pb <- ggplot(data.frame(), aes(0, 0)) + 
    geom_blank() + 
    theme_void()
  figwidth= 7.0866
  panels=4
  figheight=(figwidth*panels)/3.333
  p0=plots[[1]]; p1=plots[[2]]; p2=plots[[3]]; p3=plots[[4]]; p4=plots[[5]]; p_6=plots[[6]]
  combined_plot<-cowplot::plot_grid(p0,p1,p2,p3,p4,pb, ncol = 2,nrow=3,
                                    rel_widths = c(1,1),
                                    align = "hv",
                                    axis = "l",
                                    hspace = 0,  # Adjust horizontal spacing
                                    vspace = 0   # Adjust vertical spacing
  )
  combined_plot <- combined_plot + 
    theme(plot.background = element_rect(fill = "white"))
  final_plot <- ggdraw(combined_plot) + 
    draw_label("Time since treatment (hours)", x =0.5, y = 0.016, angle = 0, size = titlesize)
  final_plot
  ggsave(filename = paste(dir_save, "nest_components_jul.png", sep='/'), plot = final_plot, width = figwidth, height = figheight*0.6, dpi = 900, limitsize = F)
  
