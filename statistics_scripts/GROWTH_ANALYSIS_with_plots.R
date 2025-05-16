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

pw <- function(x, knot, end){
  x1 <- ifelse(x < knot, x, knot)
  x2 <- ifelse(x < knot, 0, x - knot)
  return(cbind(x1, x2, end))
}
qq_plot <- function(model, main) {
  qqnorm(resid(model), main = main)
  qqline(resid(model))
}
THEME='light'
pres_factor=5
if (THEME=='light'){
  pointsize=1
  alpha=0.6
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

DF <- read.csv('/home/ll16598/Documents/GAUDI/Data/publication/architectural_immunity_nest_analysis.csv', header = TRUE)
DF0<-subset(DF, time==0)
DF<-subset(DF, DF$name!='R1C1SK_FRI')
DF_DATA<-subset(DF, time>=24)
DF_DATA$volume<-DF_DATA$volume/1000

#Tests at final timepoint
DF_144<-subset(DF, time==144)
DF_144$volume<-DF_144$volume/1000
vol_144<-lmer(volume ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(vol_144, type=2)
summary(vol_144)

sa_144<-lmer(sqrt(surface.area+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(sa_144, type=2)
summary(sa_144)

j_144<-lmer(junctions ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(j_144, type=2)
summary(j_144)

tn_144<-lmer(sqrt(edge_num+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(tn_144, type=2)
summary(tn_144)

en_144<-lmer(sqrt(new_ne_count+0.1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(en_144, type=2)
summary(en_144)

cn_144<-lmer(log(new_cham_count+0.1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(cn_144, type=2)
summary(cn_144)

eff_144<-lmer(efficiency ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(eff_144, type=2)
summary(eff_144)

den_144<-lmer(density ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(den_144, type=2)
summary(den_144)

mod_144<-lmer(modularity ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(mod_144, type=2)
summary(mod_144)

dh_144<-lmer(log(degree_hetero+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(dh_144, type=2)
summary(dh_144)

cl_144<-lmer(sqrt(mean_clustering+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(cl_144, type=2)
summary(cl_144)

wd_144<-lmer(sqrt(weighted_diameter+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(wd_144, type=2)
summary(wd_144)

uwd_144<-lmer(sqrt(unweighted_diameter+1) ~ treatment+ (1|colony)+(1|week), data = DF_144)
car::Anova(uwd_144, type=2)
summary(uwd_144)


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

## Type-III ANOVAs ----------------------------------------------------------
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
cmult <- 1.96 
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
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Efficiency")))
p_efficiency}

p_efficiency

{
DF_DATA$log_density<-log(DF_DATA$density)
model_den <- lmer(log_density ~ time*treatment+ (1|colony/subset), data = DF_DATA)
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
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Density")))
p_density
}


{
DF_DATA$log_variable<-log(DF_DATA$weighted_diameter)
model_log <- lmer(log_variable ~ time+treatment+ (1|colony/subset)+(1|week), data = DF_DATA)
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
cmult <- 1.96
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
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Weighted diameter")))
p_diam
}

{
  DF_DATA$variable<-DF_DATA$unweighted_diameter
  model <- lmer(variable ~ time*treatment+ (1|colony/subset)+(1|week), data = DF_DATA)
  car::Anova(model, type=3)
  summary(model)
  NEW_CH_DATA<-DF_DATA
  time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
  newdat <- expand.grid(
    time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
    ,  treatment=unique(NEW_CH_DATA$treatment)
  )
  length(NEW_CH_DATA$time)
  
  newdat$variable <- predict(model,newdat,re.form=NA)
  mm <- model.matrix(terms(model),newdat)
  pvar1 <- diag(mm %*% tcrossprod(vcov(model),mm))
  cmult <- 1.96 
  newdat <- data.frame(
    newdat
    , sqrt_ci_lo = newdat$variable-cmult*sqrt(pvar1)
    , sqrt_ci_hi = newdat$variable+cmult*sqrt(pvar1)
    ,sqrt_se_lo  = newdat$variable-sqrt(pvar1)
    ,sqrt_se_hi  = newdat$variable+sqrt(pvar1)
  )
  newdat[c("variable","ci_lo","ci_hi","se_lo","se_hi")] <- newdat[c("variable","sqrt_ci_lo","sqrt_ci_hi","sqrt_se_lo","sqrt_se_hi")]
  newdat_sham     <- newdat[which(newdat$treatment=="SHAM"),]
  newdat_pathogen <- newdat[which(newdat$treatment=="PATHOGEN"),]
  NEW_CH_DATA$treatment <- factor(NEW_CH_DATA$treatment, levels = c("SHAM", "PATHOGEN"))
  
  p_udiam <- ggplot(data = NEW_CH_DATA, aes(x = time, y = variable)) +
    geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
    geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
    geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
               aes(group=treatment, color=treatment)) +
    scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
    xlab("Time since treatment (hours)")+
    ylab(expression(paste("Unweighted diameter")))
  p_udiam
}

{#Modularity
  DF_DATA$log_variable<-log(DF_DATA$modularity+1)
model_log <- lmer(log_variable ~ time*treatment+ (1|colony/subset)+(1|week), data = DF_DATA)
hist(resid(model_log))
car::Anova(model_log, type=3)
summary(model_log)
NEW_CH_DATA<-DF_DATA
time_range <- max(NEW_CH_DATA$time)-min(NEW_CH_DATA$time)
newdat <- expand.grid(
  time=seq(min(NEW_CH_DATA$time)-0.05*time_range,max(NEW_CH_DATA$time)+0.05*time_range,length.out = 100)
  ,  treatment=unique(NEW_CH_DATA$treatment)
)
length(NEW_CH_DATA$time)

newdat$log_variable <- predict(model_log,newdat,re.form=NA)
mm <- model.matrix(terms(model_log),newdat)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_log),mm))
cmult <- 1.96 
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
  geom_line(data = filter(newdat_sham), size=linwidth, color = col_sham) +
  geom_line(data = filter(newdat_pathogen), size=linwidth, color = col_path) +
  geom_point(data=NEW_CH_DATA,alpha = alpha, size=pointsize, position = position_jitterdodge(jitter.width=4.2,dodge.width = 9), 
             aes(group=treatment, color=treatment)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 144)) +
  xlab("Time since treatment (hours)")+
  ylab(expression(paste("Modularity")))#+
p_mod
}
plots=list(p_diam,p_udiam,p_efficiency, p_density, p_mod)
plots=list(p_efficiency)

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
      plot.background = element_rect(fill = "black"),
      panel.background = element_rect(fill = "black"), 
      axis.text.x = element_text(size = labsize, colour = "white"),
      plot.title=element_text(size=titlesize),
      axis.text.y = element_text(size = labsize, colour = "white"),
      axis.ticks = element_line(colour = "white"), 
      axis.title.x = element_text(size = titlesize, colour = "white"),
      axis.title.y = element_text(size = titlesize, colour = "white"),
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
    plots[[i]] = plots[[i]] +
      theme_bw()+
      theme(axis.text.x = element_text(size = labsize),
            axis.text.y = element_text(size = labsize),
            axis.title.x = element_blank(),
            plot.title=element_text(size=titlesize),
            axis.ticks = element_line(colour = "black"),
            axis.ticks.length = unit(tick_length, "cm"),
            axis.title.y = element_text(size = titlesize),
            legend.position = "none",
            plot.margin = plot_margin, 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
            axis.line = element_line(size = 0.5))+
      colScale+colScale2  
  }}
dir_save='/GAUDI/plots'
single_figwidth=1.75
single_figheight=1.5
name_list=c('udiam','diam', 'eff', 'den','mod')
name_list=c('eff')

for (i in seq_along(plots)){
  
  dir_save2=paste(dir_save, name_list[[i]],sep='/')
  ggsave(filename = paste(dir_save2, ".png", sep=''), plot = plots[[i]], width = single_figwidth, height = single_figheight, dpi = 800, limitsize = F)}
  
  
p_vol=plots[[1]]; p_eff=plots[[2]]; p_den=plots[[3]]; p_mod=plots[[4]]
p6 <- ggplot(data.frame(), aes(0, 0)) + 
  geom_blank() + 
  theme_void()

if (THEME=='presentation'){
  figwidth= 30.0866
  combined_plot<-cowplot::plot_grid(p_vol,p_eff, ncol = 2,nrow=1,
                                    rel_widths = c(1,1),
                                    align = "hv",
                                    axis = "l",
                                    hspace = 0,  
                                    vspace = 0  
  )
  combined_plot <- combined_plot + 
    theme(plot.background = element_rect(fill = "white"))
  final_plot <- ggdraw(combined_plot)
  ggsave(filename = paste(dir_save, "nest_props_pres.png", sep='/'), plot = final_plot, width = figwidth, height = figheight*1.5, dpi = 1000, limitsize = F)
}
combined_plot<-plot_grid(p_vol, p_eff,p_den, p_mod,   ncol = 2,
                         nrow = 2,
                         rel_widths = c(1,1),
                         align = "hv",
                         axis = "l",
                         hspace = 0,  
                         vspace = 0   
)
combined_plot <- combined_plot + 
  theme(plot.background = element_rect(fill = "white"))
final_plot <- ggdraw(combined_plot) + 
draw_label("Time since treatment (hours)", x =0.52, y = 0.015, angle = 0, size = titlesize)
final_plot
ggsave(filename = paste(dir_save, "nest_network.png", sep='/'), plot = final_plot, width = figwidth, height = figheight*0.4, dpi = 1000, limitsize = F)


