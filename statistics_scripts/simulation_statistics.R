rm(list=ls())

library(lme4)
library(lmerTest)
library(car)
library(plotrix)
library(scales)
library(gamm4)
library(viridis)
library(e1071)
library(splines)

FIX_TIME <-F

test_norm <- function(resids){
  print("Testing normality")
  if (length(resids)<=300){
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
  }else{
    print("More than 300 data points so using the skewness and kurtosis approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}

time_for_stats <- 22600+800*(0:100)

setwd("")#Add correctr dir
model_dir<-''

info <- read.csv("/COLONY_INFO.csv",header=T)
working_directory='/simulation_results/'


dat <- read.csv(paste0(model_dir, 'low_tun4_agg_results_ant_20_2_0.csv'),header = T,stringsAsFactors = F)
dat_distancing <-  read.csv(paste0(model_dir, 'tun4_agg_results_ant_20_2_1.5.csv'),header = T,stringsAsFactors = F)
dat_distancing_low <- read.csv(paste0(model_dir, 'tun4_agg_results_ant_20_2_1.csv'),header = T,stringsAsFactors = F)
dat_distancing_high <-  read.csv(paste0(model_dir, 'tun4_agg_results_ant_20_2_2.csv'),header = T,stringsAsFactors = F)



var='prevalence_hi'#'median_spores'
no_distancing     <- dat
low_distancing    <- dat_distancing_low
medium_distancing <- dat_distancing
high_distancing   <-dat_distancing_high


####edit tables to prepare for analyses
for (what in c("no","low","medium","high")){
  ###get right table
  table <- get(paste(what,"_distancing",sep=""))
  ###add treatment and other metadata info
  table <- merge(table,info[c("name","colony","subset","week","treatment")],all.x=T,all.y=F)
  ###keep only necessary columns for simplification
  table <- table[c("dt","name",'prevalence_hi',"median_spores","colony","subset","week","treatment")]
  ###add column encoding distancing
  table$distancing <- what
  assign(paste(what,"_distancing",sep=""),table)
}
###STATISTICS FOR SYNERGY

library(multcomp)
for (what in c("medium","high","low")){
  print(paste("Statistical analysis -",what,"distancing"))
  
  ###1. combine no_distancing and distancing in single table
  all_data <- rbind(no_distancing,get(paste(what,"_distancing",sep="")))
  all_data<-subset(all_data, dat$dt>=21600)
  time_for_stats <- 22600 + 800 * 0:109     
  all_data <- all_data[all_data$dt %in% time_for_stats, ]
  all_data <- all_data[order(all_data$name,all_data$dt),]
  all_data$prevalence_hi <- all_data$prevalence_hi/180
  
  ###2. Prepare stats - define factors with a specific order of levels, to ensure contrast matrix is ALWAYS correct
  all_data$dt         <- as.numeric(all_data$dt ) ###make sure this is seen as continuous variable
  ##have dt start at 0 and expressed in hours for scaling issues
  all_data$dt         <-  (all_data$dt   - min(all_data$dt ))/3600
  all_data$treatment  <- factor(as.character( all_data$treatment ),levels=c("SHAM","PATHOGEN"))
  all_data$distancing <- factor(as.character( all_data$distancing ),levels=c("no",what))
  
  ###4. Fit model
  if (var=='median_spores'){
  model <- lmer(sqrt(median_spores+1) ~ treatment *distancing*dt  + (1|colony/subset) + (1|week)  ,data=all_data)
  }else{
    model <- lmer(prevalence_hi ~ treatment *distancing*dt +  (1|colony/subset) + (1|week)  ,data=all_data)
    
  }
  ###5. Check residuals
  test_norm(residuals(model))
  
  ###6. Test model significance
  print(Anova(model,type=3))
  
  ###contrast matrix
  contrast_matrix <- rbind(
    "Slope_DistancingEffect_SHAM"=c(0,0,0,0,0,0,1,0)
    ,
    "Slope_DistancingEffect_PATHOGEN"=c(0,0,0,0,0,0,1,1)
    ,
    "Slope_Synergy"=c(0,0,0,0,0,0,0,1)
  )
  print(summary(glht(model,linfct=contrast_matrix),test=adjusted("BH")))
}

##MAIN
###If we want to create so there is self-isolation only in pathogen
dat_distancing_pathogen_only <- rbind(
  subset(dat,            treatment == "SHAM"),
  subset(dat_distancing, treatment == "PATHOGEN")
)

#Reassign dat as desired, depending on if self-isolation should be included
dat<-dat_distancing_pathogen_only
dat<-subset(dat, dat$dt>=21600)
time_for_stats <- 22600 + 800 * 0:109     
dat <- dat[dat$dt %in% time_for_stats, ]
dat <- dat[order(dat$name,dat$dt),]
FIX_TIME=F
if(FIX_TIME){
  dat[which(dat$dt>=26600),"dt"] <-    dat[which(dat$dt>=26600),"dt"]+2000
  dat[which(dat$dt<26600),"dt"]  <-    dat[which(dat$dt<26600),"dt"]+40
  }

dat <- merge(dat,info[c("name", "colony", "week", "subset", 'day', 'treatment')],all.x=T,all.y=F)

dat$condition<-dat$treatment
dat$condition <- factor(dat$condition , levels = c("PATHOGEN","SHAM"))


variables <- c("prevalence_hi",'median_spores')
transformations <- c("sqrt", 'power2')#,'sqrt')#power2 for main
names(variables) <- c("prevalence_hi", 'median_spores')
names(transformations) <- variables


for (variable in  variables){
  print(variable)
  
  dat$variable <- dat[,variable]
  
  if (transformations[variable]=="log10"){
    if (min(dat$variable ,na.rm=T)>0){
      dat$variable <- log10(dat$variable )
    }else{
      dat$variable <- log10(dat$variable - min(0, dat$variable ) + min( dat$variable [which(dat$variable >0)]  , na.rm=T)/sqrt(2))
    }
  }else if (transformations[variable]=="sqrt"){
    if (min(dat$variable ,na.rm=T)>0){
      dat$variable <- sqrt(dat$variable )
    }else{
      dat$variable <- sqrt(dat$variable - min(dat$variable ,na.rm=T) )
      
    }
  }else if (grepl("power",transformations[variable])){
    if (min(dat$variable ,na.rm=T)>0){
      dat$variable <- (dat$variable )^as.numeric(gsub("power","",transformations[variable]))
    }else{
      dat$variable <- (dat$variable - min(dat$variable ,na.rm=T))^as.numeric(gsub("power","",transformations[variable]))
      
    }
    
  }
  if (!all(!is.na(dat$variable ))){
    print(paste("NaN produced for", variable))
  }

    model <- lmer(variable ~ condition*dt + (1|colony/subset) + (1|week),control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)) ,data=dat)
    print(summary(model))
    test_norm(residuals(model))
    print(Anova(model), type=3)
    
}
