#GERMINATION DATA
#VIVIAN BERNAU
#31 AUGUST 2017

#Based on McNair et al 2012; Seed Science Research
#Check on Allison 2010 (Survival analysis using SAS)

#SEMIPARAMETRIC ANALYSIS
  #allows assessment of multiple covariates: quantitative & categorical; random & fixed
  #random effects are referred to as frailty effects--these are most common as subgroups within treatment groups 
      #(e.g. reps, maybe regions?)

#Hazard function: how likely it is that a seed which has not germinated by time t will germinate shortly after t.

#Set working directory and repositories
wd <- ("~/Google Drive/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(out.dir)

#read in germination data in pre-lifetab format
df <- read.csv(paste(out.dir, "/cleaned_2017-09-07.csv", sep = ""), header = T)
str(df)
df$end <- round(df$end, digits = 1)
df$run <- as.factor(df$run)

#re-order data by region
summary(df$region)
target <- c("central valleys", "ecoast", "sierra madre", "wcoast", "yucatan", "control")
library(gdata)
df$region <- reorder.factor(df$region, new.order=target)
summary(df$region)
head(df)

df$trt <- as.factor(df$trt)
str(df$trt)

library(survival)
library(survminer)
fit <- coxph(Surv(end, status) ~ region, data = df)
ftest <- cox.zph(fit)
ftest
#rho = pearson correlation with residuals and log(time) for each covariate
#p-value less than 0.05 indicates violation of proportionality

ggcoxzph(ftest) #test proportional hazards assumption

#plot residuals
ggcoxdiagnostics(fit, type = "schoenfeld", ox.scale = "time")
ggcoxdiagnostics(fit, type = "deviance", ox.scale = "linear.predictions")

fit0 <- coxph(Surv(end, status) ~ region + trt, data = df)
ftest0 <- cox.zph(fit0)
fit0
ggforest(fit0)

#ADD COVARIATES TO THE COX MODEL
#test distribution of simple model: exponential, weibull, lognormal, loglogistic
#create a cox model with covariates
  #1) trt
  #2) region
  #3) trt*region
  #4) line/region
  #5) trt*line/region
  #6) frailty(run)
survreg <- survreg(Surv(end, status) ~ trt + region + trt:region + trt:line%in%region + frailty(run) + frailty(plate%in%run), data = df1, dist = "weibull")

survreg_run <- survreg(Surv(end, status) ~ factor(trt) + factor(region)  + trt:region  + frailty(run), data = df, dist = "exponential")
survreg_run2 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df, dist = "weibull")
survreg_run3 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df, dist = "lognormal")
survreg_run4 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df, dist = "loglogistic")

anova(survreg_run,survreg_run2,survreg_run3,survreg_run4)
extractAIC(survreg_run);extractAIC(survreg_run2);extractAIC(survreg_run3);extractAIC(survreg_run4)

#Parameter estimation 
survreg_run3



