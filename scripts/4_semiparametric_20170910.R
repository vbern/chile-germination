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
df <- read.csv(paste(out.dir, "/cleaned_2017-09-17.csv", sep = ""), header = T)
str(df)
#df$end <- round(df$end, digits = 1)
df$run <- as.factor(df$run)
df<- subset(df, region != "control")
#df <- subset(df, trt != 15)
#df <- subset(df, trt != 10)
hist(df$trt)

df$cv <- as.numeric(df$region == "central valleys")
df$ecoast <- as.numeric(df$region == "ecoast")
df$wcoast <- as.numeric(df$region == "wcoast")
df$yucatan <- as.numeric(df$region == "yucatan")
df$sm <- as.numeric(df$region == "sierra maadre")

library(survival)
library(survminer)

df1 <- subset(df, status == 1)

#CHECK PROPORATIONAL HAZARDS ASSUMPTION
#define function mlogmlog() to calculate ~log(~log(S(t)))
mlogmlog <- function(y){log(-log(y))}
#estimate kaplan-meier survivor functions for each of the three sites, for all collection dates combined
fit_region <- survfit(Surv(end, status) ~ region, data=df, na.action = na.omit)
fit_landrace <- survfit(Surv(end, status) ~ landrace.name, data = df, na.action = na.omit)
#plot -log(-log(S(t))) v log(t)
plot(fit_region, fun=mlogmlog, log="x", xlab = "t [hours]", ylab = "-log(-log(S(t)))", lwd = 1.75, col = rainbow(5))
plot(fit_landrace, fun=mlogmlog, log="x", xlab = "t [hours]", ylab = "-log(-log(S(t)))", lwd = 1.75, col = rainbow(12))

fit_region <- coxph(Surv(end, status) ~ region*end, data = df)
ftest_region <- cox.zph(fit_region)
ftest_region
plot(ftest_region)

fit_trt <- coxph(Surv(end, status) ~ factor(trt), data = df)
ftest_trt <- cox.zph(fit_trt)
ftest_trt
#rho = pearson correlation with residuals and log(time) for each covariate
#p-value less than 0.05 indicates violation of proportionality

#PLOT RESIDUALS
fit0 <- coxph(Surv(end, status) ~ region + trt + region:trt, data = df)
ftest0 <- cox.zph(fit0)
fit0
ggcoxdiagnostics(fit0, type = "schoenfeld", ox.scale = "time")
ggcoxdiagnostics(fit0, type = "martingale", ox.scale = "linear.predictions")

#Plot continuous explanatory variable against martingale residuals of null cox proportional hazards model
#Helps to choose the proper functional form of continuous variable in the cox model
a <- coxph(Surv(end, status) ~ trt + trt^2, data = df)
ggcoxfunctional(a, data = df, point.col = "blue", point.alpha = 0.5)

#ADD COVARIATES TO THE COX MODEL
#test distribution of simple model: exponential, weibull, lognormal, loglogistic
#create a cox model with covariates
  #1) trt
  #2) region
  #3) trt*region
  #4) line/region
  #5) frailty(plate%in%run%in%rep)

a <- survreg(Surv(end, status) ~ trt + region + frailty(plate%in%run%in%rep), data = df, dist = "loglogistic")
b <- survreg(Surv(end, status) ~ trt + region + trt:region + frailty(plate%in%run%in%rep), data = df, dist = "loglogistic")
c <- survreg(Surv(end, status) ~ trt + region + trt:region + line%in%region + frailty(plate%in%run%in%rep), data = df, dist = "loglogistic")
#d <- survreg(Surv(end, status) ~ trt + region + trt:region + frailty(plate%in%run%in%rep), data = df, dist = "loglogistic")
anova(a, b)
anova(b, c)
extractAIC(a); extractAIC(b); extractAIC(c)

#FIT DISTRIBUTION TYPE
survreg_run <- survreg(Surv(end, status) ~ (trt) + region + (trt):region + frailty(landrace.name) + frailty(plate%in%run%in%rep) , data = df, dist = "exponential")
survreg_run2 <- survreg(Surv(end, status) ~ (trt) + region + (trt):region + landrace.name + frailty(line%in%region) + frailty(plate%in%run%in%rep), data = df, dist = "weibull")
survreg_run3 <- survreg(Surv(end, status) ~ (trt) + region + (trt):region + landrace.name +frailty(line%in%region) + frailty(plate%in%run%in%rep), data = df, dist = "lognormal")
survreg_run4 <- survreg(Surv(end, status) ~ (trt) + region + (trt):region + frailty(plate%in%run%in%rep), data = df, dist = "loglogistic")

anova(survreg_run,survreg_run2,survreg_run3,survreg_run4)
extractAIC(survreg_run);extractAIC(survreg_run2);extractAIC(survreg_run3);extractAIC(survreg_run4) #choose loglogistic distribution, lowest AIC

#PARAMETER ESTIMATION
summary(survreg_run4)
survreg_run4
out <- capture.output(survreg_run4)
cat("RegressionOutput", out, file = paste(out.dir, "/survreg", Sys.Date(), ".txt", sep = ""))

#BUILD MODEL IN ORDER OF DECREASING SIGNIFICANCE FROM PARAMETER ESTIMATION
a <- survreg(Surv(end, status) ~ trt + frailty(plate%in%run%in%rep) , data = df)
b <- survreg(Surv(end, status) ~ trt + I(trt^2) + frailty(plate%in%run%in%rep), data = df)
anova(a,b)
extractAIC(a); extractAIC(b)

c <- survreg(Surv(end, status) ~ trt + I(trt^2) + ecoast + frailty(plate%in%run%in%rep), data = df )
anova(b, c)
extractAIC(b); extractAIC(c)

d <- survreg(Surv(end, status) ~ trt + I(trt^2) + ecoast + sm + frailty(plate%in%run%in%rep), data = df)
anova(c,d)
extractAIC(c); extractAIC(d)

e <- survreg(Surv(end, status) ~ trt + I(trt^2) + ecoast + yucatan + frailty(plate%in%run%in%rep), data = df)
anova(c,e)
extractAIC(c); extractAIC(e)

f <- survreg(Surv(end, status) ~ trt + I(trt^2) + ecoast + yucatan + wcoast + frailty(plate%in%run%in%rep), data = df)
anova(e, f)
extractAIC(e); extractAIC(f)

g <- survreg(Surv(end, status) ~ trt + I(trt^2) + ecoast + yucatan + wcoast + trt:ecoast + I(trt^2):ecoast+ frailty(plate%in%run%in%rep), data = df)
anova(f, g)
extractAIC(f); extractAIC(g)

h <- survreg(Surv(end, status) ~ trt + I(trt^2) + ecoast + yucatan + wcoast + trt:ecoast + I(trt^2):ecoast + trt:wcoast + I(trt^2):wcoast + frailty(plate%in%run%in%rep), data = df)
anova(g, h)
extractAIC(g); extractAIC(h)

i <- survreg(Surv(end, status) ~ trt + I(trt^2) + ecoast + yucatan + wcoast + trt:ecoast + I(trt^2):ecoast + trt:wcoast + I(trt^2):wcoast + trt:yucatan + I(trt^2):yucatan+ frailty(plate%in%run%in%rep), data = df)
anova(h, i)
extractAIC(h); extractAIC(i)

j <- survreg(Surv(end, status) ~ trt + I(trt^2) + ecoast + yucatan + wcoast + trt:ecoast + I(trt^2):ecoast + trt:wcoast + I(trt^2):wcoast + trt:sm + I(trt^2):sm + frailty(plate%in%run%in%rep), data = df)
anova(h, j)
extractAIC(h); extractAIC(j)
