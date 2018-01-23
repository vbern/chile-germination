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
library(survival)

#define function mlogmlog() to calculate ~log(~log(S(t)))
mlogmlog <- function(y){-log(log(y))}

#estimate kaplan-meier survivor functions for each of the three sites, for all collection dates combined
fit.line <- survfit(Surv(end, status) ~ line, type = "kaplan-meier", data=df)

#plot -log(-log(S(t))) v log(t)
plot(fit.line, mark.time = F, fun=mlogmlog, log="x", xlab = "t [hours]", ylab = "log(-log(S(t)))", 
     lty = c("solid", "longdash", "dotted"), lwd = 1.75)

#CHECK FOR MULTICOLLINEARITY
    #definition:
library(HH)
#fit a generalized linear model
multicollinearitycheck <- glm(end ~ trt + region, data = df)
#check for collinearity among covariates
vif <- vif(collinearitycheck)
print(vif)

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

survreg_run <- survreg(Surv(end, status) ~ factor(trt) + factor(region)  + trt:region  + frailty(run), data = df1, dist = "exponential")
survreg_run2 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df1, dist = "weibull")
survreg_run3 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df1, dist = "lognormal")
survreg_run4 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df1, dist = "loglogistic")

anova(survreg_run,survreg_run2,survreg_run3,survreg_run4)
extractAIC(survreg_run);extractAIC(survreg_run2);extractAIC(survreg_run3);extractAIC(survreg_run4)

#Parameter estimation 
survreg_run3

#Extension to regression
surv.reg <- survreg(Surv(end, status) ~ (trt) + region + trt:region + cluster(run) + subcluster(plate), 
                    data=df1, dist="lognormal", na.action = na.omit)
summary(surv.reg)
surv.reg

out <- capture.output(surv.reg)
cat("RegressionOutput", out, file = paste(out.dir, "/survreg", Sys.Date(), ".txt", sep = ""))


cv0<- c(1,0,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0)
cv10<- c(1,10,1,0,0,0,0,1,1,1,1,1,1,10,0,0,0,0)
cv15<- c(1,15,1,0,0,0,0,1,1,1,1,1,1,15,0,0,0,0)
cv20<- c(1,20,1,0,0,0,0,1,1,1,1,1,1,20,0,0,0,0)
ec0<- c(1,0,0,1,0,0,0,1,1,1,1,1,1,0,0,0,0,0)
ec10<- c(1,10,0,1,0,0,0,1,1,1,1,1,1,0,10,0,0,0)
ec15<- c(1,15,0,1,0,0,0,1,1,1,1,1,1,0,15,0,0,0)
ec20<- c(1,20,0,1,0,0,0,1,1,1,1,1,1,0,20,0,0,0)
sm0<- c(1,0,0,0,1,0,0,1,1,1,1,1,1,0,0,0,0,0)
sm10<- c(1,10,0,0,1,0,0,1,1,1,1,1,1,0,0,10,0,0)
sm15<- c(1,15,0,0,1,0,0,1,1,1,1,1,1,0,0,15,0,0)
sm20<- c(1,20,0,0,1,0,0,1,1,1,1,1,1,0,0,20,0,0)
wc0<- c(1,0,0,0,0,1,0,1,1,1,1,1,1,0,0,0,0,0)
wc10<- c(1,10,0,0,0,1,0,1,1,1,1,1,1,0,0,0,10,0)
wc15<- c(1,15,0,0,0,1,0,1,1,1,1,1,1,0,0,0,15,0)
wc20<- c(1,20,0,0,0,1,0,1,1,1,1,1,1,0,0,0,20,0)
yuc0<- c(1,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0)
yuc10<- c(1,10,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,10)
yuc15<- c(1,15,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,15)
yuc20<- c(1,20,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,20)


t <- seq(0,550,0.5)
cv0 <- 1-pnorm((log(t)-cv0%*%surv.reg$coefficients)/surv.reg$scale)
cv10 <- 1-pnorm((log(t)-cv10%*%surv.reg$coefficients)/surv.reg$scale)
cv15<- 1-pnorm((log(t)-cv15%*%surv.reg$coeff)/surv.reg$scale)
cv20<- 1-pnorm((log(t)-cv20%*%surv.reg$coeff)/surv.reg$scale)
plot(test.cv, col = "black", lty = 1:4,
     lwd = 1.75, xlim = c(0,550)
)
#plot(t,cv0,type='l',xlab="Time",ylab="Surv Prob",ylim=c(0,1))
lines(t,cv0,lty=1)
lines(t,cv10,lty=2)
lines(t,cv15,lty=3)
lines(t,cv20,lty=4)

ec0 <- 1-pnorm((log(t)-ec0%*%surv.reg$coefficients)/surv.reg$scale)
ec10 <- 1-pnorm((log(t)-ec10%*%surv.reg$coefficients)/surv.reg$scale)
ec15<- 1-pnorm((log(t)-ec15%*%surv.reg$coeff)/surv.reg$scale)
ec20<- 1-pnorm((log(t)-ec20%*%surv.reg$coeff)/surv.reg$scale)
plot(test.ecoast, col = "red", lty = 1:4,
     lwd = 1.75, xlim = c(0,550)
)
#plot(t,ec0,type='l',xlab="Time",ylab="Surv Prob",ylim=c(0,1))
lines(t,ec0,lty=1)
lines(t,ec10,lty=2)
lines(t,ec15,lty=3)
lines(t,ec20,lty=4)

sm0 <- 1-pnorm((log(t)-sm0%*%surv.reg$coefficients)/surv.reg$scale)
sm10 <- 1-pnorm((log(t)-sm10%*%surv.reg$coefficients)/surv.reg$scale)
sm15<- 1-pnorm((log(t)-sm15%*%surv.reg$coeff)/surv.reg$scale)
sm20<- 1-pnorm((log(t)-sm20%*%surv.reg$coeff)/surv.reg$scale)
plot(test.sm, col = "green", lty = 1:4,
     lwd = 1.75, xlim = c(0,550)
)
#plot(t,sm0,type='l',xlab="Time",ylab="Surv Prob",ylim=c(0,1))
lines(t,sm0,lty=1)
lines(t,sm10,lty=2)
lines(t,sm15,lty=3)
lines(t,sm20,lty=4)

wc0 <- 1-pnorm((log(t)-wc0%*%surv.reg$coefficients)/surv.reg$scale)
wc10 <- 1-pnorm((log(t)-wc10%*%surv.reg$coefficients)/surv.reg$scale)
wc15<- 1-pnorm((log(t)-wc15%*%surv.reg$coeff)/surv.reg$scale)
wc20<- 1-pnorm((log(t)-wc20%*%surv.reg$coeff)/surv.reg$scale)
plot(test.wcoast, col = "cyan", lty = 1:4,
     lwd = 1.75, xlim = c(0,550)
)#plot(t,wc0,type='l',xlab="Time",ylab="Surv Prob",ylim=c(0,1))
lines(t,wc0,lty=1)
lines(t,wc10,lty=2)
lines(t,wc15,lty=3)
lines(t,wc20,lty=4)

yuc0 <- 1-pnorm((log(t)-yuc0%*%surv.reg$coefficients)/surv.reg$scale)
yuc10 <- 1-pnorm((log(t)-yuc10%*%surv.reg$coefficients)/surv.reg$scale)
yuc15<- 1-pnorm((log(t)-yuc15%*%surv.reg$coeff)/surv.reg$scale)
yuc20<- 1-pnorm((log(t)-yuc20%*%surv.reg$coeff)/surv.reg$scale)
plot(test.yucatan, col = "blue", lty = 1:4,
     lwd = 1.75, xlim = c(0,550)
)
#plot(t,yuc0,type='l',xlab="Time",ylab="Surv Prob",ylim=c(0,1))
lines(t,yuc0,lty=1)
lines(t,yuc10,lty=2)
lines(t,yuc15,lty=3)
lines(t,yuc20,lty=4)


