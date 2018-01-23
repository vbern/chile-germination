#Germination Analysis
#Vivian Bernau
#11 Oct 2016
#18 Oct 2016 modified with input from Onofri et al. 2010 and Fox 2001

#Based on McNair et al 2012; Seed Science Research

#My data is right-sensored.
#My data is interval, but can probably be analyzed as exact.
#Semi-parametric time to event/failure-time/reliability analysis is most statistically sound for germination data.
#Time required for germination is a random variable. 

#Survivor function: probability that the germination time is greater than t
#Hazard function: how likely it is that a seed which has not germinated by time t will germinate shortly after t.

#N (number of seeds at risk of germination) should be adjusted for seed losses which occur before the end of the experiment.
#   This in cludes seeds damaged by handling, mold growth, lost seeds, etc.  Should not exceed 5%

#Random effects are referred to as frailty effects--these are most common as subgroups within treatment groups (e.g. reps, maybe regions?)

library("survival")
##survfit() can be used to compute the Kalplan-Meier survivor function

library("gtools")
library("drc")
library("alr3")
library("m")

#Set working directory and repositories
wd <- ("~/Google Drive/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

#Plotting nonparametric estimates of survivor fucntions for different treatement groups is a useful exploratory technique.

#read in germination data in pre-lifetab format
setwd(data.dir)

###########
data1 <- read.csv("20161012_germdata_rep1_modern.csv", na.strings=c("","NA"))
data1 <- data1[1:10]
str(data1)

data2 <- read.csv("20161012_germdata_rep2_modern.csv", na.strings=c("","NA"))
data2 <- data2[1:10]
str(data2)

data3_run5 <- read.csv("20161012_germdata_run5_modern_0sAtEnd.csv", na.strings = c("", "NA"))
data3_run5 <- data3_run5[1:10]
str(data3_run5)

data <- rbind(data1, data2, data3_run5)
str(data)

#converts old data style to "modern" structure as described by McNair et al. 2012
###VERY IMPORTANT: Creates multiples of each "plate" for the number of seeds germinated at each time point
data.mod <- data[rep(row.names(data), data$germinated), 1:10]
str(data.mod)

write.csv(data.mod, file = paste(out.dir, "/data_mod_",Sys.Date(),".csv", sep = ""))

sampleID <- read.csv("seedincrease_sampleid.csv")
colnames(sampleID) <- tolower(colnames(sampleID))

data1 <- merge(data.mod, sampleID, by.x = "line", by.y = "line", all.x = T, ally.y = F)

collections <- read.csv("2013_Collections.csv")
colnames(collections) <- tolower(colnames(collections))
collections <- collections[, (colnames(collections) %in% c("sample.id", "pop", "region", "population.type", 
                                                           "landrace.name", "cultivation", "main.use"))]
data2 <- merge(data1, collections, by.x = "sampleid", by.y = "sample.id", all.x = T, all.y = F, row.names = "pedigree")

write.csv(data2, file = paste(out.dir, "/alldata_mod_", Sys.Date(), ".csv", sep = ""))

head(data2)

#must go through and edit all of the Johnny's records so that region appears

###########################

df <- read.csv(paste(data.dir, "/alldata_mod_2016-10-31.csv", sep = ""), header = T) ##rep 1 and 2, run 5
str(df)
round(df$end, digits = 1)
df$run <- as.factor(df$run)

summary(df$region)
target <- c("Johnny's", "central valleys", "ecoast", "sierra madre", "wcoast", "yucatan")

require(gdata)
df$region <- reorder.factor(df$region, new.order=target)
summary(df$region)
head(df)
require(dplyr)
df1<-df %>% arrange(region)
head(df1)
summary(df$region)

#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df1, type = "kaplan-meier")
plot(km.fit, lwd = 2, col = "blue", conf.int = T, xlab = "time (hours)", ylab = "probability of not germinating", xlim = c(0,550))

test.run.1 <-survfit(Surv(end, status) ~ run, data = df1, type = "kaplan-meier")
summary(test.run.1)
plot(test.run.1, col = rainbow(5),lwd = 1.75, xlab = "time (hours)", 
     ylab = "probability of not germinating", xlim = c(0,550), main = "Germination by Run")
legend(5, .2, c(1:5), col = rainbow(5), bty="n", cex = 1, lwd = 1.75)

df1 <- subset(df, df$run ==1 | df$run ==2 |df$run ==3| df$run ==5)
str(df1)

test.trt.1 <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", conf.type = "log-log")
plot(test.trt.1, col = "black", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), main = "Germination by PEG concentration")
legend(5, .1, c("  0% PEG", "10% PEG", "15% PEG", "20% PEG"), lty = 1:4,
       col = "black", lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

test.region.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier", conf.type = "log-log")
plot(test.region.1, col = c("grey","black", "red", "green", "cyan", "blue"), 
     lwd = 1.40, xlim = c(0,550), 
     main = "Germination by region"
)
legend(5, .35, c("Johnny's", "Central Valleys", "E Coast", "Sierra Madre de la Sur", "W Coast", "Yucatan" ), 
       col =  c("grey","black", "red", "green", "cyan", "blue"),lwd = 1.40, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

test.landrace.1 <- survfit(Surv(end, status) ~ landrace.name, data = df1, type = "kaplan-meier")
plot(test.landrace.1, col = rainbow(18), lwd = 1.75, xlim = c(0,550))
legend(5, .5, landraces, 
       col =  rainbow(18), lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

test.shelf.1 <- survfit(Surv(end, status) ~ shelf, data = df1, type = "kaplan-meier")
plot(test.shelf.1, col = rainbow(9), lwd = 1.4, xlim = c(0,550))
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

test.cultivation.1 <- survfit(Surv(end, status) ~ cultivation, data = df1, type = "kaplan-meier")
plot(test.cultivation.1, col = rainbow(5), lty = 1:5, lwd = 1.75, xlim = c(0,550))
legend(5, .5, c("backyard", "forest", "milpa", "plantation", "NA"), 
       col = rainbow(5), lty = 1:5, lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

#define function mlogmlog() to calculate -log(-log(S(t))) checks proportional hazards (PH) assumption
mlogmlog <- function(y){-log(-log(y))}

#this model is reasonably robust, so only worry about clear departures from proportionality, 
#as indiecated by decisive crossing of the functions for two or more coariates in the diagnostic plot
plot(test.trt.1, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     col = c("green", "blue", "yellow", "red")) #looks good!

plot(test.run.1, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     lty = c("solid", "longdash", "solid", "longdash"), col = c("grey", "grey", "black","black")) #hazard lines cross

plot(test.region.1, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     col = c("grey","black", "red", "green", "cyan", "blue")) #hazard lines converge and cross

#log rank test
survdiff_20_logrank <- survdiff(Surv(end, status)~region, data = df1, subset = {trt==20  & region!= "wcoast"}, rho = 0) #log rank test
print(survdiff_20_logrank)

survdiff_cv <- survdiff(Surv(end, status)~trt, data = df1, subset = {region == "central valleys"}, rho = 0)
print(survdiff_cv)

survdiff_ecoast <- survdiff(Surv(end, status)~trt, data = df1, subset = {region == "ecoast"}, rho = 0)
print(survdiff_ecoast)

survdiff_wcoast <- survdiff(Surv(end, status)~trt, data = df1, subset = {region == "wcoast"}, rho = 0)
print(survdiff_wcoast)

survdiff_sm <- survdiff(Surv(end, status)~trt, data = df1, subset = {region == "sierra madre"}, rho = 0)
print(survdiff_sm)

survdiff_yucatan <- survdiff(Surv(end, status)~trt, data = df1, subset = {region == "yucatan"}, rho = 0)
print(survdiff_yucatan)

setwd(out.dir)

#Creating composite figures of with one frame for each [PEG]
png(filename="allconcentrations.png", width = 700, height = 1400, res = 120)
plot.new()
par(mfrow=c(4,1),mar = c(0,4,0,4)+.1)
par(oma=c(5,0,3,0))

test.20.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier", subset = (trt ==20))
plot(test.20.1, col = c("grey","black", "red", "green", "cyan", "blue"), lty = 4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
    )
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "20% PEG", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.15.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier", subset = (trt ==15))
plot(test.15.1, col = c("grey","black", "red", "green", "cyan", "blue"), lty = 3,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "15% PEG", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.10.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier", subset = (trt ==10))
plot(test.10.1, col = c("grey","black", "red", "green", "cyan", "blue"), lty = 2, 
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "10% PEG", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.0.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier", subset = (trt ==0))
plot(test.0.1, col = c("grey","black", "red", "green", "cyan", "blue"), lty = 1,
     lwd = 1.75, xlim = c(0,550), xlab = "time (hours)")
text(450,.9, "0% PEG", cex = 2)
axis(4)
title(line = 2, ylab = "probability of not germinating")
legend(5, .35, c("Johnny's", "Central Valleys", "E Coast", "Sierra Madre de la Sur", "W Coast", "Yucatan" ), 
       col =  c("grey","black", "red", "green", "cyan", "blue"),lwd = 1.40, bty = "n", cex = 1)
mtext(text="time (hours)",side=1,line=3,outer=TRUE)
dev.off()

#Creating composite figure with 1 frame for each region
png(filename="allregions.png", width = 700, height = 1400, res = 120)
plot.new()
par(mfrow=c(5,1),mar = c(0,4,0,4)+.1)
par(oma=c(5,0,3,0))

test.cv <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="central valleys"))
plot(test.cv, col = "black", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "Central Valleys", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.ecoast <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="ecoast"))
plot(test.ecoast, col = "red", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "E Coast", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.sm <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="sierra madre"))
plot(test.sm, col = "green", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "Sierra Madre", cex = 2)
text(450,.8, "del Sur", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.wcoast <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="wcoast"))
plot(test.wcoast, col = "cyan", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "W Coast", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.yucatan <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="yucatan"))
plot(test.yucatan, col = "blue", lty = 1:4,
     lwd = 1.75, xlim = c(0,550))
text(450,.9, "Yucatan", cex = 2)
legend(20, .7, c("  0% PEG", "10% PEG", "15% PEG", "20% PEG"), lty = 1:4,
       col = "black", lwd = 1.75, bty = "n", cex = 1.5)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating", cex.axis = 1.5)
axis(4)
mtext(text="time (hours)",side=1,line=3,outer=TRUE)
dev.off()

#survreg <- survreg(Surv(end, status) ~ trt + region + trt:region + trt:line%in%region + frailty(run) + frailty(plate%in%run), data = df1, dist = "weibull")

survreg_run <- survreg(Surv(end, status) ~ factor(trt) + factor(region)  + trt:region  + frailty(run), data = df1, dist = "exponential")
survreg_run2 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df1, dist = "weibull")
survreg_run3 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df1, dist = "lognormal")
survreg_run4 <- survreg(Surv(end, status) ~ factor(trt) + factor(region) + trt:region + frailty(run), data = df1, dist = "loglogistic")

anova(survreg_run,survreg_run2,survreg_run3,survreg_run4)
extractAIC(survreg_run);extractAIC(survreg_run2);extractAIC(survreg_run3);extractAIC(survreg_run4)

#Parameter estimation 
##Question: Is this where you estimate the degree of your continuous parameters?
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

