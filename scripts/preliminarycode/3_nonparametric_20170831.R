#GERMINATION DATA
#VIVIAN BERNAU
#31 AUGUST 2017

#Based on McNair et al 2012; Seed Science Research

#NONPARAMETRIC ANALYSIS
  #1) Characterizing pattern of germination within groups --> Kaplan-Meier test: survfit()
  #2) Comparing patterns of germination in groups --> Fleming-Harrington test(?): survdiff()

#My data is right-sensored.
#My data is interval, but can probably be analyzed nonparametrically as exact for as long as plates with large losses are removed.

#Survivor function: probability that the germination time is greater than t

#Set working directory and repositories
wd <- ("~/Google Drive/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(out.dir)

#read in germination data in pre-lifetab format
df <- read.csv(paste(out.dir, "/alldata_mod_2017-09-04.csv", sep = ""), header = T)
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
#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
par(mfrow=c(1,1),mar = c(0,4,0,4)+.1)
km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df, type = "kaplan-meier")
plot(km.fit, lwd = 2, col = "blue", conf.int = T, xlab = "time (hours)", ylab = "probability of not germinating", xlim = c(0,550))

test.region <- survfit(Surv(end, status) ~ region, data = df, subset = {trt == 0}, type = "kaplan-meier", conf.type = "log-log")
jpeg('test_region_20170831.jpg')
plot(test.region, col = c("grey","black", "red", "green", "cyan", "blue"), 
     lwd = 1.40, xlim = c(0,550), 
     main = "Germination by region")
legend(5, .35, c("Johnny's", "Central Valleys", "E Coast", "Sierra Madre de la Sur", "W Coast", "Yucatan" ), 
       col =  c("grey","black", "red", "green", "cyan", "blue"),lwd = 1.40, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")
dev.off()
survdiff(Surv(end,status)~region, data = df, subset = {trt == 0}, na.omit, rho = 0)

#Creating composite figures of with one frame for each [PEG]
png(filename="allconcentrations.png", width = 700, height = 1400, res = 120)
plot.new()
par(mfrow=c(4,1),mar = c(0,4,0,4)+.1)
par(oma=c(5,0,3,0))

test.20<- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier", subset = (trt ==20))
plot(test.20, col = c("grey","black", "red", "green", "cyan", "blue"), lty = 4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "20% PEG", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.15 <- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier", subset = (trt ==15))
plot(test.15, col = c("grey","black", "red", "green", "cyan", "blue"), lty = 3,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "15% PEG", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.10 <- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier", subset = (trt ==10))
plot(test.10, col = c("grey","black", "red", "green", "cyan", "blue"), lty = 2, 
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "10% PEG", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.0 <- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier", subset = (trt ==0))
plot(test.0, col = c("grey","black", "red", "green", "cyan", "blue"), lty = 1,
     lwd = 1.75, xlim = c(0,550), xlab = "time (hours)")
text(450,.9, "0% PEG", cex = 2)
axis(4)
title(line = 2, ylab = "probability of not germinating")
legend(5, .35, c("Johnny's", "Central Valleys", "E Coast", "Sierra Madre de la Sur", "W Coast", "Yucatan" ), 
       col =  c("grey","black", "red", "green", "cyan", "blue"),lwd = 1.40, bty = "n", cex = 1)
mtext(text="time (hours)",side=1,line=3,outer=TRUE)
dev.off()

#Creating composite figure of Kaplan-Meier survivor estimates with 1 frame for each region
png(filename="allregions.png", width = 700, height = 1400, res = 120)
plot.new()
par(mfrow=c(5,1),mar = c(0,4,0,4)+.1)
par(oma=c(5,0,3,0))

test.cv <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="central valleys"))
plot(test.cv, col = "black", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "Central Valleys", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.ecoast <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="ecoast"))
plot(test.ecoast, col = "red", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "E Coast", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.sm <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="sierra madre"))
plot(test.sm, col = "green", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "Sierra Madre", cex = 2)
text(450,.8, "del Sur", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.wcoast <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="wcoast"))
plot(test.wcoast, col = "cyan", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), xaxt = "n"
)
axis(4)
axis(1, at=c(100,200,300,400,500), labels=c("","","","",""))
text(450,.9, "W Coast", cex = 2)
title(line = 2, ylab = "probability of not germinating")

test.yucatan <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="yucatan"))
plot(test.yucatan, col = "blue", lty = 1:4,
     lwd = 1.75, xlim = c(0,550))
text(450,.9, "Yucatan", cex = 2)
legend(20, .7, c("  0% PEG", "10% PEG", "15% PEG", "20% PEG"), lty = 1:4,
       col = "black", lwd = 1.75, bty = "n", cex = 1.5)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating", cex.axis = 1.5)
axis(4)
mtext(text="time (hours)",side=1,line=3,outer=TRUE)
dev.off()

par(mfrow=c(1,1),mar = c(0,4,0,4)+.1)

#COMPARISON OF SURVIVOR FUNCTIONS
#log rank tests -- HAVING PROBLEMS HERE
        #Error in survdiff.fit(y, groups, strata.keep, rho) : 
        #NA/NaN/Inf in foreign function call (arg 5)
survdiff_cv <- survdiff(Surv(end, status)~factor(trt), data = df, subset = {region == "central valleys"}, rho = 0)
print(survdiff_cv)

survdiff_ecoast <- survdiff(Surv(end, status)~trt, data = df, subset = {region == "ecoast"}, rho = 0)
print(survdiff_ecoast)

survdiff_wcoast <- survdiff(Surv(end, status)~trt, data = df, subset = {region == "wcoast"}, rho = 0)
print(survdiff_wcoast)

survdiff_sm <- survdiff(Surv(end, status)~trt, data = df, subset = {region == "sierra madre"}, rho = 0)
print(survdiff_sm)

survdiff_yucatan <- survdiff(Surv(end, status)~trt, data = df, subset = {region == "yucatan"}, rho = 0)
print(survdiff_yucatan)