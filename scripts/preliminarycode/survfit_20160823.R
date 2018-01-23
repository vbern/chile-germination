#Germination Analysis
#Vivian Bernau
#26 July 2016

#Based on McNair et al 2012; Seed Science Research

#My data is right-sensored.
#My data is interval, but can probably be analyzed as exact.
#Fully parametric time to event/failure-time/reliability analysis is most statistically sound for germination data.
#Time required for germination is a random variable. 

#Survivor function: probability that the germination time is greater than t
#Hazard function: how likely it is that a seed which has not germinated by time t will germinate shortly after t.

#N (number of seeds at risk of germination) should be adjusted for seed losses which occur before the end of the experiment.
#   This in cludes seeds damaged by handling, mold growth, lost seeds, etc.  Should not exceed 5%

#Random effects are referred to as frailty effects--these are most common as subgroups within treatment groups (e.g. reps, maybe regions?)

library("KMsurv") 
    ##lifetab() can be used to estimate the survivor fuction

library("survival")
    ##survfit() can be used to compute the Kalplan-Meier survivor function
    ##survdiff() can be used to to compare the survival patterns for two or more groups (with data assumed to be exact)

library("gtools")

#Set working directory and repositories
wd <- ("~/Google Drive/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

#Plotting nonparametric estimates of survivor fucntions for different treatement groups is a useful exploratory technique.

#read in germination data in pre-lifetab format
setwd(data.dir)



#############################

data1 <- read.csv("20160728_germdata_rep1_modern.csv", na.strings=c("","NA"))
data1[data1==Inf] <- 542.5

#data2 <- read.csv("20160728_germdata_rep2_modern.csv", na.strings=c("","NA"))
#data2[data2==Inf] <- 497.5

data <- rbind(data1,data2)

str(data)

#converts old data style to "modern" structure as described by McNair et al. 2012
data.mod <- data[rep(row.names(data), data$germinated), 1:11]

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
#################

df <- read.csv(paste(out.dir, "/alldata_mod_2016-08-01.csv", sep = ""), header = T) ##rep 1 and 2
str(df)
df$trt <- as.factor(df$trt)
df$run <- as.factor(df$run)

df$trt0 <- as.numeric(df$trt == 0)
df$trt10 <- as.numeric(df$trt == 10)
df$trt15 <- as.numeric(df$trt == 15)
df$trt20 <- as.numeric(df$trt == 20)

df$wcoast <- as.numeric(df$region == "wcoast")
df$ecoast <- as.numeric(df$region == "ecoast")
df$cv <- as.numeric(df$region == "central valleys")
df$yucatan <- as.numeric(df$region == "yucatan")
df$sierramadre <- as.numeric(df$region =="sierra madre")

df$run1 <- as.numeric(df$run == "1")
df$run2 <- as.numeric(df$run == "2")
df$run3 <- as.numeric(df$run == "3")
df$run4 <- as.numeric(df$run == "4")

################################

#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df, type = "kaplan-meier")
plot(km.fit, lwd = 2, col = "blue", conf.int = T, xlab = "time (hours)", ylab = "probability of not germinating", xlim = c(0,550))
print(summary(km.fit))

test.run <-survfit(Surv(end, status) ~ run, data = df, type = "kaplan-meier")
print(test.run)
plot(test.run, lty = c("solid", "longdash", "solid", "longdash"), col = c("grey", "grey", "black","black"),lwd = 1.75, xlab = "time (hours)", 
     ylab = "probability of not germinating", xlim = c(0,550), main = "Germination by run")

diff.run <- survdiff(Surv(end, status) ~ run, data = df)
print(diff.run)

test.trt <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier")
plot(test.trt, col = c("green", "blue", "yellow", "red"), lwd = 1.75, xlab = "time (hours)", 
     ylab = "probability of not germinating", xlim = c(0,550), main = "Germination by PEG concentration", )

test.region <- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier")
plot(test.region, col = c("black", "red", "green", "cyan", "blue"), lwd = 1.75, xlab = "time (hours)", 
     ylab = "probability of not germinating", xlim = c(0,550), main = "Germination by region")
legend(10, .3, c("central valleys", "ecoast", "sierra madre", "wcoast", "yucatan" ), 
       col =  c("black", "red", "green", "cyan", "blue"), pch = 19, bty = "n")

#define function mlogmlog() to calculate -log(-log(S(t))) checks proportional hazards (PH) assumption
mlogmlog <- function(y){-log(-log(y))}

#this model is reasonably robust, so only worry about clear departures from proportionality, 
#as indiecated by decisive crossing of the functions for two or more coariates in the diagnostic plot
plot(test.trt, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     col = c("green", "blue", "yellow", "red")) #looks good!

plot(test.run, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     lty = c("solid", "longdash", "solid", "longdash"), col = c("grey", "grey", "black","black")) #looks good!

plot(test.region, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     col = c("black", "red", "green", "cyan", "blue")) #looks good!

#create a cox model with two covariates
cox.region.trt <- coxph(Surv(end, status) ~ region + trt, data = df)
print(cox.region.trt)

cox.region.trt.regiontrt <- coxph(Surv(end,status) ~ region + trt + region:trt, data = df)
print(cox.region.trt.regiontrt)

cox.reg.run.trt <- coxph(Surv(end,status) ~ region + trt + region:trt + run, data = df)
print(cox.reg.run.trt)


######################
#ANALYSIS OF DATA FROM REP 1

df1 <- subset(df, rep ==1)

#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df, type = "kaplan-meier")
plot(km.fit, lwd = 2, col = "blue", conf.int = T, xlab = "time (hours)", ylab = "probability of not germinating", xlim = c(0,550))
print(summary(km.fit))

test.run.1 <-survfit(Surv(end, status) ~ run, data = df1, type = "kaplan-meier")
print(test.run.1)
plot(test.run.1, lty = c("solid", "longdash"), col = c("grey","black"),lwd = 1.75, xlab = "time (hours)", 
     ylab = "probability of not germinating", xlim = c(0,550), main = "Germination by run")

test.trt.1 <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier")
plot(test.trt.1, col = "black", lty = 1:4,
    lwd = 1.75, xlim = c(0,550), 
    #main = "Germination by PEG concentration"
    )
legend(5, .1, c("  0% PEG", "10% PEG", "15% PEG", "20% PEG"), lty = 1:4,
       col = "black", lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")


test.region.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier")
plot(test.region.1, col = c("black", "red", "green", "cyan", "blue"), lwd = 1.75, xlim = c(0,550), 
     #main = "Germination by region"
     )
legend(5, .15, c("Central Valleys", "E Coast", "Sierra Madre de la Sur", "W Coast", "Yucatan" ), 
       col =  c("black", "red", "green", "cyan", "blue"), lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

test.landrace.1 <- survfit(Surv(end, status) ~ landrace.name, data = df1, type = "kaplan-meier")
plot(test.landrace.1, col = rainbow(18), lwd = 1.75, xlim = c(0,550))
legend(5, .5, landraces, 
       col =  rainbow(18), lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")



#define function mlogmlog() to calculate -log(-log(S(t))) checks proportional hazards (PH) assumption
mlogmlog <- function(y){-log(-log(y))}

#this model is reasonably robust, so only worry about clear departures from proportionality, 
#as indiecated by decisive crossing of the functions for two or more coariates in the diagnostic plot
plot(test.trt.1, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     col = c("green", "blue", "yellow", "red")) #looks good!

plot(test.run.1, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     lty = c("solid", "longdash", "solid", "longdash"), col = c("grey", "grey", "black","black")) #looks good!

plot(test.region.1, mark.time = F, fun = mlogmlog, log = "x", xlab = "t (hours)", ylab = "-log(-log(S(t)))", 
     col = c("black", "red", "green", "cyan", "blue")) #looks good!


#create a cox model with two covariates
#cox.region.trt.1 <- coxph(Surv(end, status) ~ region + trt, data = df1)
#print(cox.region.trt.1)

cox.trt.1 <- coxph(Surv(end,status) ~ trt, data = df1)
print(summary(cox.trt.1) )

cox.trt.region.1 <- coxph(Surv(end,status) ~ trt + region, data = df1)
print(summary(cox.trt.region.1))

cox.trt.region.x.1 <- coxph(Surv(end,status) ~ trt + region + trt:region, data = df1)
print(summary(cox.trt.region.x.1))


##########################
# Building final model using forward selection

cox.0 <- coxph(Surv(end,status) ~ trt0, data = df1)
cox.0 #significant p =2.3e-12
cox.1 <- coxph(Surv(end,status) ~ trt10, data = df1)
cox.1 #significant p = 3.7e-10
cox.2 <- coxph(Surv(end,status) ~ trt15, data = df1)
cox.2 #NS
cox.3 <- coxph(Surv(end,status) ~ trt20, data = df1)
cox.3 #significant p = <2e-16
cox.4 <- coxph(Surv(end,status) ~ cv, data = df1)
cox.4 #significant p = 2.1e-10
cox.5 <- coxph(Surv(end,status) ~ sierramadre, data = df1)
cox.5 #NS
cox.6 <- coxph(Surv(end,status) ~ ecoast, data = df1)
cox.6 #significant p = <2e-16
cox.7 <- coxph(Surv(end, status) ~ wcoast, data = df1)
cox.7 #significant p = <2e-16
cox.8 <- coxph(Surv(end, status) ~ yucatan, data = df1)
cox.8 #significant p = .00014

cox.run1 <- coxph(Surv(end,status) ~ run1, data = df1)
cox.run1

cox.run2 <- coxph(Surv(end, status) ~ run2, data = df1)
cox.run2

#add covariates in order of significance 
cox.9 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + run2, data = df1)
cox.9 #all significant
cox.10 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + trt0 + run2, data = df1)
cox.10 #all significant
#cox.11 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + trt0 + cv, data = df1)
#cox.11 # ecoast NS
#cox.12 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + trt0 + trt10, data = df1)
#cox.12
#cox.13 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + trt0 +trt10 + yucatan, data = df1)
#cox.13

cox.11 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + trt0 + run2 + ecoast:trt20, data = df1)
cox.11 #NS
cox.12 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + trt0 + run2 + wcoast:trt20, data = df1)
cox.12 #NS
cox.13 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + trt0 + run2 + ecoast:trt0, data = df1)
cox.13 #NS
cox.14 <- coxph(Surv(end,status) ~ ecoast + wcoast + trt20 + trt0 + run2 + wcoast:trt0, data = df1)
cox.14 #NS

print(summary(cox.10))
