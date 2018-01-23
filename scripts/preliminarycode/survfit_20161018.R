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

df <- read.csv(paste(out.dir, "/alldata_mod_2016-10-12.csv", sep = ""), header = T) ##rep 1 and 2, run 5
str(df)
round(df$end, digits = 1)
df$trt <- as.factor(df$trt)
df$run <- as.factor(df$run)

#PEG = 0 ommitted
df$trt10 <- as.numeric(df$trt == 10)
df$trt15 <- as.numeric(df$trt == 15)
df$trt20 <- as.numeric(df$trt == 20)

#Johnny's ommitted
df$wcoast <- as.numeric(df$region == "wcoast")
df$ecoast <- as.numeric(df$region == "ecoast")
df$cv <- as.numeric(df$region == "central valleys")
df$yucatan <- as.numeric(df$region == "yucatan")
df$sierramadre <- as.numeric(df$region =="sierra madre")

#plantation ommitted
df$backyard <- as.numeric(df$cultivation == "Backyard")
df$forest<- as.numeric(df$cultivation == "Forest")
df$milpa<- as.numeric(df$cultivation == "Milpa")

#ANALYSIS OF DATA FROM Runs 1, 2, 5

#df1 <- df
#str(df1)
df1 <- subset(df, df$run ==1 | df$run ==2 | df$run ==5)
str(df1)

#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df1, type = "kaplan-meier")
plot(km.fit, lwd = 2, col = "blue", conf.int = T, xlab = "time (hours)", ylab = "probability of not germinating", xlim = c(0,550))

test.run.1 <-survfit(Surv(end, status) ~ run, data = df1, type = "kaplan-meier")
summary(test.run.1)
plot(test.run.1, col = rainbow(5),lwd = 1.75, xlab = "time (hours)", 
     ylab = "probability of not germinating", xlim = c(0,550), main = "Germination by Run")
legend(5, .2, c(1:5), col = rainbow(5), bty="n", cex = 1, lwd = 1.75)

test.trt.1 <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier")
plot(test.trt.1, col = "black", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), main = "Germination by PEG concentration")
legend(5, .1, c("  0% PEG", "10% PEG", "15% PEG", "20% PEG"), lty = 1:4,
       col = "black", lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

test.region.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier")
plot(test.region.1, col = c("black", "red", "grey", "green", "cyan", "blue"), 
     lwd = 1.40, xlim = c(0,550), 
     main = "Germination by region"
)
legend(5, .2, c("Central Valleys", "E Coast", "Johnny's", "Sierra Madre de la Sur", "W Coast", "Yucatan" ), 
       col =  c("black", "red", "grey", "green", "cyan", "blue"),lwd = 1.40, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

test.landrace.1 <- survfit(Surv(end, status) ~ landrace.name, data = df1, type = "kaplan-meier")
plot(test.landrace.1, col = rainbow(18), lwd = 1.75, xlim = c(0,550))
legend(5, .5, landraces, 
       col =  rainbow(18), lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

test.shelf.1 <- survfit(Surv(end, status) ~ shelf, data = df1, type = "kaplan-meier")
plot(test.shelf.1, col = rainbow(9), lwd = 1.4, xlim = c(0,550))
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
     col = c("black", "red", "grey","green", "cyan", "blue")) #hazard lines converge and cross

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


test.20.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier", subset = (trt ==20))
plot(test.20.1, col = c("black", "red", "grey", "green", "cyan", "blue"),
     lwd = 1.75, xlim = c(0,550)
    # main = "Germination by PEG concentration: 20% PEG"
    )
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")
legend(5, .2, c("Central Valleys", "E Coast", "Johnny's", "Sierra Madre del Sur", "W Coast", "Yucatan" ), 
       col =  c("black", "red", "green", "grey", "cyan", "blue"),lwd = 1.40, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")


test.0.1 <- survfit(Surv(end, status) ~ region, data = df1, type = "kaplan-meier", subset = (trt ==0))
plot(test.0.1, col = c("black", "red", "grey", "green", "cyan", "blue"),
     lwd = 1.75, xlim = c(0,550),
     #main = "Germination by PEG concentration: 0% PEG"
     )
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")

png(filename="allregions.tiff", width = 900, height = 1800, res = 120)
plot.new()
par(mfrow=c(5,1),mar = c(0,4,0,4)+0.1) 

test.cv <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="central valleys"))
plot(test.cv, col = "black", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), axes = F)
title(ylab = "probability of not germinating", cex.axis = 1.5)
axis(4)
axis(2)
box()

test.sm <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="sierra madre"))
plot(test.sm, col = "green", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), axes = F)
title(ylab = "probability of not germinating", cex.axis = 1.5)
axis(4)
axis(2)
box()

test.ecoast <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="ecoast"))
plot(test.ecoast, col = "red", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), axes = F)
title(ylab = "probability of not germinating", cex.axis = 1.5)
axis(4)
axis(2)
box()

test.wcoast <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="wcoast"))
plot(test.wcoast, col = "cyan", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), axes = F)
title(ylab = "probability of not germinating", cex.axis = 1.5)
axis(4)
axis(2)
box()

test.yucatan <- survfit(Surv(end, status) ~ trt, data = df1, type = "kaplan-meier", subset = (region =="yucatan"))
plot(test.yucatan, col = "blue", lty = 1:4,
     lwd = 1.75, xlim = c(0,550), axes = F)
legend(20, .7, c("  0% PEG", "10% PEG", "15% PEG", "20% PEG"), lty = 1:4,
       col = "black", lwd = 1.75, bty = "n", , cex = 1.5)
legend(20, .4, c("Central Valleys", "E Coast", "Sierra Madre del Sur", "W Coast", "Yucatan" ), 
       col =  c("black", "red", "green", "cyan", "blue"),lwd = 1.40, bty = "n", cex = 1.5)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating", cex.axis = 1.5)
axis(4)
axis(2)
box()
dev.off()

survreg <- survreg(Surv(end, status) ~ trt + region + trt:region + frailty(run), data = df1)
survreg

survreg_plate <- survreg(Surv(end, status) ~ trt + region + trt:region + frailty(run) + frailty(plate%in%run), data = df1)
survreg_plate

survreg_plate_line <- survreg(Surv(end, status) ~ trt + region + trt:region + frailty(run) + frailty(plate%in%run) + frailty(line%in%region), data = df1)
survreg_plate_line

test.cultivation <- survfit(Surv(end, status) ~ cultivation, data = df1, type = "kaplan-meier")
plot(test.cultivation, col = "black", lty = 1:4, lwd = 1.75, xlim = c(0,550))
legend(5, .5, c("backyard", "forest", "milpa", "plantation"), 
       col = "black", lty = 1:4, lwd = 1.75, bty = "n", cex = 1)
title(xlab = "time (hours)", line = 2, ylab = "probability of not germinating")


