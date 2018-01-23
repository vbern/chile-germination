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

#Set working directory and repositories
wd <- ("~/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

#Plotting nonparametric estimates of survivor fucntions for different treatement groups is a useful exploratory technique.

#read in germination data in pre-lifetab format
setwd(data.dir)

data1 <- read.csv("20160728_germdata_rep1_modern.csv", na.strings=c("","NA"))
data2 <- read.csv("20160728_germdata_rep2_modern.csv", na.strings=c("","NA"))

data <- smartbind(data1,data2)

str(data)

data.mod <- data[rep(row.names(data), data$germinated), 1:11]

write.csv(data.mod, file = paste(out.dir, "/data_mod_",Sys.Date(),".csv", sep = ""))

df <- read.csv(paste(out.dir, "/data_mod_2016-07-28.csv", sep = ""), header = T) ##rep 1 and 2
str(df)

km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df, type = "kaplan-meier")
plot(km.fit, lwd = 2, col = "blue", conf.int = T, xlab = "time (hours)", ylab = "probability of not germinating")
print(summary(km.fit))
