#Germination Analysis
#Vivian Bernau
#11 Oct 2016
#18 Oct 2016 modified with input from Onofri et al. 2010 and Fox 2001

#Based on McNair et al 2012; Seed Science Research
#Updated based on Onofri et al. 2010, Weed Research

#My data is right-sensored.
#My data is interval, but can probably be analyzed as exact.
#Semi-parametric time to event/failure-time/reliability analysis is most statistically sound for germination data.
#Time required for germination is a random variable. 

#Survivor function: probability that the germination time is greater than t
#Hazard function: how likely it is that a seed which has not germinated by time t will germinate shortly after t.

#N (number of seeds at risk of germination) should be adjusted for seed losses which occur before the end of the experiment.
#   This in cludes seeds damaged by handling, mold growth, lost seeds, etc.  Should not exceed 5%

#Random effects are referred to as frailty effects--these are most common as subgroups within treatment groups (e.g. reps, maybe regions?)

library("gtools")
library("drc")

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

data3 <- read.csv("20170503_germdata_run6_modern.csv", na.strings = c("", "NA"))
data3 <- data3[1:10]
str(data3)

data4 <- read.csv("20170503_germdata_rep4_modern.csv", na.strings = c("", "NA"))
data4 <- data4[1:10]
str(data4)

data <- rbind(data1, data2, data3, data4)
str(data)
data.1 <- data[c(1:6, 8:10)]
data.na <- na.omit(data.1)
summary(data.na)
str(data.na)
#converts old data style to "modern" structure as described by McNair et al. 2012
###VERY IMPORTANT: Creates multiples of each "plate" for the number of seeds germinated at each time point
times = data.na$germinated
data.mod <- data.na[rep(row.names(data.na), times), 1:9]
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

df <- read.csv(paste(out.dir, "/alldata_mod_2017-05-04.csv", sep = ""), header = T) ##rep 1 and 2, run 5
str(df)
round(df$end, digits = 1)
df$run <- as.factor(df$run)

setwd(out.dir)

summary(df$region)
target <- c("central valleys", "ecoast", "sierra madre", "wcoast", "yucatan", "control")

library(gdata)
df$region <- reorder.factor(df$region, new.order=target)
summary(df$region)
head(df)
library(dplyr)
df1<-df %>% arrange(region)
head(df1)
summary(df$region)

library("survival")
#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df1, type = "kaplan-meier")
plot(km.fit, lwd = 2, col = "blue", conf.int = T, xlab = "time (hours)", ylab = "probability of not germinating", xlim = c(0,550))

test.run.1 <-survfit(Surv(end, status) ~ run, data = df1, type = "kaplan-meier")
#summary(test.run.1)
jpeg('test_run_1_20170504.jpg')
plot(test.run.1, col = rainbow(8),lwd = 1.75, xlab = "time (hours)", 
     ylab = "probability of not germinating", xlim = c(0,550), main = "Germination by Run")
legend(5, .4, c(1:8), col = rainbow(8), bty="n", cex = 1, lwd = 1.75)
dev.off()
#survdiff(Surv(end, status) ~ factor(run), data = df)

df1 <- subset(df, df$run ==1 | df$run ==2 |df$run ==3| df$run ==5) #why this subset?
str(df1)
