#OARDC Poster
#Germ Data Analysis (Rep 1)
# Vivian Bernau

# 12 February 2016

library(drc)
library(RCurl)
library(lme4)

wd <- ("~/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

data <- read.csv(paste(out.dir, "/data_2016-03-12.csv", sep = ""))

data <- data[with(data, order(type, trt, region, run, plate, start)),]

c1 <- subset(data, type == "C1")
c1.run1 <- subset(c1, run==1)
c1.run2 <- subset(c1, run ==2)
c2 <- subset(data, type == "C2")
c2.run1 <-subset(c2, run ==1)
c2.run2 <-subset(c2, run ==2)
e3 <- subset (data, type == "E3")
e3.run1 <- subset(e3, run ==1)
e3.run2 <- subset(e3, run ==2)
write.csv(c1, file = paste(out.dir, "/c1_", Sys.Date(), ".csv", sep = ""))
write.csv(c2, file = paste(out.dir, "/c2_", Sys.Date(), ".csv", sep = ""))
write.csv(e3, file = paste(out.dir, "/e3_", Sys.Date(), ".csv", sep = ""))
c1 <- read.csv(file = paste(out.dir, "/c1_", "2016-03-12", ".csv", sep = ""))
c2 <- read.csv(file = paste(out.dir, "/c2_", Sys.Date(), ".csv", sep = ""))
c3 <- read.csv(file = paste(out.dir, "/e3_", Sys.Date(), ".csv", sep = ""))

c1.run1 <- subset(c1, run==1)
c1.run2 <- subset(c1, run ==2)
c2.run1 <-subset(c2, run ==1)
c2.run2 <-subset(c2, run ==2)
e3.run1 <- subset(e3, run ==1)
e3.run2 <- subset(e3, run ==2)

c1.LL23 <-drm(germinated ~ start + end, region:factor(trt), data = c1, 
              fct = LL2.3(names = c("Slope", "Max", "T50")), type = "event", na.action = na.omit)
C1.summary <-summary(c1.LL23)

c2 <- c2[with(c2, order(type, trt, region, run, plate, start)),]
c2.LL23 <-drm(germinated ~ start + end, region:factor(trt), data = c2, 
              fct = LL2.3(names = c("Slope", "Max", "T50")), type = "event", na.action = na.omit)
C2.summary <- summary(c2.LL23)

e3 <- e3[with(e3, order(type, trt, region, run, plate, start)),]
e3.LL23 <- drm(germinated ~ start + end, region:factor(trt), data = e3, 
                fct = LL2.3(names = c("Slope", "Max", "T50")), type = "event", na.action = na.omit)
E3.summary <- summary(e3.LL23)

results <- merge.list(c(c1,c2,e3), list())

results.df <- read.csv(paste(out.dir,"/2016-03-13_results.csv",sep= ""))
max.df <- subset(results.df, Name = "Max")
t50.df <- subset(results.df, Name = "T50")
slope.df <- subset(results.df, Name = "Slope")

max.lmer <- lmer(Estimate ~ type + peg + 1|region + 1|run + 1|type:region + 
           type:peg + 1|type:run + peg|run + type:peg|run + type:peg|run, data = max.df)
