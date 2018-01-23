# Graphs of Germ Data
# Vivian Bernau

# 13 February 2016

#rep = random
#block = random
#shelf(block) = random
#peg = fixed
#line = random
#pegXline = random
#error = random

#GENERAL MODEL STATEMENT:
# lmer(y ~ peg + (1 + peg | line) + (1 | block / shelf))

library(drc)

wd <- ("~/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(data.dir)
#list.files()
#data0 <- read.csv("germ_20160213c.csv")

#file.copy("~/RFiles/Chiles_OSU/Seed_Increase/seedincrease_sampleid.csv", data.dir)
#file.copy("~/RFiles/Chiles_OSU/Seed_Increase/2013_Collections.csv", data.dir)

#sampleID <- read.csv("seedincrease_sampleid.csv")
#colnames(sampleID) <- tolower(colnames(sampleID))

#data1 <- merge(data0, sampleID, by.x = "line", by.y = "line", all.x = T, ally.y = F)

#collections <- read.csv("2013_Collections.csv")
#colnames(collections) <- tolower(colnames(collections))
#collections <- collections[, (colnames(collections) %in% c("sample.id", "pop", "region", "population.type", 
                                                           "landrace.name", "cultivation", "main.use"))]
#data2 <- merge(data1, collections, by.x = "sampleid", by.y = "sample.id", all.x = T, all.y = F, row.names = "pedigree")

#data <- data2[with(data2, order(trt,plate, run, start)), ]
#write.csv(data, file = paste(out.dir, "/data.csv", sep = ""))

data <- read.csv(paste(out.dir, "/data.csv", sep = ""))
head(data)
tail(data)

data<- data[with(data, order(trt, plate, run, start)),]

data.LL3 <- drm(germinated~start+end, data = data, fct = LL.3(names=c("Slope","Max","T50")), 
                type = "event")
summary(data.LL3)  # showing a summmary of the model fit (including parameter estimates)
png(paste(out.dir, "/data_", Sys.Date(),".png", sep = ""), width=1080, height=720, res=120)
plot(data.LL3, log = "", ylim = c(0,1), xlim = c(0,500), ylab = "Proportion of germinated seed", xlab = "Time (h)",
    legendPos = c(20,.40))
dev.off()

#Curves fit by PEG level
data<- data[with(data, order(trt, plate, run, start)),]
peg.LL3 <- drm(germinated ~ start+end, trt, data=data, fct=LL.3(names=c("Slope","Max","T50")), 
               type="event")
summary(peg.LL3)
png(paste(out.dir, "/pegll3_", Sys.Date(),".png", sep = ""), width=1080, height=720, res=120)
plot(peg.LL3, log="", ylim=c(0,1), xlim=c(0,500), ylab="Proportion of germinated seed", xlab="Time (h)",
     legendPos=c(60,1), main = "Germination by PEG level")
dev.off()

#The ED5,ED50,and ED95 for the the germination
#ED(peg.LL3,c(5,50,95),interval="delta")

#Coefficient of variation for ED50 (not working)
#T50<-ED(peg.drc,c(50),display=FALSE)$EDdisplay[,2]*100/ED(peg.drc,c(50),display=FALSE)$EDdisplay[,1]
#T50

compParm(peg.LL3, "Slope")
compParm(peg.LL3, "T50")
compParm(peg.LL3, "Max")

library(plyr)

#Curves fitted by region of origin
data <- data[with(data, order(region, plate, run, start)),]
region.LL3 <- drm(germinated ~ start+end, region, data=data, fct=LL.3(names=c("Slope","Max","T50")), 
               type="event", na.action = na.omit)
summary(region.LL3)
png(paste(out.dir, "/regionll3_col_", Sys.Date(), ".png", sep = ""), width=1080, height=720, res=120)
plot(region.LL3, log="", ylim=c(0,1), xlim=c(0,500), col = c("black", "red", "grey", "green", "cyan", "blue"), 
     ylab="Proportion of germinated seed", xlab="Time (h)", main = "Germination by Region of Origin",
     legendPos=c(120,1))
dev.off()

compParm(region.LL3, "Slope")
compParm(region.LL3, "T50")
compParm(region.LL3, "Max")

#two-parameter log-logistic curves: region and PEG
#germinated ~ peg + (1 + peg | line) + (1 | run / shelf))
#germinated ~ peg + (1 + peg | region / line) + (1 | shelf) 
data <- data[with(data, order(region, trt, plate, run, start)),]
region.LL3 <- drm(germinated ~ start + end, (trt + (1 + trt|region/line)), 
                data = data, fct = LL.3(names=c("Slope","Max","T50")), type = "event", na.action=na.omit)
plot(region.LL3, ylim=c(0, 1), legendPos=c(2.5,1.5))  # plotting the fitted curves and the data
summary(germLL.2)  # showing the parameter estimates