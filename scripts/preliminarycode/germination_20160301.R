# Graphs of Germ Data
# Vivian Bernau

# 13 February 2016

library(drc)

wd <- ("~/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(data.dir)
list.files()
data0 <- read.csv("20160312_germdata_rep1.csv")

#file.copy("~/RFiles/Chiles_OSU/Seed_Increase/seedincrease_sampleid.csv", data.dir)
#file.copy("~/RFiles/Chiles_OSU/Seed_Increase/2013_Collections.csv", data.dir)

sampleID <- read.csv("seedincrease_sampleid.csv")
colnames(sampleID) <- tolower(colnames(sampleID))

data1 <- merge(data0, sampleID, by.x = "line", by.y = "line", all.x = T, ally.y = F)

collections <- read.csv("2013_Collections.csv")
colnames(collections) <- tolower(colnames(collections))
collections <- collections[, (colnames(collections) %in% c("sample.id", "pop", "region", "population.type", 
                                                           "landrace.name", "cultivation", "main.use"))]
data2 <- merge(data1, collections, by.x = "sampleid", by.y = "sample.id", all.x = T, all.y = F, row.names = "pedigree")

data <- data2[with(data2, order(trt,plate, run, start.time)), ]
write.csv(data, file = paste(out.dir, "/data_",Sys.Date(),".csv", sep = ""))

rep1 <- read.csv(paste(out.dir, "/data_2016-03-12.csv", sep = "")) ##rep 1
#rep2 <- read.csv(paste(out.dir, "/data_2016-XX-XX.csv", sep = "")) ##rep 2

data <- rbind(rep1 #, rep2
)

#write.csv(data, file = paste(out.dir, "/alldata_", Sys.Date(), ".csv", sep = ""))

head(data)
tail(data)

#Curves fit by PEG level
data<- data[with(data, order(trt, plate, run, start)),]
peg.LL3 <- drm(germinated ~ start+end, trt, data=data, fct=LL.3(names=c("Slope","Max","T50")), 
               type="event", na.action = na.omit)
summary(peg.LL3)
png(paste(out.dir, "/pegll3_", Sys.Date(),".png", sep = ""), width=1080, height=720, res=120)
plot(peg.LL3, log="", ylim=c(0,1), xlim=c(0,500), ylab="Proportion of germinated seed", xlab="Time (h)",
     legendPos=c(60,1), main = "Germination by PEG level (LL.3)")
dev.off()

peg.LL23 <- drm(germinated ~ start+end, trt, data=data, fct=LL2.3(names=c("Slope","Max","T50")), 
               type="event", na.action = na.omit)
summary(peg.LL23)
png(paste(out.dir, "/pegll23_", Sys.Date(),".png", sep = ""), width=1080, height=720, res=120)
plot(peg.LL23, log="", ylim=c(0,1), xlim=c(0,500), ylab="Proportion of germinated seed", xlab="Time (h)",
     legendPos=c(60,1), main = "Germination by PEG level (LL2.3)")
dev.off()

compParm(peg.LL3, "Slope")
compParm(peg.LL3, "T50")
compParm(peg.LL3, "Max")

library(plyr)

#Curves fitted by region of origin
data <- data[with(data, order(region, run, plate, start)),]
region.LL3 <- drm(germinated ~ start+end, region, data=data, fct=LL.3(names=c("Slope","Max","T50")), 
               type="event", na.action = na.omit)
summary(region.LL3)
png(paste(out.dir, "/regionll3_col_", Sys.Date(), ".png", sep = ""), width=1260, height=720, res=120)
plot(region.LL3, log="", ylim=c(0,1), xlim=c(0,500), col = c("black", "red", "grey", "green", "cyan", "blue"), 
     ylab="Proportion of germinated seed", xlab="Time (h)", main = "Germination by Region of Origin (LL.3)",
     legendPos=c(120,1))
dev.off()

region.LL23 <- drm(germinated ~ start+end, region, data=data, fct=LL2.3(names=c("Slope","Max","T50")), 
                  type="event", na.action = na.omit)
summary(region.LL23)
png(paste(out.dir, "/regionll23_col_", Sys.Date(), ".png", sep = ""), width=1080, height=720, res=120)
plot(region.LL23, log="", ylim=c(0,1), xlim=c(0,500), col = c("black", "red", "grey", "green", "cyan", "blue"), 
     ylab="Proportion of germinated seed", xlab="Time (h)", main = "Germination by Region of Origin (LL2.3)",
     legendPos=c(120,1))
dev.off()


# Curves fit by type (control v. experimental)
data <- data[with(data, order(type, run, plate, start)),]
type.LL3 <- drm(germinated ~ start+end, type, data=data, fct=LL.3(names=c("Slope","Max","T50")), 
                  type="event", na.action = na.omit)
summary(type.LL3)

type.LL23 <- drm(germinated ~ start+end, type, data=data, fct=LL2.3(names=c("Slope","Max","T50")), 
                   type="event", na.action = na.omit)
summary(type.LL23)

#two-parameter log-logistic curves: region and PEG 
data <- data[with(data, order(region, trt, run, plate, start)),]
data.LL3 <- drm(germinated ~ region + trt, 
                data = data, fct = LL.3(names=c("Slope","Max","T50")), 
                type = "event", na.action=na.omit, ursa(fixed = trt))
plot(region.LL3, ylim=c(0, 1), legendPos=c(2.5,1.5))  # plotting the fitted curves and the data
summary(germLL.3)  # showing the parameter estimates