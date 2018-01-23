#OARDC Poster
#Germ Data Analysis (Rep 1)
# Vivian Bernau

# 12 February 2016

library(drc)
library(RCurl)
library(lme4)
library(gtools)

wd <- ("~/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(data.dir)

data1 <- read.csv("20160728_germdata_rep1.csv", na.strings=c("","NA"))
data2 <- read.csv("20160728_germdata_rep2.csv", na.strings=c("","NA"))
 data2$notes <- NA
identical(names(data1[[1]]), names(data2[[2]]) )
str(data1)
str(data2)
data <- smartbind(data1, data2)

sampleID <- read.csv("seedincrease_sampleid.csv")
colnames(sampleID) <- tolower(colnames(sampleID))

df1 <- merge(data, sampleID, by.x = "line", by.y = "line", all.x = T, ally.y = F)

collections <- read.csv("2013_Collections.csv")
colnames(collections) <- tolower(colnames(collections))
collections <- collections[, (colnames(collections) %in% c("sample.id", "pop", "region", "population.type", 
                                                           "landrace.name", "cultivation", "main.use"))]
df2 <- merge(df1, collections, by.x = "sampleid", by.y = "sample.id", all.x = T, all.y = F, row.names = "pedigree")
write.csv(df2, file = paste(out.dir, "/data_",Sys.Date(),".csv", sep = ""))

df <- read.csv(paste(out.dir, "/data_2016-07-28.csv", sep = "")) ##rep 1 and 2

c1 <- subset(df, type == "C1")
c2 <- subset(df, type == "C2")
e3 <- subset (df, type == "E3")

c1.LL23 <-drm(germinated~start+end, factor(trt), data = c1,
              fct = LL2.3(names = c("Slope", "Max", "T50")), type = "event", na.action = na.omit)

C1.summary <- summary(c1.LL23)

c2.LL23 <-drm(germinated ~ start + end, trt, data = c2, 
              fct = LL2.3(names = c("Slope", "Max", "T50")), type = "event", na.action = na.omit)
C2.summary <- summary(c2.LL23)

e3.LL23 <- drm(germinated ~ start + end, trt, data = e3, 
                fct = LL2.3(names = c("Slope", "Max", "T50")), type = "event", na.action = na.omit)
E3.summary <- summary(e3.LL23)

results.df <- read.csv(paste(out.dir,"/2016-07-28_results.csv",sep= ""))
max.df <- subset(results.df, Name = "Max")
t50.df <- subset(results.df, Name = "T50")
slope.df <- subset(results.df, Name = "Slope")

max.lmer <- lmer(Estimate ~ type + peg + 1|type/region, data = max.df)
max.lmer
summary(max.lmer)
