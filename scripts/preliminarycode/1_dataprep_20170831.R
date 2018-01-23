#GERMINATION DATA
#VIVIAN BERNAU
#31 AUGUST 2017

wd <- ("~/Google Drive/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(data.dir)

#load data
data1 <- read.csv("rep1.csv", na.strings=c("","NA"))
data1 <- data1[1:10]
str(data1)

data2 <- read.csv("rep2.csv", na.strings=c("","NA"))
data2 <- data2[1:10]
str(data2)

data3 <- read.csv("rep3.csv", na.strings=c("","NA"))
data3 <- data3[1:10]
str(data3)

data4 <- read.csv("rep4.csv", na.strings=c("","NA"))
data4 <- data4[1:10]
str(data4)

#merge datasets
data <- rbind(data1, data2, data3, data4)
str(data)
df <- data[c(1:6, 8:10)]
df.na <- na.omit(df)
summary(df.na)
str(df.na)

#convert old data style to "modern" structure as described by McNair et al. 2012
###VERY IMPORTANT: Creates multiples of each "plate" for the number of seeds germinated at each time point
times = df.na$number
df.mod <- df.na[rep(row.names(df.na), times), 1:9]
str(df.mod)
write.csv(df.mod, file = paste(out.dir, "/new_data_mod_",Sys.Date(),".csv", sep = ""))

#merge germination data with descriptive datasets
df.mod <- read.csv(file = paste(out.dir, "/new_data_mod_",Sys.Date(),".csv", sep = ""))
sampleID <- read.csv("seedincrease_sampleid.csv")
colnames(sampleID) <- tolower(colnames(sampleID))

df.sample <- merge(df.mod, sampleID, by.x = "line", by.y = "line", all.x = T, ally.y = F)

collections <- read.csv("2013_Collections.csv")
colnames(collections) <- tolower(colnames(collections))
collections <- collections[, (colnames(collections) %in% c("sample.id", "pop", "region", "population.type", 
                                                           "landrace.name", "cultivation", "main.use"))]
df.collections <- merge(df.sample, collections, by.x = "sampleid", by.y = "sample.id", all.x = T, all.y = F, row.names = "pedigree")
df.collections$loc_id <- NULL

head(df.collections)
tail(df.collections)
str(df.collections)

write.csv(df.collections, file = paste(out.dir, "/new_alldata_mod_",Sys.Date(),".csv", sep = ""))

#LAST STEP: must go through and edit all of the Johnny's records so that region appears

