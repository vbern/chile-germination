#GERMINATION DATA
#VIVIAN BERNAU
#31 AUGUST 2017

#PLOT AND CLEAN DATA

#N (number of seeds at risk of germination) should be adjusted for seed losses which occur before the end of the experiment.
#   This in cludes seeds damaged by handling, mold growth, lost seeds, etc.  Should not exceed 5%

wd <- ("~/Google Drive/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(out.dir)

#read in germination data in pre-lifetab format
df <- read.csv(paste(out.dir, "/new_alldata_mod_2017-09-17.csv", sep = ""), header = T)
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

attach(df)
library(dplyr)
df$uniqueplate <- paste(run,"_",plate)
df$uniqueplate <-gsub(" ","",df$uniqueplate)
str(df$uniqueplate)

#check that total germinates !>10
summary <- aggregate(df$viable, by=list(df$uniqueplate), FUN=sum)
x <- strsplit(summary$Group.1,"_")
y <- as.data.frame(x)
z <- (t(y))
a <- as.data.frame(z)
summary$rep <- a$V1
summary$plate <- a$V2
rm("x", "y", "z", "a")
hist(summary$x)

#check that 1 + 0 !>10
x <- subset(df, df$status == 0)
y <- subset(df, df$status ==1)
x$status <- 1
z <- rbind(x, y)
summary2 <- aggregate(z$status, by=list(z$uniqueplate), FUN=sum)
hist(summary2$x)
a <- subset(summary2, x>10)

x <- subset.data.frame(summary, x >= 9)
hist(x$x)

clean <- df[df$uniqueplate %in% x$Group.1,]
str(clean)
write.csv(clean, file = paste(out.dir, "/cleaned_", Sys.Date(), ".csv", sep = ""))
