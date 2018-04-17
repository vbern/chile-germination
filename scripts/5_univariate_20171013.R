#' ---
#' title: "Germination Analysis: univariate"
#' author: "Vivian Bernau"
#' date: "29 September 2017"
#' output: word_document
#' ---

#UNIVARIATE ANALYSIS
#total viable
#total germ
#t50
#lag time
#germination rate (100/t50)
#germination uniformity (t75-t25)

#Set working directory and repositories
wd <- ("~/Google Drive/RFiles/chile-germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(out.dir)

#read in germination data in pre-lifetab format
#df <- read.csv(paste(out.dir, "/cleaned_2017-09-17.csv", sep = ""), header = T)
df <- read.csv(paste(out.dir, "/cleaned8_2017-10-13.csv", sep = ""), header = T)
str(df)

#extract values for viable seed and germinated seed per plate
viable <- aggregate(df$viable, by=list(df$uniqueplate), FUN=sum)
colnames(viable) <- c("uniqueplate", "viable")
germ <- aggregate(df$status, by = list(df$uniqueplate), FUN = sum)
colnames(germ) <- c("uniqueplate", "germ")
summary <- merge(viable, germ, by.x = "uniqueplate")
rm(germ, viable)

#calculate percent germ and percent viable
summary$perc_germ <- summary$germ/summary$viable
summary$perc_notgerm <- 1-summary$perc_germ
summary$perc_viable <- summary$viable/10

#calculate delay to first germ (first data point per plate)
delay <- aggregate(df$end ~ df$uniqueplate, FUN = min)
colnames(delay) <- c("uniqueplate", "delay")
summary$delay <- delay$delay

#calculate germination rate
#at which time point is the sum of end > germ50?
summary$germ50 <- .5*summary$viable

df <- df[order(df$end, (df$plate), (df$run)),]

plates <- split(df, df$uniqueplate)

for(i in seq_along(plates)){
    x <- cumsum(plates[[i]][,"status"])
    plates[[i]]$cumsum <- x
}
plates[[2]]$cumsum #confirm loop was successful

min(plates[[1]]$end[plates[[1]]$cumsum >= summary$germ50[1]])

n <- nrow(summary)
out50 <- vector("list", n)
suppressWarnings(
  for(i in seq_along(plates)){
  x<- min(plates[[i]]$end[
    plates[[i]]$cumsum >= summary$germ50[[i]]])
  out50[[i]] <- x
})
summary$t50 <- t(as.data.frame(out50))
  
#calculate germation uniformity
summary$germ25 <- .25*summary$viable

out25 <- vector("list", n)
suppressWarnings(
  for(i in seq_along(plates)){
  x<- min(plates[[i]]$end[
    plates[[i]]$cumsum >= summary$germ25[[i]]])
  out25[[i]] <- x
})
summary$t25 <- t(as.data.frame(out25))

summary$germ75 <- .75*summary$viable
out75 <- vector("list", n)
suppressWarnings(
  for(i in seq_along(plates)){
  x<- min(plates[[i]]$end[
    plates[[i]]$cumsum >= summary$germ75[[i]]])
  out75[[i]] <- x
})
summary$t75 <- t(as.data.frame(out75))

#calculuate germ rate and uniformity
summary$rategerm <- 1/summary$t50
summary$uniform <- summary$t75-summary$t25

#merge summary statistics with descriptive datasets
sum <- merge(summary, df, by = "uniqueplate", all.x = T)
sum = sum[!duplicated(sum$uniqueplate),]

sum$run <- as.factor(sum$run)
sum$rep <- as.factor(sum$rep)
sum$trt <- as.factor(sum$trt)

is.na(sum) <- sapply(sum, is.infinite)
sum <- sum[!(is.na(sum$landrace.name)) | !(is.na(sum$region)) | !(is.na(sum$cultivation)),]

sum0 <- subset(sum, sum$perc_germ!=0)
remove <- c("X.1", "X.2", "X", "viable.y", "number", "status")
sum1 <- sum0[ , !(names(sum0) %in% remove)]
write.csv(sum1, file = paste(data.dir, "/summarydata_2018-01-23.csv", sep = ""))

#distribution plots
plots <- subset(sum[,c("trt", "perc_notgerm", "delay", "uniform", "rategerm")])
str(plots)
plots$rategerm <- as.vector((plots$rategerm))
plots$uniform <- as.vector((plots$uniform))
str(plots)
plots <- na.omit(plots)

suppressMessages(
library(ggplot2,reshape2))
x <- melt(plots, id = "trt")
ggplot(x, aes(x=value)) + geom_histogram() + facet_wrap(~variable, scales = "free_x")
ggplot(x, aes(x=variable, y = value)) + geom_boxplot(aes(fill=trt)) + facet_wrap( ~ variable, scales="free")

suppressMessages(
library(lme4,lmerTest,multcomp,afex,lattice,pbkrtest)
)

options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                         check.nobs.vs.nlev = "warning",
                         check.nobs.vs.nRE = "ignore",
                         check.nlev.gtreq.5 = "warning",
                         check.nlev.gtr.1 = "warning"))

#####
#CORRELATION OF UNIVARIATE VARBIABLES
#####
library(Hmisc)
library(corrplot)

rcorr <- rcorr(as.matrix(sum[c("perc_germ", "rategerm", "uniform", "delay")]))

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flat_corr <- flattenCorrMatrix(rcorr$r, rcorr$P)

corrplot(rcorr$r, type="upper", order="hclust", hclust.method = "centroid",
         p.mat = rcorr$P, sig.level = 0.05, insig = "blank", 
         tl.col = "black")
title(main ="Pearson Correlation Matrix", sub = "size and color of circles related to
      the correlation strength, 
      correlations where p<0.05 left blank")



