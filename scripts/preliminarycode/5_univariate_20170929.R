#GERMINATION DATA
#VIVIAN BERNAU
#29 September 2017

#UNIVARIATE ANALYSIS
#total viable
#total germ
#t50
#lag time
#germination rate (100/t50)
#germination uniformity (t75-t25)

#Set working directory and repositories
wd <- ("~/Google Drive/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(out.dir)

#read in germination data in pre-lifetab format
df <- read.csv(paste(out.dir, "/cleaned_2017-09-17.csv", sep = ""), header = T)
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
for(i in seq_along(plates)){
  x<- min(plates[[i]]$end[
    plates[[i]]$cumsum >= summary$germ50[[i]]])
  out50[[i]] <- x
}
summary$t50 <- t(as.data.frame(out50))
  
#calculate germation uniformity
summary$germ25 <- .25*summary$viable

out25 <- vector("list", n)
for(i in seq_along(plates)){
  x<- min(plates[[i]]$end[
    plates[[i]]$cumsum >= summary$germ25[[i]]])
  out25[[i]] <- x
}
summary$t25 <- t(as.data.frame(out25))

summary$germ75 <- .75*summary$viable
out75 <- vector("list", n)
for(i in seq_along(plates)){
  x<- min(plates[[i]]$end[
    plates[[i]]$cumsum >= summary$germ75[[i]]])
  out75[[i]] <- x
}
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

#distribution plots
plots <- subset(sum[,c("perc_germ", "delay", "uniform", "rategerm")])
str(plots)
plots$rategerm <- as.vector((plots$rategerm))
plots$uniform <- as.vector((plots$uniform))
str(plots)
plots <- na.omit(plots)

library(ggplot2)
library(reshape2)
ggplot(melt(plots),aes(x=value)) + geom_histogram() + facet_wrap(~variable, scales = "free_x") + stat_function(fun=dnorm)


library(lme4)
library(lmerTest)

options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                         check.nobs.vs.nlev = "warning",
                         check.nobs.vs.nRE = "ignore",
                         check.nlev.gtreq.5 = "warning",
                         check.nlev.gtr.1 = "warning"))

#####
#RATE OF GERMINATION
#####
rategerm.a <- glmer(rategerm ~ trt + (0+1|region) + (0+trt|region) + (0+1|landrace.name%in%region) + (0+trt|landrace.name%in%region)
                   + (1|rep) + (1|run%in%rep), data = sum)
summary(rategerm.a)
step(rategerm.a)
rategerm.a2 <- glmer(rategerm ~ perc_germ ~ trt + (0 + 1 | region) + (0 + trt | region) + (0 + 1 | landrace.name %in% region), data = sum)


#####
#UNIFORMITY
#####
uniform.a <- lmer(uniform ~ trt + (0+1|region) + (0+trt|region) + (0+1|landrace.name%in%region) + (0+trt|landrace.name%in%region)
                  + (1|rep) + (1|run%in%rep), data = sum)
summary(uniform.a)

#####
#DELAY
#####
delay.a <- lmer(delay ~ trt + (0+1|region) + (0+trt|region) + (0+1|landrace.name%in%region) + (0+trt|landrace.name%in%region) 
                + (1|rep) + (1|run%in%rep), data = sum)
summary(delay.a)
step(delay.a)
delay.a2 <- lmer(delay ~ trt + (0 + 1 | region) + (0 + trt | region) + (0 + 1 | landrace.name %in% region) + (1 | rep), data = sum)
summary(delay.a2)
lsmeansLT(delay.a2, test.effs = "trt", adjust = "tukey")
rand(delay.a2)

#percgerm is a percent value, so the data is right skewed. We account for this by using glmer with a guassian or logistic linking function.
perc_germ.a <- glmer(perc_germ ~ trt + (0+1|region) + (0+trt|region) + (0+1|landrace.name%in%region) + (0+trt|landrace.name%in%region)
          + (1|rep) + (1|run%in%rep), data = sum, family = gaussian(link = "identity"))
summary(perc_germ.a)
step(perc_germ.a)

perc_germ.a2 <- lmer(perc_germ ~ trt + (0+1|region) + (0+trt|region) + (0+1|landrace.name%in%region) + (0+trt|landrace.name%in%region)
                     + (1|rep) + (1|run%in%rep), data = sum)
summary(perc_germ.a)
step(perc_germ.a)

perc_germ.a2 <- lmer(perc_germ ~ trt + (0 + 1 | region) + (0 + trt | region) + (0 + 1 | landrace.name %in% region), data = sum, REML = F, 
           control = lmerControl(check.nlev.gtr.1 = "ignore"), contrasts = list(trt = "contr.SAS"))
summary(perc_germ.a2)
anova(perc_germ.a2)
rand(perc_germ.a2)
lsmeans(a2,spec="landrace.name")

perc_germ.b <- lmer(perc_germ ~ trt + (0+1|region) + (0+trt|cultivation) + (0 + 1|landrace.name%in%cultivation) + (0+trt|landrace.name%in%cultivation) + (1|rep) + (1|run%in%rep), 
          control=lmerControl(check.nlev.gtr.1 = "ignore"), data = sum, REML = F )
step(perc_germ.b)
perc_germ.b2 <- lmer(perc_germ ~)
anova(perc_germ.b2)
lsmeansLT(perc_germ.b2)
rand(perc_germ.b2)

