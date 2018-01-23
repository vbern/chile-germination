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

is.na(sum) <- sapply(sum, is.infinite)
sum <- sum[!(is.na(sum$landrace.name)) | !(is.na(sum$region)) | !(is.na(sum$cultivation)),]

sum0 <- subset(sum, sum$perc_germ!=0)

#write.csv(summary, file = paste(data.dir, "/summarydata_2018-01-23.csv", sep = ""))

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

plot.new()
par(mfrow=c(2,2))
X = na.omit(plots$delay)
h<-hist(X, breaks=20, col="GREY", xlab="delay", 
        main="Histogram of observations")
xfit<-seq(min(X),max(X),length=40)
yfit<-dnorm(xfit,mean=mean(X),sd=sd(X))
yfit <- yfit*diff(h$mids[1:2])*length(X)
lines(xfit, yfit, col="4", lwd=2)

X = na.omit(plots$uniform)
h<-hist(X, breaks=20, col="GREY", xlab="uniformity", 
        main="Histogram of observations")
xfit<-seq(min(X),max(X),length=40)
yfit<-dnorm(xfit,mean=mean(X),sd=sd(X))
yfit <- yfit*diff(h$mids[1:2])*length(X)
lines(xfit, yfit, col="4", lwd=2)

X = na.omit(plots$rategerm)
h<-hist(X, breaks=20, col="GREY", xlab="rate of germination", 
        main="Histogram of observations")
xfit<-seq(min(X),max(X),length=40)
yfit<-dnorm(xfit,mean=mean(X),sd=sd(X))
yfit <- yfit*diff(h$mids[1:2])*length(X)
lines(xfit, yfit, col="4", lwd=2)

X = na.omit(plots$perc_germ)
h<-hist(X, breaks=20, col="GREY", xlab="rate of germination", 
        main="Histogram of observations")
xfit<-seq(min(X),max(X),length=40)
yfit<-dnorm(xfit,mean=mean(X),sd=sd(X))
yfit <- yfit*diff(h$mids[1:2])*length(X)
lines(xfit, yfit, col="4", lwd=2)


#test for homogeneity of data (p>0.05 == data is homogenous)
bartlett.test(sum$delay~sum$trt)
plot(sum$trt, sum$delay, ylab = "delay")

bartlett.test(sum$uniform~sum$trt)
plot(sum$trt, sum$uniform, ylab = "uniform")

bartlett.test(sum$rategerm~sum$trt)
plot(sum$trt, sum$rategerm, ylab = "rategerm")

bartlett.test(sum$perc_germ~sum$trt)
plot(sum$trt, sum$perc_germ, ylab = "percnotgerm")

library(lme4)
library(lmerTest)
library(multcomp)
library(afex)
library(lattice)
library(pbkrtest)

options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                         check.nobs.vs.nlev = "warning",
                         check.nobs.vs.nRE = "ignore",
                         check.nlev.gtreq.5 = "warning",
                         check.nlev.gtr.1 = "warning"))

#a1 <- analysis by region, without including cultivation
#a2 <- analysis by region, cultivation included
#b1 <- analysis by cultivation, region not included
#b2 <- analysis by cultivation, region included

#####
#DELAY
#####
delay.a1 <- lmer(delay ~  trt + (1|region) + (trt|region) + (1|region:landrace.name) + (1|cultivation) + (0+trt|region:landrace.name)
                + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(delay.a1)
anova(delay.a1) #trt significant at p = 0.0047
anova(delay.a1, ddf = "Kenward-Roger")
rand(delay.a1) #trt:region:landrace significant at p<0.001, run within rep significnat at p = 0.008

plot(residuals(delay.a1))
hist(residuals(delay.a1))
plot(fitted(delay.a1), residuals(delay.a1))
qqnorm(resid(delay.a1)) #check normal distribution of residuals
qqline(resid(delay.a1)) #heavy tails
ranef<- ranef(delay.a1, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(delay.a1) > 200)

delay.a2 <- lmer(delay ~  trt + (1|region) + (trt|region) + (1|region:landrace.name) + (1|cultivation) + (0+trt|region:landrace.name)
                + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(delay.a2)
anova(delay.a2) #trt significant at p = 0.0047
anova(delay.a2, ddf = "Kenward-Roger")
rand(delay.a2) #trt:region:landrace significant at p<0.001, run within rep significnat at p = 0.008

plot(residuals(delay.a2))
hist(residuals(delay.a2))
plot(fitted(delay.a2), residuals(delay.a2))
qqnorm(resid(delay.a2)) #check normal distribution of residuals
qqline(resid(delay.a2)) #heavy tails
ranef<- ranef(delay.a2, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(delay.a2) > 200)


delay.b1 <- lmer(delay ~ trt + (1|cultivation) + (trt|cultivation) + (1|cultivation:landrace.name) +(1|region)+ (0+trt|cultivation:landrace.name)
                + (1|rep) + (1|rep:run), data = sum, REML = T, na.action = na.omit)
summary(delay.b1)
anova(delay.b1) #trt not significant
rand(delay.b1) #trt:cultivation:landrace significant at p < 0.001

plot.new()
plot(residuals(delay.b1))
hist(residuals(delay.b1))
plot(fitted(delay.b1), residuals(delay.b1))
qqnorm(resid(delay.b1)) #check normal distribution of residuals
qqline(resid(delay.b1)) #heavy tails
ranef<- ranef(delay.b1, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(delay.b1) > 300)
which(residuals(delay.b1) < -200)

delay.b2 <- lmer(delay ~ trt + (1|cultivation) + (trt|cultivation) + (1|cultivation:landrace.name) +(1|region)+ (0+trt|cultivation:landrace.name)
                + (1|rep) + (1|rep:run), data = sum, REML = T, na.action = na.omit)
summary(delay.b2)
anova(delay.b2) #trt not significant
rand(delay.b2) #trt:cultivation:landrace significant at p < 0.001

plot.new()
plot(residuals(delay.b2))
hist(residuals(delay.b2))
plot(fitted(delay.b2), residuals(delay.b2))
qqnorm(resid(delay.b2)) #check normal distribution of residuals
qqline(resid(delay.b2)) #heavy tails
ranef<- ranef(delay.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(delay.b2) > 300)
which(residuals(delay.b2) < -200)


#####
#RATE OF GERMINATION
#####
summary(sum$rategerm)

rategerm.a1 <- lmer(t50 ~ trt + (1|region) + (trt|region) + (1|region:landrace.name) + (0 + trt|region:landrace.name)
                   + (1|rep) + (1|rep:run), data = sum, REML = T, na.action = na.omit)
summary(rategerm.a1)
anova(rategerm.a1) #trt significant (p < 0.0099)
anova(rategerm.a1, ddf = "Kenward-Roger") #trt not significant, p =.5891
rand(rategerm.a1) #run nested in rep highly significant

par(mfrow=c(2,2))
plot(residuals(rategerm.a1))
hist(residuals(rategerm.a1))
plot(fitted(rategerm.a1), residuals(rategerm.a1))
qqnorm(resid(rategerm.a1)) #check normal distribution of residuals
qqline(resid(rategerm.a1)) #heavy tails present
ranef<- ranef(rategerm.a1, condVar = T); ranef #extract conditional means
dotplot(ranef)
#which(residuals(rategerm.a1) > 200)
#which(residuals(rategerm.a1) < -200)
#sum[c("2768","3537","8270","8330","8760","8900","8989","9068","9828"),]

rategerm.a2 <- lmer(t50 ~ trt + (1|region) + (trt|region) + (1|region:landrace.name)  + (1|cultivation) + (0 + trt|region:landrace.name)
                    + (1|rep) + (1|rep:run), data = sum, REML = T, na.action = na.omit)
summary(rategerm.a2)
anova(rategerm.a2) #trt significant (p < 0.0099)
anova(rategerm.a2, ddf = "Kenward-Roger") #trt not significant, p =.5891
rand(rategerm.a2) #run nested in rep highly significant

par(mfrow=c(2,2))
plot(residuals(rategerm.a2))
hist(residuals(rategerm.a2))
plot(fitted(rategerm.a2), residuals(rategerm.a2))
qqnorm(resid(rategerm.a2)) #check normal distribution of residuals
qqline(resid(rategerm.a2)) #heavy tails present
ranef<- ranef(rategerm.a2, condVar = T); ranef #extract conditional means
dotplot(ranef)
#which(residuals(rategerm.a2) > 200)
#which(residuals(rategerm.a2) < -200)
#sum[c("2768","3537","8270","8330","8760","8900","8989","9068","9828"),]


rategerm.b1 <- lmer(t50 ~ trt + (1|cultivation) + (trt|cultivation) + (1|cultivation:landrace.name) + (1|region) + 
                     (0+trt|cultivation:landrace.name)
                   + (1|rep) + (1|rep:run), data = sum, REML = T, na.action = na.omit)
summary(rategerm.b1)
anova(rategerm.b1, ddf = "Kenward-Roger") #trt highly significant (p < 0.0001)
rand(rategerm.b1) #run in rep highly significant, trt:cultivation significant at p = 0.09

plot(residuals(rategerm.b1))
hist(residuals(rategerm.b1))
plot(fitted(rategerm.b1), residuals(rategerm.b1))
qqnorm(resid(rategerm.b1)) #check normal distribution of residuals
qqline(resid(rategerm.b1))
ranef<- ranef(rategerm.b1, condVar = T) #extract conditional means
ranef
dotplot(ranef)
#which(residuals(rategerm.b1) > 200)
#which(residuals(rategerm.b1) < -180)

rategerm.b2 <- lmer(t50 ~ trt + (1|cultivation) + (trt|cultivation) + (1|cultivation:landrace.name) + (1|region) + 
                     (0+trt|cultivation:landrace.name)
                   + (1|rep) + (1|rep:run), data = sum, REML = T, na.action = na.omit)
summary(rategerm.b2)
anova(rategerm.b2, ddf = "Kenward-Roger") #trt highly significant (p < 0.0001)
rand(rategerm.b2) #run in rep highly significant, trt:cultivation significant at p = 0.09

plot(residuals(rategerm.b2))
hist(residuals(rategerm.b2))
plot(fitted(rategerm.b2), residuals(rategerm.b2))
qqnorm(resid(rategerm.b2)) #check normal distribution of residuals
qqline(resid(rategerm.b2))
ranef<- ranef(rategerm.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)
#which(residuals(rategerm.b) > 200)
#which(residuals(rategerm.b) < -180)

#####
#UNIFORMITY
#####
summary(sum$uniform)
uniform.a1 <- lmer(uniform ~ trt + (1|region) + (trt|region) + (1|region:landrace.name) +(1|cultivation) + (trt|region:landrace.name)
                  + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(uniform.a1)
anova(uniform.a1) #error in calculation
anova(uniform.a1, ddf = "Kenward-Roger") #trt not significant
rand(uniform.a1) #run within rep significant at p =.006, trt:region:landrace significant at p = 0.005

plot(residuals(uniform.a1)) #more positive than negative residuals
hist(residuals(uniform.a1))
plot(fitted(uniform.a1), residuals(uniform.a1))
qqnorm(resid(uniform.a1)) #check normal distribution of residuals
qqline(resid(uniform.a1)) #indicates right skew
ranef<- ranef(uniform.a1, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(uniform.a1) > 200)

uniform.a2 <- lmer(uniform ~ trt + (1|region) + (trt|region) + (1|region:landrace.name) +(1|cultivation) + (trt|region:landrace.name)
                  + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(uniform.a2)
anova(uniform.a2) #error in calculation
anova(uniform.a2, ddf = "Kenward-Roger") #trt not significant
rand(uniform.a2) #run within rep significant at p =.006, trt:region:landrace significant at p = 0.005

plot(residuals(uniform.a2)) #more positive than negative residuals
hist(residuals(uniform.a2))
plot(fitted(uniform.a2), residuals(uniform.a2))
qqnorm(resid(uniform.a2)) #check normal distribution of residuals
qqline(resid(uniform.a2)) #indicates right skew
ranef<- ranef(uniform.a2, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(uniform.a2) > 200)


uniform.b1 <- lmer(uniform ~ trt + (1|cultivation) + (trt|cultivation) + (1|cultivation:landrace.name) + (1|region) + (0+trt|cultivation:landrace.name)
                  + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(uniform.b1)
anova(uniform.b1) #error in calculation
anova(uniform.b1, ddf = "Kenward-Roger") #trt not significant
rand(uniform.b1) #cultivation:landrace significant at p = 0.003, run within rep significant at p = 0.02

plot(residuals(uniform.b1))
hist(residuals(uniform.b1))
plot(fitted(uniform.b1), residuals(uniform.b1))
qqnorm(resid(uniform.b1)) #check normal distribution of residuals
qqline(resid(uniform.b1)) #indicates right skew
ranef<- ranef(uniform.b1, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(uniform.b1) > 300)

uniform.b2 <- lmer(uniform ~ trt + (1|cultivation) + (trt|cultivation) + (1|cultivation:landrace.name) + (1|region) + (0+trt|cultivation:landrace.name)
                  + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(uniform.b2)
anova(uniform.b2) #error in calculation
anova(uniform.b2, ddf = "Kenward-Roger") #trt not significant
rand(uniform.b2) #cultivation:landrace significant at p = 0.003, run within rep significant at p = 0.02

plot(residuals(uniform.b2))
hist(residuals(uniform.b2))
plot(fitted(uniform.b2), residuals(uniform.b2))
qqnorm(resid(uniform.b2)) #check normal distribution of residuals
qqline(resid(uniform.b2)) #indicates right skew
ranef<- ranef(uniform.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(uniform.b2) > 300)

#####
#TOTAL PRECENT GERM
#####

#percgerm is a percent value, so the data is left skewed. We account for this by using glmer with a guassian or logistic linking function.
perc_germ.a1 <- glmer(perc_germ ~ trt + (1|region) + (trt|region) + (1|cultivation) + (1|region:landrace.name)  + (0+trt|region:landrace.name)
                     + (1|rep) + (1|rep:run), data = sum, family = gaussian(link = "identity"), REML = T)
summary(perc_germ.a1)
anova(perc_germ.a1) #trt not significant
rand(perc_germ.a1) #region significant at p= 0.1, trt:region:landrace significant at p < 0.0001, run within rep significant at p < 0.0001

plot(residuals(perc_germ.a1))
hist(residuals(perc_germ.a1))
plot(fitted(perc_germ.a1), residuals(perc_germ.a1))
qqnorm(resid(perc_germ.a1)) #check normal distribution of residuals
qqline(resid(perc_germ.a1)) #uneven tails, heavy
ranef<- ranef(perc_germ.a1, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(perc_germ.a1) < -.6)

perc_germ.a2 <- glmer(perc_germ ~ trt + (1|region) + (trt|region) + (1|cultivation) + (1|region:landrace.name)  + (0+trt|region:landrace.name)
                     + (1|rep) + (1|rep:run), data = sum, family = gaussian(link = "identity"), REML = T)
summary(perc_germ.a2)
anova(perc_germ.a2) #trt not significant
rand(perc_germ.a2) #region significant at p= 0.1, trt:region:landrace significant at p < 0.0001, run within rep significant at p < 0.0001

plot(residuals(perc_germ.a2))
hist(residuals(perc_germ.a2))
plot(fitted(perc_germ.a2), residuals(perc_germ.a2))
qqnorm(resid(perc_germ.a2)) #check normal distribution of residuals
qqline(resid(perc_germ.a2)) #uneven tails, heavy
ranef<- ranef(perc_germ.a2, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(perc_germ.a2) < -.6)

perc_germ.b1 <- glmer(perc_germ ~ trt + (1|cultivation) + (1|region) + (trt|cultivation) + (1|cultivation:landrace.name) + (0+trt|cultivation:landrace.name)
                     + (1|rep) + (1|rep:run), data = sum, family = gaussian(link = "identity"), REML = T)
summary(perc_germ.b1)
anova(perc_germ.b1) #trt not significant
rand(perc_germ.b1) #trt:cultivation significant at p= 0.03, trt:cultivation:landrace singificant at p<0.0001, run within rep significant at p<0.0001

plot(residuals(perc_germ.b1))
hist(residuals(perc_germ.b1))
plot(fitted(perc_germ.b1), residuals(perc_germ.b1))
qqnorm(resid(perc_germ.b1)) #check normal distribution of residuals
qqline(resid(perc_germ.b1))
ranef<- ranef(perc_germ.b1, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(perc_germ.b1) < -.6)

perc_germ.b2 <- glmer(perc_germ ~ trt + (1|cultivation) + (1|region) + (trt|cultivation) + (1|cultivation:landrace.name) + (0+trt|cultivation:landrace.name)
                     + (1|rep) + (1|rep:run), data = sum, family = gaussian(link = "identity"), REML = T)
summary(perc_germ.b2)
anova(perc_germ.b2) #trt not significant
rand(perc_germ.b2) #trt:cultivation significant at p= 0.03, trt:cultivation:landrace singificant at p<0.0001, run within rep significant at p<0.0001

plot(residuals(perc_germ.b2))
hist(residuals(perc_germ.b2))
plot(fitted(perc_germ.b2), residuals(perc_germ.b2))
qqnorm(resid(perc_germ.b2)) #check normal distribution of residuals
qqline(resid(perc_germ.b2))
ranef<- ranef(perc_germ.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)
which(residuals(perc_germ.b2) < -.6)


#####
#ALL-RANDOM MODELS FOR CALCULATING HERITABILITY
#####
#VG = VG/(VG + VGT/t + Ve/r)
#Where VG is equal to your genetic variance, 
#VGT is equal to your genotype x treatment variance, 
#t is your number of treatments, Ve is your residual variance,
#and r is your number of reps.
###
delay.h <- lmer(delay ~  (1|trt) + (1|landrace.name) + (1|trt:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(delay.h)
#3587.4/(3487.4 + 917.9/4 +4427.5/4) #0.74369
3628.6/(3628.6 + 920.44/4 + 4427.06/4) #.7308

#delay.h2 <- lmer(delay ~  (1|trt) + (1|sampleid) + (1|trt:sampleid) + (1|rep) + (1|rep:run), data = sum, REML = T)
#summary(delay.h2)
#2442.55/(2442.55 + 871.8/4 +3767.97/4) #0.6780

uniformity.h <- lmer(uniform ~  (1|trt) + (1|landrace.name) + (1|trt:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(uniformity.h)
#832.25/(832.25 + 287.26/4 +3102.8/4) #0.495
877.63/(877.63+285.83/4+3101.85/4) #.5089

#uniformity.h2 <- lmer(uniform ~  (1|trt) + (1|sampleid) + (1|trt:sampleid) + (1|rep) + (1|rep:run), data = sum, REML =T)
#summary(uniformity.h2)
#754.38/(754.38 + 487.51/4 +2696.66/4) #0.4865

rate.h <- lmer(rategerm ~  (1|trt) + (1|landrace.name) + (1|trt:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(rate.h)
#3.885e-06/(3.885e-06 + 1.575e-07/4 + 3.166e-06/4) #0.8238
3.956e-06/(3.956e-06+1.583e-07/4+3.166e-06/4) #.8264

#rate.h2 <- lmer(rategerm ~  (1|trt) + (1|sampleid) + (1|trt:sampleid) + (1|rep) + (1|rep:run), data = sum, REML = T)
#summary(rate.h2)
#2.939e-06/(2.939e-06+2.869e-07/4+2.262e-06/4) #0.8218

percgerm.h <- lmer(perc_germ ~  (1|trt) + (1|landrace.name) + (1|trt:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(percgerm.h)
#1.636e-02/(1.636e-02+1.291e-02/4+3.268e-02/4) #0.5894
0.01649/(0.01649+0.01292/4+0.03268/4) #.5912

#percgerm.h2 <- lmer(perc_germ ~  (1|trt) + (1|sampleid) + (1|trt:sampleid) + (1|rep) + (1|rep:run), data = sum, REML = T)
#summary(percgerm.h2)
#9.173e-03/(9.173e-03+1.472e-02/4+2.614e-02/4) #0.4731

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
