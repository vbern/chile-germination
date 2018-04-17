#' ---
#' title: "Germination Analysis: univariate subsets"
#' author: "Vivian Bernau"
#' date: "30 January 2018"
#' output: word_document
#' ---

#UNIVARIATE ANALYSIS
#total viable
#total germ
#t50
#lag time
#germination rate (100/t50)
#germination uniformity (t75-t25)
#####
#Analysis of Subsets
#####

wd <- ("~/Google Drive/RFiles/chile-germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

library(lme4)
library(lmerTest)
library(lattice)

par(mfrow=c(2,2))

sum <- read.csv(file = paste(data.dir, "/summarydata_2018-01-23.csv", sep = ""), na.strings=c(""," ","NA"), header =T)
a <- subset(sum, !is.na(landrace.name)); summary(a$region)
b <- subset(a, !is.na(region)); summary(b$region)
c <- subset(b, !is.na(cultivation)); summary(c$region)
sum <- c

tusta <- subset(sum, landrace.name == "Tusta")
str(tusta)

tusta.1a <- lmer(delay ~  0 + trt*region + trt*cultivation + trt*cultivation*region
                 + (1|rep) + (1|rep:run)
                 ,na.action = na.omit, data = tusta, REML = F)
summary(tusta.1a)
anova(tusta.1a) 
rand(tusta.1a)

plot(residuals(tusta.1a))
hist(residuals(tusta.1a))
plot(fitted(tusta.1a), residuals(tusta.1a))
qqnorm(resid(tusta.1a))
qqline(resid(tusta.1a))
#ranef<- ranef(tusta.1a, condVar = T) 
#ranef
#dotplot(ranef)

tusta.2a <- lmer(rategerm ~  0 + trt*region + trt*cultivation + trt*cultivation*region
                 + (1|rep) + (1|rep:run)
                 ,na.action = na.omit, data = tusta, REML = F)
summary(tusta.2a)
anova(tusta.2a) 
rand(tusta.2a) 

plot(residuals(tusta.2a))
hist(residuals(tusta.2a))
plot(fitted(tusta.2a), residuals(tusta.2a))
qqnorm(resid(tusta.2a))
qqline(resid(tusta.2a))
#ranef<- ranef(tusta.2a, condVar = T) 
#ranef
#dotplot(ranef)


tusta.3a <- lmer(uniform ~  0 + trt*region + trt*cultivation + trt*cultivation*region
                 + (1|rep) + (1|rep:run)
                 ,na.action = na.omit, data = tusta, REML = F)
summary(tusta.3a)
anova(tusta.3a) 
rand(tusta.3a) 

plot(residuals(tusta.3a))
hist(residuals(tusta.3a))
plot(fitted(tusta.3a), residuals(tusta.3a))
qqnorm(resid(tusta.3a))
qqline(resid(tusta.3a))
#ranef<- ranef(tusta.3a, condVar = T); ranef
#dotplot(ranef)

tusta$not10germ <- tusta$perc_germ
tusta$not10germ <- as.numeric(sub(1, NA, tusta$not10germ))
tusta$did10germ <- tusta$perc_germ
tusta$did10germ <- tusta$did10germ==1
tusta$did10germ <- tusta$did10germ*1

tusta.4a <- lmer(did10germ ~ 0 + trt*region + trt*cultivation + trt*cultivation*region
                 + (1|rep) + (1|rep:run)
                 ,na.action = na.omit, data = tusta, REML = F)
summary(tusta.4a)
anova(tusta.4a) 
rand(tusta.4a) 

plot(residuals(tusta.4a))
hist(residuals(tusta.4a))
plot(fitted(tusta.4a), residuals(tusta.4a))
qqnorm(resid(tusta.4a))

tusta.4b <- lmer(not10germ ~ 0 + trt*region + trt*cultivation + region*cultivation + trt*cultivation*region
                 + (1|rep) + (1|rep:run)
                 ,na.action = na.omit, data = tusta, REML = F)
summary(tusta.4b)
anova(tusta.4b) 
rand(tusta.4b) 

plot(residuals(tusta.4b))
hist(residuals(tusta.4b))
plot(fitted(tusta.4b), residuals(tusta.4b))
qqnorm(resid(tusta.4b))
qqline(resid(tusta.4a))
ranef<- ranef(tusta.4a, condVar = T) 
ranef
dotplot(ranef)