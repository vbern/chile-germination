#' ---
#' title: "Germination Analysis: univariate of full dataset"
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
install.packages("effects", repos="http://R-Forge.R-project.org", type="source") 
library(effects)
library(multcomp)


par(mfrow=c(2,2))

sum <- read.csv(file = paste(data.dir, "/summarydata_2018-01-23.csv", sep = ""), na.strings=c(""," ","NA"), header =T)
a <- subset(sum, !is.na(landrace.name)); summary(a$region)
b <- subset(a, !is.na(region)); summary(b$region)
c <- subset(b, !is.na(cultivation)); summary(c$region)
sum <- c

sum$trt <- as.character(sum$trt)

#distribution plots
plots <- subset(sum[,c("trt", "perc_notgerm", "delay", "uniform", "rategerm")])
str(plots)
plots$rategerm <- as.vector((plots$rategerm))
plots$uniform <- as.vector((plots$uniform))
str(plots)
plots <- na.omit(plots)

library(ggplot2)
library(reshape2)

x <- as.data.frame(melt(plots, id = "trt"))
ggplot(x, aes(x=value)) + geom_histogram() + facet_wrap(~variable, scales = "free_x")
ggplot(x, aes(x=variable, y = value, fill=trt)) + geom_boxplot() + facet_wrap( ~ variable, scale="free")

options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                         check.nobs.vs.nlev = "warning",
                         check.nobs.vs.nRE = "ignore",
                         check.nlev.gtreq.5 = "warning",
                         check.nlev.gtr.1 = "warning"))

#####
#Delay
#####

delay.a <- lmer(delay ~ 0 + trt*region + (1|region:landrace.name) +(trt|region:landrace.name) 
                            + (1|rep) + (1|rep:run) + (1|cultivation)
                            ,na.action = na.omit, data = sum, REML = F)
summary(delay.a)
anova(delay.a)
rand(delay.a)
options(max.print = 999999999)
detach("package:lmerTest", unload = T)
library(lsmeans)
library(multcompView)
lsm.delay.a <- lsmeans(delay.a, pairwise~trt, adjust="tukey")
cld(lsm.delay.a)

library(lmerTest)
plot(residuals(delay.a))
hist(residuals(delay.a))
plot(fitted(delay.a), residuals(delay.a))
qqnorm(resid(delay.a)) 
qqline(resid(delay.a)) 
ranef<- ranef(delay.a, condVar = T); ranef
dotplot(ranef)


delay.b <- lmer(delay ~ 0 + trt*cultivation + (1|landrace.name) + (trt|landrace.name) 
                 + (1|rep) + (1|rep:run) + (1|region)
                 ,na.action = na.omit, data = sum, REML = F)
summary(delay.b)
anova(delay.b)
rand(delay.b) 
delay.b.pairs <- lsmeansLT(delay.b, c("trt", "cultivation"))
pairs(delay.b.pairs)
lsm.delay.b <- lsmeans(delay.b, pairwise~trt:cultivation, adjust="tukey")
cld(lsm.delay.b)

plot(residuals(delay.b))
hist(residuals(delay.b))
plot(fitted(delay.b), residuals(delay.b))
qqnorm(resid(delay.b)) #check normal distribution of residuals
qqline(resid(delay.b)) #heavy tails
ranef<- ranef(delay.b, condVar = T) #extract conditional means
ranef
dotplot(ranef)

#####
#RATE OF GERMINATION
#####
summary(sum$rategerm)

rategerm.a <- lmer(100/t50 ~ 0 + trt*region + (1|region:landrace.name) +(trt|region:landrace.name) 
                                + (1|rep) + (1|rep:run) + (1|cultivation)
                                ,na.action = na.omit, data = sum, REML = F)
summary(rategerm.a)
library(lmerTest)
anova(rategerm.a)
rand(rategerm.a)
detach("package:lmerTest")
lsm.rategerm.a <- lsmeans(rategerm.a, pairwise~trt:region, adjust="tukey")
cld(lsm.rategerm.a)

par(mfrow=c(2,2))
plot(residuals(rategerm.a))
hist(residuals(rategerm.a))
plot(fitted(rategerm.a), residuals(rategerm.a))
qqnorm(resid(rategerm.a)) #check normal distribution of residuals
qqline(resid(rategerm.a)) #heavy tails present
ranef<- ranef(rategerm.a, condVar = T); ranef #extract conditional means
dotplot(ranef)

rategerm.b <- lmer(100/t50 ~ 0+trt*cultivation + (1|landrace.name) +(trt|landrace.name) 
                    + (1|rep) + (1|rep:run) + (1|region), 
                    data = sum, REML = F, na.action = na.omit)
summary(rategerm.b)
library(lmerTest)
anova(rategerm.b)
rand(rategerm.b)
detach("package:lmerTest")
lsm.rategerm.b<-lsmeans(rategerm.b, pairwise~trt:cultivation, adjust="tukey")
cld(lsm.rategerm.b)

plot(residuals(rategerm.b))
hist(residuals(rategerm.b))
plot(fitted(rategerm.b), residuals(rategerm.b))
qqnorm(resid(rategerm.b)) #check normal distribution of residuals
qqline(resid(rategerm.b))
ranef<- ranef(rategerm.b, condVar = T) #extract conditional means
ranef
dotplot(ranef)

#####
#UNIFORMITY
#####
summary(sum$uniform)

uniform.a <- lmer(uniform ~ 0 + trt*region + (1|region:landrace.name) +(trt|region:landrace.name) 
                   + (1|rep) + (1|rep:run)# + (1|cultivation)
                   ,na.action = na.omit, data = sum, REML = F)
summary(uniform.a)
anova(uniform.a)
rand(uniform.a)

par(mfrow=c(2,2))
plot(residuals(uniform.a))
hist(residuals(uniform.a))
plot(fitted(uniform.a), residuals(uniform.a))
qqnorm(resid(uniform.a)) #check normal distribution of residuals
qqline(resid(uniform.a)) #heavy tails present
ranef<- ranef(uniform.a, condVar = T); ranef #extract conditional means
dotplot(ranef)

uniform.b <- lmer(uniform ~ 0 + trt*cultivation + (1|landrace.name) + (trt|landrace.name) 
                  + (1|rep) + (1|rep:run) + (1|region), 
                  data = sum, REML = F, na.action = na.omit)
summary(uniform.b)
anova(uniform.b)
rand(uniform.b)
detach("package:lmerTest", unload = T)
lsmeans(uniform.b)

plot(residuals(uniform.b))
hist(residuals(uniform.b))
plot(fitted(uniform.b), residuals(uniform.b))
qqnorm(resid(uniform.b)) #check normal distribution of residuals
qqline(resid(uniform.b))
ranef<- ranef(uniform.b, condVar = T) #extract conditional means
ranef
dotplot(ranef)

#####
#TOTAL PRECENT GERM
#####
detach(lmerTest)
library(lsmeans)
summary(sum$perc_notgerm)

par(mfrow=c(1,1))
hist(sum$perc_germ)
sum$not10germ <- sum$perc_germ
  sum$not10germ <- as.numeric(sub(1, NA, sum$not10germ))
sum$did10germ <- sum$perc_germ
  sum$did10germ <- sum$did10germ==1
  sum$did10germ <- sum$did10germ*1
  
did10germ.a <- glmer(did10germ ~ 0 + trt*region + (1|region:landrace.name) + (trt|region:landrace.name)  
                     + (1|rep) + (1|rep:run) + (1|cultivation)
                     ,family = binomial(link = "logit"), na.action = na.omit, data = sum)
summary(did10germ.a)
anova(did10germ.a)
did10germ.a.landrace <- glmer(did10germ ~ 0 + trt*region + (1|rep) + (1|rep:run) + (1|cultivation)
                               ,family = binomial(link = "logit"), na.action = na.omit, data = sum)
anova(did10germ.a, did10germ.a.landrace)
did10germ.a.run <- glmer(did10germ ~ 0 + trt*region + (1|region:landrace.name) +(trt|region:landrace.name)  
                          + (1|rep) + (1|cultivation)
                          ,family = binomial(link = "logit"), na.action = na.omit, data = sum)
anova(did10germ.a.run, did10germ.a)
did10germ.a.cultivation <- glmer(did10germ ~ 0 + trt*region + (1|region:landrace.name)+ (trt|region:landrace.name)  
                                  + (1|rep) + (1|rep:run)
                                  ,family = binomial(link = "logit"), na.action = na.omit, data = sum)
anova(did10germ.a.cultivation, did10germ.a)
did10germ.a.rep <- glmer(did10germ ~ 0 + trt*region + (1|region:landrace.name)+ (trt|region:landrace.name)  
                                + (1|run) + (1|cultivation)
                                ,family = binomial(link = "logit"), na.action = na.omit, data = sum)
anova(did10germ.a.rep, did10germ.a)
did10germ.a.landracetrt <- glmer(did10germ ~ 0 + trt*region + (1|region:landrace.name) 
                     + (1|rep) + (1|rep:run) + (1|cultivation)
                     ,family = binomial(link = "logit"), na.action = na.omit, data = sum)
anova(did10germ.a.landracetrt, did10germ.a)

lsmeans(did10germ.a, c("trt", "region"), type = "response")
lsmeans(did10germ.a, "trt", type = "response")
lsmeans(did10germ.a, "region", type = "response")

par(mfrow=c(2,2))
plot(residuals(did10germ.a))
hist(residuals(did10germ.a))
plot(fitted(did10germ.a), residuals(did10germ.a))
qqnorm(resid(did10germ.a)) #check normal distribution of residuals
qqline(resid(did10germ.a)) #heavy tails present
ranef<- ranef(did10germ.a, condVar = T)
ranef #extract conditional means
dotplot(ranef)


did10germ.b <- glmer(did10germ ~ 0 + trt*cultivation + (trt|landrace.name) 
                    + (1|rep) + (1|rep:run) + (1|region)
                    ,family = binomial(link = "logit"), na.action = na.omit, data = sum)
summary(did10germ.b)
anova(did10germ.b)

did10germ.b.landrace <- glmer(did10germ ~ 0 + trt*cultivation + (1|rep) + (1|rep:run)
                              ,family = binomial(link = "logit"), data = sum)
anova(did10germ.b, did10germ.b.landrace)
did10germ.b.run <- glmer(did10germ ~ 0 + trt*cultivation + (1|landrace.name) +(trt|landrace.name)  
                         + (1|rep)
                         ,family = binomial(link = "logit"), data = sum)
anova(did10germ.b.run, did10germ.b)
did10germ.b.rep <- glmer(did10germ ~ 0 + trt*cultivation + (1|landrace.name)+ (trt|landrace.name)  
                         + (1|run)
                         ,family = binomial(link = "logit"),  data = sum)
anova(did10germ.b.rep, did10germ.b)
did10germ.b.landracetrt <- glmer(did10germ ~ 0 + trt*cultivation + (1|landrace.name) 
                                 + (1|rep) + (1|rep:run)
                                 ,family = binomial(link = "logit"), data = sum)
anova(did10germ.b.landracetrt, did10germ.b)

lsmeans(did10germ.b, c("trt", "cultivation"), type = "response")
lsmeans(did10germ.b, "trt", type = "response")
lsmeans(did10germ.b, "cultivation", type = "response")

plot(residuals(did10germ.b))
hist(residuals(did10germ.b))
plot(fitted(did10germ.b), residuals(did10germ.b))
qqnorm(resid(did10germ.b)) #check normal distribution of residuals
qqline(resid(did10germ.b))
ranef<- ranef(did10germ.b, condVar = T) #extract conditional means
ranef
dotplot(ranef)


not10germ.a <- lmer(not10germ ~ 0 + trt*region + (1|region:landrace.name) +(trt|region:landrace.name) 
                    + (1|rep) + (1|rep:run) + (1|cultivation), 
                    data = sum, REML = F, na.action = na.omit)
summary(not10germ.a)
anova(not10germ.a)
rand(not10germ.a)
lsm.not10.a<- lsmeans(not10germ.a, pairwise~trt:region, adjust="tukey")
cld(lsm.not10.a)

plot(residuals(not10germ.a))
hist(residuals(not10germ.a))
plot(fitted(not10germ.a), residuals(not10germ.a))
qqnorm(resid(not10germ.a)) #check normal distribution of residuals
qqline(resid(not10germ.a))
ranef<- ranef(not10germ.a, condVar = T) #extract conditional means
ranef
dotplot(ranef)


not10germ.b <- lmer(not10germ ~ 0 + trt*cultivation + (1|landrace.name) +(trt|landrace.name) 
                    + (1|rep) + (1|rep:run) #+ (1|region)
                    ,data = sum, REML = F, na.action = na.omit)

summary(not10germ.b)
anova(not10germ.b)
rand(not10germ.b)
lsm.not10.b<- lsmeans(not10germ.b, pairwise~trt:cultivation, adjust="tukey")
cld(lsm.not10.b)

plot(residuals(not10germ.b))
hist(residuals(not10germ.b))
plot(fitted(not10germ.b), residuals(not10germ.b))
qqnorm(resid(not10germ.b)) #check normal distribution of residuals
qqline(resid(not10germ.b))
ranef<- ranef(not10germ.b, condVar = T) #extract conditional means
ranef
dotplot(ranef)


#http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf