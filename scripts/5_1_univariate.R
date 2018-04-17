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
library(effects)
library(pbkrtest)

par(mfrow=c(2,2))

sum <- read.csv(file = paste(data.dir, "/summarydata_2018-01-23.csv", sep = ""), na.strings=c(""," ","NA"), header =T)
a <- subset(sum, !is.na(landrace.name)); summary(a$region)
b <- subset(a, !is.na(region)); summary(b$region)
c <- subset(b, !is.na(cultivation)); summary(c$region)
sum <- c

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

#a1 <- analysis by region, without cultivation
#a2 <- analysis by region, cultivation included
#b1 <- analysis by cultivation, region not included
#b2 <- analysis by cultivation, region included

#####
#Delay
#####

# delay.a1 <- lmer(delay ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = F)
# summary(delay.a1)
# rand(delay.a1)
# 
# plot(residuals(delay.a1))
# hist(residuals(delay.a1))
# plot(fitted(delay.a1), residuals(delay.a1))
# qqnorm(resid(delay.a1))
# qqline(resid(delay.a1))
# ranef<- ranef(delay.a1, condVar = T)
# ranef
# dotplot(ranef)


delay.a2 <- lmer(delay ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run) + (1|cultivation),
                 na.action = na.omit, data = sum, REML = F)
summary(delay.a2)
rand(delay.a2) 

plot(residuals(delay.a2))
hist(residuals(delay.a2))
plot(fitted(delay.a2), residuals(delay.a2))
qqnorm(resid(delay.a2)) 
qqline(resid(delay.a2)) 
ranef<- ranef(delay.a2, condVar = T) 
ranef
dotplot(ranef)

# delay.b1 <- lmer(delay ~ trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
#                  data = sum, REML = F, na.action = na.omit)
# summary(delay.b1)
# rand(delay.b1) 
# 
# plot(residuals(delay.b1))
# hist(residuals(delay.b1))
# plot(fitted(delay.b1), residuals(delay.b1))
# qqnorm(resid(delay.b1)) 
# qqline(resid(delay.b1)) 
# ranef<- ranef(delay.b1, condVar = T) 
# ranef
# dotplot(ranef)
# which(residuals(delay.b1) > 300)


delay.b2 <- lmer(delay ~  trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run) + (1|region), 
                 data = sum, REML = F, na.action = na.omit)
summary(delay.b2)
rand(delay.b2) 

plot(residuals(delay.b2))
hist(residuals(delay.b2))
plot(fitted(delay.b2), residuals(delay.b2))
qqnorm(resid(delay.b2)) #check normal distribution of residuals
qqline(resid(delay.b2)) #heavy tails
ranef<- ranef(delay.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)


delay.c <- lmer(delay ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + 
                   (1+trt|cultivation) +(1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
                na.action = na.omit, data = sum, REML = F)

summary(delay.c)
rand(delay.c)

plot(residuals(delay.c))
hist(residuals(delay.c))
plot(fitted(delay.c), residuals(delay.c))
qqnorm(resid(delay.c)) #check normal distribution of residuals
qqline(resid(delay.c)) #heavy tails
ranef<- ranef(delay.c, condVar = T) #extract conditional means
ranef
dotplot(ranef)

#####
#RATE OF GERMINATION
#####
summary(sum$rategerm)

rategerm.a2 <- lmer(rategerm ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run) + (1|cultivation), 
                    na.action = na.omit, data = sum, REML = F)
summary(rategerm.a2)
rand(rategerm.a2)

par(mfrow=c(2,2))
plot(residuals(rategerm.a2))
hist(residuals(rategerm.a2))
plot(fitted(rategerm.a2), residuals(rategerm.a2))
qqnorm(resid(rategerm.a2)) #check normal distribution of residuals
qqline(resid(rategerm.a2)) #heavy tails present
ranef<- ranef(rategerm.a2, condVar = T); ranef #extract conditional means
dotplot(ranef)



# rategerm.b1 <- lmer(rategerm ~ trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
#                     data = sum, REML = F, na.action = na.omit)
# summary(rategerm.b1)
# rand(rategerm.b1)
# 
# plot(residuals(rategerm.b1))
# hist(residuals(rategerm.b1))
# plot(fitted(rategerm.b1), residuals(rategerm.b1))
# qqnorm(resid(rategerm.b1)) #check normal distribution of residuals
# qqline(resid(rategerm.b1))
# ranef<- ranef(rategerm.b1, condVar = T) #extract conditional means
# ranef
# dotplot(ranef)
# which(residuals(rategerm.a1) > .006)


rategerm.b2 <- lmer(rategerm ~ trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run) + (1|region), 
                    data = sum, REML = F, na.action = na.omit)
summary(rategerm.b2)
rand(rategerm.b2)

plot(residuals(rategerm.b2))
hist(residuals(rategerm.b2))
plot(fitted(rategerm.b2), residuals(rategerm.b2))
qqnorm(resid(rategerm.b2)) #check normal distribution of residuals
qqline(resid(rategerm.b2))
ranef<- ranef(rategerm.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)

rategerm.c <- lmer(rategerm ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + 
                  (1+trt|cultivation) +(1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
                na.action = na.omit, data = sum, REML = F)

summary(rategerm.c)
rand(rategerm.c)

plot(residuals(rategerm.c))
hist(residuals(rategerm.c))
plot(fitted(rategerm.c), residuals(rategerm.c))
qqnorm(resid(rategerm.c)) #check normal distribution of residuals
qqline(resid(rategerm.c))

#run again without outliers <-- no change
x2 <- as.vector(which(residuals(rategerm.a2) > .005))
sum_rategerm <- sum[-x2,]

rategerm.a2.new <- lmer(rategerm ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run) + (1|cultivation) 
                        ,na.action = na.omit,data = sum_rategerm, REML = F)
summary(rategerm.a2.new)
rand(rategerm.a2.new) 


plot(residuals(rategerm.a2.new))
hist(residuals(rategerm.a2.new))
plot(fitted(rategerm.a2.new), residuals(rategerm.a2.new))
qqnorm(resid(rategerm.a2.new)) 
qqline(resid(rategerm.a2.new)) 
ranef<- ranef(rategerm.a2.new, condVar = T) 
ranef
dotplot(ranef)


rategerm.b2.new <- lmer(rategerm ~  trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run) + (1|region) 
                     ,data = sum_rategerm, REML = F, na.action = na.omit)
summary(rategerm.b2.new)
rand(rategerm.b2.new) 

plot(residuals(rategerm.b2.new))
hist(residuals(rategerm.b2.new))
plot(fitted(rategerm.b2.new), residuals(rategerm.b2.new))
qqnorm(resid(rategerm.b2.new)) #check normal distribution of residuals
qqline(resid(rategerm.b2.new)) #heavy tails
ranef<- ranef(rategerm.b2.new, condVar = T) #extract conditional means
ranef
dotplot(ranef)


rategerm.c.new <- lmer(rategerm ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + 
                      (1+trt|cultivation) +(1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
                    na.action = na.omit, data = sum_rategerm, REML = F)

summary(rategerm.c.new)
rand(rategerm.c.new)

plot(residuals(rategerm.c.new))
hist(residuals(rategerm.c.new))
plot(fitted(rategerm.c.new), residuals(rategerm.c.new))
qqnorm(resid(rategerm.c.new)) #check normal distribution of residuals
qqline(resid(rategerm.c.new)) #heavy tails
ranef<- ranef(rategerm.c.new, condVar = T) #extract conditional means
ranef
dotplot(ranef)


#####
#UNIFORMITY
#####
summary(sum$uniform)

# uniform.a1 <- lmer(uniform ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run), 
#               na.action = na.omit, data = sum, REML = F)
# summary(uniform.a1)
# rand(uniform.a1) 
# 
# par(mfrow=c(2,2))
# plot(residuals(uniform.a1))
# hist(residuals(uniform.a1))
# plot(fitted(uniform.a1), residuals(uniform.a1))
# qqnorm(resid(uniform.a1)) #check normal distribution of residuals
# qqline(resid(uniform.a1)) #heavy tails present
# ranef<- ranef(uniform.a1, condVar = T); ranef #extract conditional means
# dotplot(ranef)


uniform.a2 <- lmer(uniform ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run) + (1|cultivation), 
                  na.action = na.omit, data = sum, REML = F)
summary(uniform.a2)
rand(uniform.a2)

par(mfrow=c(2,2))
plot(residuals(uniform.a2))
hist(residuals(uniform.a2))
plot(fitted(uniform.a2), residuals(uniform.a2))
qqnorm(resid(uniform.a2)) #check normal distribution of residuals
qqline(resid(uniform.a2)) #heavy tails present
ranef<- ranef(uniform.a2, condVar = T); ranef #extract conditional means
dotplot(ranef)


# uniform.b1 <- lmer(uniform ~ trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
#                     data = sum, REML = F, na.action = na.omit)
# summary(uniform.b1)
# rand(uniform.b1)
# 
# plot(residuals(uniform.b1))
# hist(residuals(uniform.b1))
# plot(fitted(uniform.b1), residuals(uniform.b1))
# qqnorm(resid(uniform.b1)) #check normal distribution of residuals
# qqline(resid(uniform.b1))
# ranef<- ranef(uniform.b1, condVar = T) #extract conditional means
# ranef
# dotplot(ranef)


uniform.b2 <- lmer(uniform ~ trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run) #+ (1|region)
                    ,data = sum, REML = F, na.action = na.omit)
summary(uniform.b2)
rand(uniform.b2)

plot(residuals(uniform.b2))
hist(residuals(uniform.b2))
plot(fitted(uniform.b2), residuals(uniform.b2))
qqnorm(resid(uniform.b2)) #check normal distribution of residuals
qqline(resid(uniform.b2))
ranef<- ranef(uniform.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)


uniform.c <- lmer(uniform ~ trt + (1+trt|region) + (1+trt|region:landrace.name) +
                     (1+trt|cultivation) +(1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run),
                   na.action = na.omit, data = sum, REML = F)
summary(uniform.c)
rand(uniform.c)

plot(residuals(uniform.c))
hist(residuals(uniform.c))
plot(fitted(uniform.c), residuals(uniform.c))
qqnorm(resid(uniform.c)) #check normal distribution of residuals
qqline(resid(uniform.c))

#run again without outliers <-- no change
which(residuals(uniform.a2) > 300)
x3 <- as.vector(580)
sum_uniform <- sum[-x3,]


uniform.a2.new <- lmer(uniform ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run) + (1|cultivation) 
                       ,na.action = na.omit,data = sum_uniform, REML = F)
summary(uniform.a2.new)
rand(uniform.a2.new) 

plot(residuals(uniform.a2.new))
hist(residuals(uniform.a2.new))
plot(fitted(uniform.a2.new), residuals(uniform.a2.new))
qqnorm(resid(uniform.a2.new)) 
qqline(resid(uniform.a2.new)) 
ranef<- ranef(uniform.a2.new, condVar = T) 
ranef
dotplot(ranef)


uniform.b2.new <- lmer(uniform ~  trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run) #+ (1|region) 
                        ,data = sum_uniform, REML = F, na.action = na.omit)
summary(uniform.b2.new)
rand(uniform.b2.new) 

plot(residuals(uniform.b2.new))
hist(residuals(uniform.b2.new))
plot(fitted(uniform.b2.new), residuals(uniform.b2.new))
qqnorm(resid(uniform.b2.new)) #check normal distribution of residuals
qqline(resid(uniform.b2.new)) #heavy tails
ranef<- ranef(uniform.b2.new, condVar = T) #extract conditional means
ranef
dotplot(ranef)


uniform.c.new <- lmer(uniform ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + 
                         (1+trt|cultivation) +(1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
                       na.action = na.omit, data = sum_uniform, REML = F)

summary(uniform.c.new)
rand(uniform.c.new)

plot(residuals(uniform.c.new))
hist(residuals(uniform.c.new))
plot(fitted(uniform.c.new), residuals(uniform.c.new))
qqnorm(resid(uniform.c.new)) #check normal distribution of residuals
qqline(resid(uniform.c.new)) #heavy tails
ranef<- ranef(uniform.c.new, condVar = T) #extract conditional means
ranef
dotplot(ranef)

#####
#TOTAL PRECENT GERM
#####

summary(sum$perc_notgerm)

par(mfrow=c(1,1))
hist(sum$perc_germ)
sum$not10germ <- sum$perc_germ
  sum$not10germ <- as.numeric(sub(1, NA, sum$not10germ))
sum$did10germ <- sum$perc_germ
  sum$did10germ <- sum$did10germ==1
  sum$did10germ <- sum$did10germ*1
  
hist(sum$not10germ)

did10germ.a2 <- lmer(did10germ ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run) + (1|cultivation), 
                   na.action = na.omit, data = sum, REML = F)
summary(did10germ.a2)
rand(did10germ.a2)

par(mfrow=c(2,2))
plot(residuals(did10germ.a2))
hist(residuals(did10germ.a2))
plot(fitted(did10germ.a2), residuals(did10germ.a2))
qqnorm(resid(did10germ.a2)) #check normal distribution of residuals
qqline(resid(did10germ.a2)) #heavy tails present
ranef<- ranef(did10germ.a2, condVar = T)
ranef #extract conditional means
dotplot(ranef)


did10germ.b2 <- lmer(did10germ ~ trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run) + (1|region), 
                   data = sum, REML = F, na.action = na.omit)
summary(did10germ.b2)
rand(did10germ.b2)

plot(residuals(did10germ.b2))
hist(residuals(did10germ.b2))
plot(fitted(did10germ.b2), residuals(did10germ.b2))
qqnorm(resid(did10germ.b2)) #check normal distribution of residuals
qqline(resid(did10germ.b2))
ranef<- ranef(did10germ.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)


did10germ.c <- lmer(did10germ ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + 
                    (1+trt|cultivation) +(1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
                  na.action = na.omit, data = sum, REML = F)

summary(did10germ.c)
rand(did10germ.c)

not10germ.a2 <- lmer(not10germ ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + (1|rep) + (1|rep:run) + (1|cultivation), 
                     na.action = na.omit, data = sum, REML = F)
summary(not10germ.a2)
rand(not10germ.a2)

par(mfrow=c(2,2))
plot(residuals(not10germ.a2))
hist(residuals(not10germ.a2))
plot(fitted(not10germ.a2), residuals(not10germ.a2))
qqnorm(resid(not10germ.a2)) #check normal distribution of residuals
qqline(resid(not10germ.a2)) #heavy tails present
ranef<- ranef(not10germ.a2, condVar = T); ranef #extract conditional means
dotplot(ranef)


not10germ.b2 <- lmer(not10germ ~ trt + (1+trt|cultivation) + (1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run) + (1|region), 
                     data = sum, REML = F, na.action = na.omit)
summary(not10germ.b2)
rand(not10germ.b2)

plot(residuals(not10germ.b2))
hist(residuals(not10germ.b2))
plot(fitted(not10germ.b2), residuals(not10germ.b2))
qqnorm(resid(not10germ.b2)) #check normal distribution of residuals
qqline(resid(not10germ.b2))
ranef<- ranef(not10germ.b2, condVar = T) #extract conditional means
ranef
dotplot(ranef)


not10germ.c <- lmer(not10germ ~ trt + (1+trt|region) + (1+trt|region:landrace.name) + 
                      (1+trt|cultivation) +(1+trt|cultivation:landrace.name) + (1|rep) + (1|rep:run), 
                    na.action = na.omit, data = sum, REML = F)

summary(not10germ.c)
rand(not10germ.c)

plot(residuals(not10germ.c))
hist(residuals(not10germ.c))
plot(fitted(not10germ.c), residuals(not10germ.c))
qqnorm(resid(not10germ.c)) #check normal distribution of residuals
qqline(resid(not10germ.c))
ranef<- ranef(not10germ.c, condVar = T) #extract conditional means
ranef
dotplot(ranef)

#http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf