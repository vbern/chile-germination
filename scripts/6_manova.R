#' ---
#' title: "Germination Analysis: multivariate"
#' author: "Vivian bernau"
#' date: "23 January 2018"
#' output: word_document
#' ---

#MANOVA
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

setwd(data.dir)

germ <- read.csv("summarydata_2018-01-23.csv", header = T)

library(dplyr)

germ_melt <- melt(germ, idv.vars = c(1:2, 16:32), measure.vars = c("delay", "rategerm", "uniform", "perc_germ"))
summary(germ_melt$variable)
summary(germ_melt$perc_notgerm)


# Omit rows with missings in response variable:
germ_melt <- germ_melt[!is.na(germ_melt$value),]

library(nlme)

manova1 <- lme(value ~ 0 + variable,
              random = 1 + variable|trt + trt|region + 1|region:landrace.name + 0+trt|region:landrace.name + 1|rep + 1|rep:run
              ,data=germ_melt)

manova1 <- lme(value ~ 0 + variable,
               random = list(region = pdDiag(~trt), )
               ,data = germ_melt)
summary(manova1) 
anova(manova1)
random.effects(manova1)


test1 <- manova(cbind(delay, rategerm, uniform, perc_germ) ~ trt, data = germ)
summary(test1)

test2 <- lme(value ~ variable*trt, data = germ_melt)
