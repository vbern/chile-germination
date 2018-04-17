
#' ---
#' title: "Germination Analysis: heritability"
#' author: "Vivian Bernau"
#' date: "29 September 2017"
#' output: word_document
#' ---


wd <- ("~/Google Drive/RFiles/chile-germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

library(lme4)

par(mfrow=c(2,2))

sum <- read.csv(file = paste(data.dir, "/summarydata_2018-01-23.csv", sep = ""), na.strings=c(""," ","NA"), header =T)
a <- subset(sum, !is.na(landrace.name)); summary(a$region)
b <- subset(a, !is.na(region)); summary(b$region)
c <- subset(b, !is.na(cultivation)); summary(c$region)
sum <- c

sum$trt <- as.character(sum$trt)

sum$not10germ <- sum$perc_germ
sum$not10germ <- as.numeric(sub(1, NA, sum$not10germ))
sum$did10germ <- sum$perc_germ
sum$did10germ <- sum$did10germ==1
sum$did10germ <- sum$did10germ*1

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
Var <- as.data.frame(VarCorr(delay.h))
row.names(Var) = Var$grp
Var["landrace.name","vcov"]/(Var["landrace.name","vcov"] + Var["trt:landrace.name","vcov"]/4 +attr(VarCorr(delay.h), "sc")/4)
repeatability(sum$delay, sum$landrace.name, line.repeatability = T)
repeatability(sum$delay, sum$sampleid, line.repeatability = T)
repeatability(sum$delay, sum$line, line.repeatability = T)

uniformity.h <- lmer(uniform ~  (1|trt) + (1|landrace.name) + (1|trt:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(uniformity.h)
Var <- as.data.frame(VarCorr(uniformity.h))
row.names(Var) = Var$grp
Var["landrace.name","vcov"]/(Var["landrace.name","vcov"] + Var["trt:landrace.name","vcov"]/4 +attr(VarCorr(delay.h), "sc")/4)
repeatability(sum$uniform, sum$landrace.name, line.repeatability = T)
repeatability(sum$uniform, sum$sampleid, line.repeatability = T)
repeatability(sum$uniform, sum$line, line.repeatability = T)

rate.h <- lmer(rategerm ~  (1|trt) + (1|landrace.name) + (1|trt:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(rate.h)
Var <- as.data.frame(VarCorr(rate.h))
row.names(Var) = Var$grp
Var["landrace.name","vcov"]/(Var["landrace.name","vcov"] + Var["trt:landrace.name","vcov"]/4 +attr(VarCorr(rate.h), "sc")/4)
repeatability(sum$rategerm, sum$landrace.name, line.repeatability = T)
repeatability(sum$rategerm, sum$sampleid, line.repeatability = T)
repeatability(sum$rategerm, sum$line, line.repeatability = T)

did10germ.h <- lmer(did10germ ~ (1|trt) + (1|landrace.name) + (1|trt:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(did10germ.h)
Var <- as.data.frame(VarCorr(did10germ.h))
row.names(Var) = Var$grp
Var["landrace.name","vcov"]/(Var["landrace.name","vcov"] + Var["trt:landrace.name","vcov"]/4 +attr(VarCorr(did10germ.h), "sc")/4)
repeatability(sum$did10germ, sum$landrace.name, line.repeatability = T)
repeatability(sum$did10germ, sum$sampleid, line.repeatability = T)
repeatability(sum$did10germ, sum$line, line.repeatability = T)

not10germ.h <- lmer(not10germ ~ (1|trt) + (1|landrace.name) + (1|trt:landrace.name) + (1|rep) + (1|rep:run), data = sum, REML = T)
summary(not10germ.h)
Var <- as.data.frame(VarCorr(not10germ.h))
row.names(Var) = Var$grp
Var["landrace.name","vcov"]/(Var["landrace.name","vcov"] + Var["trt:landrace.name","vcov"]/4 +attr(VarCorr(not10germ.h), "sc")/4)
repeatability(sum$did10germ, sum$landrace.name, line.repeatability = T)
repeatability(sum$did10germ, sum$sampleid, line.repeatability = T)
repeatability(sum$did10germ, sum$line, line.repeatability = T)
