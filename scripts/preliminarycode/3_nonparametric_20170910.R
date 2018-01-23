#GERMINATION DATA
#VIVIAN BERNAU
#6 Sept 2017
#Based on McNair et al 2012; Seed Science Research

#Using SURVMINER

#NONPARAMETRIC ANALYSIS
#1) Characterizing pattern of germination within groups --> Kaplan-Meier test: survfit()
#2) Comparing patterns of germination in groups --> Fleming-Harrington test(?): survdiff()

#My data is right-sensored.
#My data is interval, but can probably be analyzed nonparametrically as exact for as long as plates with large losses are removed.

#Survivor function: probability that the germination time is greater than t

#Set working directory and repositories
wd <- ("~/Google Drive/RFiles/Chiles_OSU/3_Germination/")

src.dir <- paste(wd,"scripts", sep = "")
data.dir <- paste (wd,"data", sep = "")
out.dir <- paste(wd, "output", sep ="")

setwd(out.dir)

#read in germination data in pre-lifetab format
#df10 <- read.csv(paste(out.dir, "/cleaned10_2017-09-07.csv", sep = ""), header = T)
df <- read.csv(paste(out.dir, "/cleaned_2017-09-10.csv", sep = ""), header = T)
str(df)
df$end <- round(df$end, digits = 1)
df$run <- as.factor(df$run)

df <- subset(df, region!="control")
df$cv <- as.numeric(df$region == "central valleys")
df$ecoast <- as.numeric(df$region == "ecoast")
df$wcoast <- as.numeric(df$region == "wcoast")
df$yucatan <- as.numeric(df$region == "yucatan")
df$sm <- as.numeric(df$region == "sierra madre")

library(survival)
library(survminer)
#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df, type = "kaplan-meier")
ggsurvplot(km.fit, conf.int = T)

test.region <- survfit(Surv(end, status) ~ region, data = df, subset = {trt == 0}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.region, data =df, conf.int = T, pval = T, palette = "Dark2", risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under Water Control")

#Creating composite figures of with one frame for each [PEG]
# List of ggsurvplots
test.20<- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier",  subset = (trt == 20))
test.15 <- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier", subset = (trt ==15))
test.10 <- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier", subset = (trt ==10))
test.0 <- survfit(Surv(end, status) ~ region, data = df, type = "kaplan-meier", subset = (trt ==0))

pegplots <- list()
pegplots[[1]] <- ggsurvplot(test.20, data = df, pval = T, conf.int = T, palette = "Dark2", legend = "none", xlab = "Time (h)", title = "20% [PEG]", ggtheme = theme_grey())
pegplots[[2]] <- ggsurvplot(test.15, data = df, pval = T, conf.int = T, palette = "Dark2", legend = "none", xlab = "Time (h)", title = "15% [PEG]", ggtheme = theme_grey())
pegplots[[3]] <- ggsurvplot(test.10, data = df, pval = T, conf.int = T, palette = "Dark2", legend = "none", xlab = "Time (h)", title = "10% [PEG]", ggtheme = theme_grey())
pegplots[[4]] <- ggsurvplot(test.0, data = df, pval = T, conf.int = T, palette = "Dark2", legend = "none", xlab = "Time (h)", title = "0% [PEG]" , ggtheme = theme_grey())

# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(pegplots, title = "Germination under different [PEG]",
                    ncol = 2, nrow = 2)
res <- arrange_ggsurvplots(pegplots, print = FALSE)
ggsave("allconcentrations_20170908.png", res)

#Creating composite figure of Kaplan-Meier survivor estimates with 1 frame for each region
#png(filename="allregions20170906.png", width = 700, height = 1400, res = 120)
test.cv <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="central valleys"))
test.ecoast <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="ecoast"))
test.sm <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="sierra madre"))
test.wcoast <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="wcoast"))
test.yucatan <- survfit(Surv(end, status) ~ trt, data = df, type = "kaplan-meier", subset = (region =="yucatan"))

library(RColorBrewer)
regplots <- list()
regplots[[1]] <- ggsurvplot(test.cv, data = df, pval = T, title = "Central Valleys", conf.int = T, legend = "none", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
regplots[[2]] <- ggsurvplot(test.ecoast, data = df, pval = T, title = "E Coast", conf.int = T, legend = "none", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
regplots[[3]] <- ggsurvplot(test.sm, data = df, pval = T, title = "Sierra Madre de la Sur", conf.int = T, legend = "none", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
regplots[[4]] <- ggsurvplot(test.wcoast, data = df, pval = T, title = "W Coast", conf.int = T, legend = "none", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
regplots[[5]] <- ggsurvplot(test.yucatan, data = df, pval = T, title = "Yucatan", conf.int = T, legend = "right", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())

# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(regplots, print = TRUE, title = "germination from different regions",
                    ncol = 2, nrow = 3)
res <- arrange_ggsurvplots(regplots, print = FALSE)
ggsave("allconcentrations_20170907.png", res)

#COMPARISON OF SURVIVOR FUNCTIONS
library(RColorBrewer)
ggsurvplot(test.cv, data = df, pval = T, title = "Central Valleys", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
df_cv <- subset(df, region=="central valleys")
survdiff_cv <- pairwise_survdiff(Surv(end, status)~trt, data = df_cv, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_cv)

#no significant difference between treatments
#survdiff_ecoast <- pairwise_survdiff(Surv(end, status)~trt, data = df, subset = {region == "ecoast"}, rho = 0)
#print(survdiff_ecoast)

ggsurvplot(test.wcoast, data = df, pval = T, title = "W Coast", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
df_wcoast <- subset(df, region == "wcoast")
survdiff_wcoast <-pairwise_survdiff(Surv(end, status)~trt, data = df_wcoast, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_wcoast)

ggsurvplot(test.sm, data = df, pval = T, title = "Sierra Madre de la Sur", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
df_sm <- subset(df, region == "sierra madre")
survdiff_sm <- pairwise_survdiff(Surv(end, status)~trt, data = df_sm, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_sm)

#no significant differences between treatments
#survdiff_yucatan <- survdiff(Surv(end, status)~trt, data = df10, subset = {region == "yucatan"}, rho = 0)
#print(survdiff_yucatan)

ggsurvplot(test.0, data = df, pval = T, title = "0% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
df_0 <- subset(df, trt==0)
survdiff_0 <- pairwise_survdiff(Surv(end, status)~region, data = df_0, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_0)
surv_median(test.0)

ggsurvplot(test.10, data = df, pval = T, title = "10% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
df_10 <- subset(df, trt==10)
survdiff_10 <- pairwise_survdiff(Surv(end, status)~region, data = df_10, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_10)
surv_median(test.10) 

ggsurvplot(test.15, data = df, pval = T, title = "15% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
df_15 <- subset(df, trt==15)
survdiff_15 <- pairwise_survdiff(Surv(end, status)~region, data = df_15, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_15)
surv_median(test.15) 

ggsurvplot(test.20, data = df, pval = T, title = "20% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
df_20 <- subset(df, trt==20)
survdiff_20 <- pairwise_survdiff(Surv(end, status)~region, data = df_20, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_20)
surv_median(test.20)

#EXTRACT .25 AND .75 QUANTILES
data.frame(q = c(.25, .5, .75),
                  km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_0, cv == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_0, wcoast == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_0, sm == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_0, ecoast == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_0, yucatan == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_10, cv == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_10, wcoast == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_10, sm == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_10, ecoast == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_10, yucatan == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_15, cv == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_15, wcoast == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_15, sm == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_15, ecoast == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_15, yucatan == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_20, cv == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_20, wcoast == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_20, sm == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_20, ecoast == 1 ), type = "kaplan-meier")))
data.frame(q = c(.25, .5, .75),
           km = quantile(survfit(Surv(end,status) ~ 1, data = subset(df_20, yucatan == 1 ), type = "kaplan-meier")))

#TESTING DIFFERENCES ACROSS CULTIVATION SYSTEMS
test.cult.0 <- survfit(Surv(end, status) ~ cultivation, data = df, type = "kaplan-meier", subset = (trt ==0))
ggsurvplot(test.cult.0, data = df, pval = T, title = "0% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
survdiff_cult_0 <- pairwise_survdiff(Surv(end, status)~cultivation, data = df_0, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_cult_0)

test.cult.20 <- survfit(Surv(end, status) ~ cultivation, data = df, type = "kaplan-meier", subset = (trt ==20))
ggsurvplot(test.cult.20, data = df, pval = T, title = "20% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
survdiff_cult_20 <- pairwise_survdiff(Surv(end, status)~cultivation, data = df_20, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_cult_20)

#TESTING DIFFERENCES ACROSS DOMESTICATION LEVEL
test.domestication.0 <- survfit(Surv(end, status) ~ population.type, data = df, type = "kaplan-meier", subset = (trt ==0))
ggsurvplot(test.domestication.0, data = df, pval = T, title = "0% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
survdiff_domestication_0 <- pairwise_survdiff(Surv(end, status)~population.type, data = df_0, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_domestication_0)

test.domestication.20 <- survfit(Surv(end, status) ~ population.type, data = df, type = "kaplan-meier", subset = (trt ==20))
ggsurvplot(test.domestication.20, data = df, pval = T, title = "20% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
survdiff_domestication_20 <- pairwise_survdiff(Surv(end, status)~population.type, data = df_20, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_domestication_20)


#TESTING DIFFERENT LANDRACES
df_0 <- subset(df, trt==0)
summary(df$landrace.name)

test.lr.0 <- survfit(Surv(end, status) ~ landrace.name, data = df, type = "kaplan-meier", subset = (trt ==0))
ggsurvplot(test.lr.0, data = df, pval = T, title = "0% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", ggtheme = theme_grey())
survdiff_lr_0 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_0, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_lr_0)

df_ecoast <- subset(df, region =="ecoast")
df_ecoast_0 <- subset(df_ecoast, trt == 0)
test.landrace.ecoast <- survfit(Surv(end, status) ~ landrace.name, data = df_ecoast, subset = {trt == 0}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.landrace.ecoast, data =df, conf.int = T, pval = T,  xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under Water Control")
survdiff_lr_ecoast_0 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_ecoast_0,p.adjust.method = "bonferroni", rho = 0)
print(survdiff_lr_ecoast_0)

df_cv <- subset(df, region =="central valley")
#test.landrace.cv <- survfit(Surv(end, status) ~ landrace.name, data = df_cv, subset = {trt == 0}, type = "kaplan-meier", conf.type = "log-log")
#ggsurvplot(test.landrace.cv, data =df, conf.int = T, pval = T, risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
#           surv.median.line = "hv", title = "Germination under Water Control")
df_wcoast <- subset(df, region =="wcoast")
df_wcoast_0 <- subset(df_wcoast, trt ==0)
test.landrace.wcoast <- survfit(Surv(end, status) ~ landrace.name, data = df_wcoast, subset = {trt == 0}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.landrace.wcoast, data =df, conf.int = T, pval = T,  xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under Water Control")
survdiff_lr_wcoast_0 <- pairwise_survdiff(Surv(end, status) ~ landrace.name, data = df_wcoast, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_lr_wcoast_0)
df_yucatan <- subset(df, region =="yucatan")
df_yucatan_0 <- subset(df_yucatan, trt ==0)
test.landrace.yucatan <- survfit(Surv(end, status) ~ landrace.name, data = df_yucatan, subset = {trt == 0}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.landrace.yucatan, data =df, conf.int = T, pval = T, xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under Water Control")
survdiff_lr_yucatan_0 <- pairwise_survdiff(Surv(end, status) ~ landrace.name, data = df_yucatan, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_lr_yucatan_0)
df_sm <- subset(df, region =="sierra madre")
#df_sm_0 <- subset(df_sm, trt == "trt")
#test.landrace.sm <- survfit(Surv(end, status) ~ landrace.name, data = df_sm, subset = {trt == 0}, type = "kaplan-meier", conf.type = "log-log")
#ggsurvplot(test.landrace.sm, data =df, conf.int = T, pval = T, risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
#           surv.median.line = "hv", title = "Germination under Water Control")

df_20 <- subset(df, trt==20)
test.lr.20 <- survfit(Surv(end, status) ~ landrace.name, data = df, type = "kaplan-meier", subset = (trt ==20))
ggsurvplot(test.lr.20, data = df, pval = T, title = "20% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", ggtheme = theme_grey())
survdiff_lr_20 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_20, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_lr_20)

df_ecoast_20 <- subset(df_ecoast, trt == 20)
test.landrace.ecoast.20 <- survfit(Surv(end, status) ~ landrace.name, data = df_ecoast, subset = {trt == 20}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.landrace.ecoast.20, data =df, conf.int = T, pval = T, xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under 20% PEG")
survdiff_lr_ecoast_20 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_ecoast_20, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_lr_ecoast_20)
#test.landrace.cv.20 <- survfit(Surv(end, status) ~ landrace.name, data = df_cv, subset = {trt == 20}, type = "kaplan-meier", conf.type = "log-log")
#ggsurvplot(test.landrace.cv.20, data =df, conf.int = T, pval = T, risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
#           surv.median.line = "hv", title = "Germination under 20% PEG")
df_wcoast_20 <- subset(df_wcoast, trt == 20)
test.landrace.wcoast.20 <- survfit(Surv(end, status) ~ landrace.name, data = df_wcoast, subset = {trt == 20}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.landrace.wcoast.20, data =df, conf.int = T, pval = T, xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under 20% PEG")
survdiff_lr_wcoast_20 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_wcoast_20, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_lr_wcoast_20)
#test.landrace.yucatan.20 <- survfit(Surv(end, status) ~ landrace.name, data = df_yucatan, subset = {trt == 20}, type = "kaplan-meier", conf.type = "log-log")
#ggsurvplot(test.landrace.yucatan.20, data =df, conf.int = T, pval = T, risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
#           surv.median.line = "hv", title = "Germination under 20% PEG"")
#test.landrace.sm.20 <- survfit(Surv(end, status) ~ landrace.name, data = df_sm, subset = {trt == 20}, type = "kaplan-meier", conf.type = "log-log")
#ggsurvplot(test.landrace.sm.20, data =df, conf.int = T, pval = T, risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
#           surv.median.line = "hv", title = "Germination under 20% PEG")

df_15 <- subset(df, trt==15)
test.lr.15 <- survfit(Surv(end, status) ~ landrace.name, data = df, type = "kaplan-meier", subset = (trt ==15))
ggsurvplot(test.lr.15, data = df, pval = T, risk.table = T, title = "15% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", 
           ggtheme = theme_grey())
survdiff_lr_15 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_15, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_lr_15)

df_ecoast_15 <- subset(df_ecoast, trt == 15)
test.landrace.ecoast.15 <- survfit(Surv(end, status) ~ landrace.name, data = df_ecoast, subset = {trt == 15}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.landrace.ecoast.15, data =df, conf.int = T, pval = T,  xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under 15% PEG")
survdiff_lr_ecoast_15 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_ecoast_15, p.adjust.method = "bonferroni",rho = 0)
print(survdiff_lr_ecoast_15)
#test.landrace.cv.15 <- survfit(Surv(end, status) ~ landrace.name, data = df_cv, subset = {trt == 15}, type = "kaplan-meier", conf.type = "log-log")
#ggsurvplot(test.landrace.cv.15, data =df, conf.int = T, pval = T, risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
#           surv.median.line = "hv", title = "Germination under 15% PEG")
df_wcoast_15 <- subset(df_wcoast, trt == 15)
test.landrace.wcoast.15 <- survfit(Surv(end, status) ~ landrace.name, data = df_wcoast, subset = {trt == 15}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.landrace.wcoast.15, data =df, conf.int = T, pval = T,  xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under 15% PEG")
survdiff_lr_wcoast_15 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_wcoast_15, p.adjust.method = "bonferroni", rho = 0)
print(survdiff_lr_wcoast_15)

df_yucatan_15 <- subset(df_yucatan, trt ==15)
test.landrace.yucatan.15 <- survfit(Surv(end, status) ~ landrace.name, data = df_yucatan, subset = {trt == 15}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.landrace.yucatan.15, data =df, conf.int = T, pval = T, xlab = "Time (h)", 
           surv.median.line = "hv", title ="Germination under 15% PEG")
survdiff_lr_yucatan_15 <- pairwise_survdiff(Surv(end, status) ~ landrace.name, data = df_yucatan_15, p.adjust.method = "bonferroni", rho = 0)
print(survdiff_lr_yucatan_15)
#test.landrace.sm.15 <- survfit(Surv(end, status) ~ landrace.name, data = df_sm, subset = {trt == 15}, type = "kaplan-meier", conf.type = "log-log")
#ggsurvplot(test.landrace.sm.15, data =df, conf.int = T, pval = T, risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
#           surv.median.line = "hv", title = "Germination under Water Control")