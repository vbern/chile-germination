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
df <- read.csv(paste(out.dir, "/cleaned_2017-09-07.csv", sep = ""), header = T)
str(df)
df$end <- round(df$end, digits = 1)
df$run <- as.factor(df$run)

df <- subset(df, region!="control")

#re-order data by region
summary(df$region)
target <- c("central valleys", "ecoast", "sierra madre", "wcoast", "yucatan")
library(gdata)
df$region <- reorder.factor(df$region, new.order=target)
summary(df$region)
df$trt <- as.factor(df$trt)
head(df)

library(survival)
library(survminer)
#fit and plot Kaplan-Meier survivor function, includes 95% confidence intervals
km.fit <-survfit(Surv(df$end, df$status) ~ 1, data = df, type = "kaplan-meier")
ggsurvplot(km.fit, conf.int = T)

test.region <- survfit(Surv(end, status) ~ region, data = df, subset = {trt == 0}, type = "kaplan-meier", conf.type = "log-log")
ggsurvplot(test.region, data =df, conf.int = T, pval = T, palette = "Dark2", risk.table = T, risk.table.col = "strata", xlab = "Time (h)", 
           surv.median.line = "hv", title = "Germination under Water Control", legend.labs = target)

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
library(devtools)
devtools::install_github("kassambara/survminer", build_vignettes = FALSE)
library(RColorBrewer)
ggsurvplot(test.cv, data = df, pval = T, title = "Central Valleys", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
df_cv <- subset(df, region=="central valleys")
survdiff_cv <- pairwise_survdiff(Surv(end, status)~trt, data = df_cv, rho = 0)
print(survdiff_cv)

#no significant difference between treatments
#survdiff_ecoast <- pairwise_survdiff(Surv(end, status)~trt, data = df, subset = {region == "ecoast"}, rho = 0)
#print(survdiff_ecoast)

ggsurvplot(test.wcoast, data = df, pval = T, title = "W Coast", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
df_wcoast <- subset(df, region == "wcoast")
survdiff_wcoast <-pairwise_survdiff(Surv(end, status)~trt, data = df_wcoast, rho = 0)
print(survdiff_wcoast)

ggsurvplot(test.sm, data = df, pval = T, title = "Sierra Madre de la Sur", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = rev(brewer.pal(4,"Spectral")), ggtheme = theme_grey())
df_sm <- subset(df, region == "sierra madre")
survdiff_sm <- pairwise_survdiff(Surv(end, status)~trt, data = df_sm, rho = 0)
print(survdiff_sm)

#no significant differences between treatments
#survdiff_yucatan <- survdiff(Surv(end, status)~trt, data = df10, subset = {region == "yucatan"}, rho = 0)
#print(survdiff_yucatan)

ggsurvplot(test.0, data = df, pval = T, title = "0% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
df_0 <- subset(df, trt==0)
survdiff_0 <- pairwise_survdiff(Surv(end, status)~region, data = df_0, rho = 0)
print(survdiff_0)
surv_median(test.0)

ggsurvplot(test.10, data = df, pval = T, title = "10% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
df_10 <- subset(df, trt==10)
survdiff_10 <- pairwise_survdiff(Surv(end, status)~region, data = df_10, rho = 0)
print(survdiff_10)
surv_median(test.10) 


ggsurvplot(test.15, data = df, pval = T, title = "15% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
df_15 <- subset(df, trt==15)
survdiff_15 <- pairwise_survdiff(Surv(end, status)~region, data = df_15, rho = 0)
print(survdiff_15)
surv_median(test.15) 


ggsurvplot(test.20, data = df, pval = T, title = "20% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
df_20 <- subset(df, trt==20)
survdiff_20 <- pairwise_survdiff(Surv(end, status)~region, data = df_20, rho = 0)
print(survdiff_20)
surv_median(test.20)

#TESTING DIFFERENCES ACROSS CULTIVATION SYSTEMS
test.cult.0 <- survfit(Surv(end, status) ~ cultivation, data = df, type = "kaplan-meier", subset = (trt ==0))
ggsurvplot(test.cult.0, data = df, pval = T, risk.table = T, title = "0% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
survdiff_cult_0 <- pairwise_survdiff(Surv(end, status)~cultivation, data = df_0, rho = 0)
print(survdiff_cult_0)

test.cult.20 <- survfit(Surv(end, status) ~ cultivation, data = df, type = "kaplan-meier", subset = (trt ==20))
ggsurvplot(test.cult.20, data = df, pval = T, risk.table = T, title = "20% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
survdiff_cult_20 <- pairwise_survdiff(Surv(end, status)~cultivation, data = df_20, rho = 0)
print(survdiff_cult_20)

#TESTING DIFFERENCES ACROSS DOMESTICATION LEVEL
test.domestication.0 <- survfit(Surv(end, status) ~ population.type, data = df, type = "kaplan-meier", subset = (trt ==0))
ggsurvplot(test.domestication.0, data = df, pval = T, risk.table = T, title = "0% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
survdiff_domestication_0 <- pairwise_survdiff(Surv(end, status)~population.type, data = df_0, rho = 0)
print(survdiff_domestication_0)

test.domestication.20 <- survfit(Surv(end, status) ~ population.type, data = df, type = "kaplan-meier", subset = (trt ==20))
ggsurvplot(test.domestication.20, data = df, pval = T, risk.table = T, title = "20% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", palette = "Dark2", ggtheme = theme_grey())
survdiff_domestication_20 <- pairwise_survdiff(Surv(end, status)~population.type, data = df_20, rho = 0)
print(survdiff_domestication_20)


#TESTING DIFFERENT LANDRACES
df_0 <- subset(df, trt==0)
summary(df$landrace.name)
target <- c()
df$landrace.name <- reorder.factor(df$landrace.name)
summary(df$landrace.name)

summary(df$landrace.name, by = region)

 test.lr.0 <- survfit(Surv(end, status) ~ landrace.name, data = df, type = "kaplan-meier", subset = (trt ==0))
ggsurvplot(test.lr.0, data = df, pval = T, risk.table = T, title = "0% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", ggtheme = theme_grey())
survdiff_lr_0 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_0, rho = 0)
print(survdiff_lr_0)

df_20 <- subset(df, trt==20)
test.lr.20 <- survfit(Surv(end, status) ~ landrace.name, data = df, type = "kaplan-meier", subset = (trt ==20))
ggsurvplot(test.lr.20, data = df, pval = T, risk.table = T, title = "20% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", ggtheme = theme_grey())
survdiff_lr_20 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_20, rho = 0)
print(survdiff_lr_20)

df_15 <- subset(df, trt==15)
test.lr.15 <- survfit(Surv(end, status) ~ landrace.name, data = df, type = "kaplan-meier", subset = (trt ==15))
ggsurvplot(test.lr.15, data = df, pval = T, risk.table = T, title = "15% PEG", conf.int = T, legend = "bottom", xlab = "Time (h)", ggtheme = theme_grey())
survdiff_lr_15 <- pairwise_survdiff(Surv(end, status)~landrace.name, data = df_15, rho = 0)
print(survdiff_lr_15)

df_cv <- subset(df,region == "central valleys") #landraces == Chile de Agua, Taviche, Tusta
df_sm <- subset(df, region == "sierra madre") #landraces ==Tusta
df_wcoast <- subset(df, region == "wcoast") #landraces == Costeño Amarillo, Costeño Rojo, Piquin
df_ecoast <- subset(df, region == "ecoast") #landraces == Chigole, Chile Bolita, Chile de Monte, Costeno Rojo, Guajillo, Guina Dahni, Mareno, Mirasol, Payaso, Tusta
df_yucatan <- subset(df, region == "yucatan") #landraces == Dulce, Pardito
summary(df_yucatan$landrace.name)
