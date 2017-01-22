setwd("C:/Users/Andrea Mack/Desktop/mack_hub/writing_project/data")
x <- read.csv("SEDataForAnalysisSpSurveyFinal.csv")

head(x)
dim(x)
summary(x)

cbind(as.character(x$EvalWater), x$SeWater)
#nn -> no data
cbind(as.character(sx$PlantType),as.character(x$SeSed),as.character(x$SeRoots))

# what is the difference bw xcoord and coords.x1, likewise ycoord and coords.x2?
plot(x$xcoord~x$coords.x1)
ratio1 <- x$xcoord/x$coords.x1
summary(ratio1)

ratio2 <- x$ycoord/x$coords.x2
ratio2 <- x$ycoord/x$coords.x2
summary(ratio2)
# these are the same, so let's drop coords.x1 and coords.x2
drop <- names(x) %in% c("coords.x1", "coords.x2")
x <- x[,!drop]
dim(x)

summary(x)

##############
# codebook #
##############
#X same as siteID, will delete
# siteID = unique site identifier (2 samples per site? -- yes see protocol p.12)
#X.x?
# xcoord = x coordinate
# ycoord = y coordinate

plot(mdcaty.x~mdcaty.y, data = x)
# mdcaty.x = type of flooding
#wgt
# stratum.x = zone
#EvalStatus.x = all NotEval, what is this? I will delete it
#EvalReason = all NA's I deleted
drop <- names(x) %in% c("X", "EvalReason", "EvalStatus.x")
x <- x[,!drop]
dim(x)
summary(x)
#wgt
#init.wgt
#what is the difference?
#which was used to calculate
#adj weight for root, sediment, bugs, water 

plot(x$wgt~x$init.wgt)
summary(x$adj.wgt.sed/x$wgt)
summary(x$adj.wgt.sed/x$init.wgt)
summary(x)

pairs(x$OBJECTID_1~x$OBJECTID.x+x$OBJECTID.y)
apply(data.frame(cbind(x$OBJECTID_1, x$OBJECTID.x, x$OBJECTID.y)), 2, summary)
#OBJECTID_1
#OBJECTID.x
#OBJECTID.y
#don't know what these are

plot(x$FID_bnl_Se~x$FID_bnl_ba)

#FID_bnl_Se 
#FID_bnl_ba
#both of these appear to be factors, one with 5 levels (se) the other with 8 levels (ba)

plot(x$Shape_Leng~x$Shape_Le_1)
plot(x$Area~x$Shape_Area)
plot(x$Shape_Leng*x$Shape_Le_1~x$Shape_Area)

#Shape_Leng
#Shape_Le_1
#Area
#Shape_Area
#no idea relationship

plot(x$Acres~x$Hectares)
# Acres
# Hectares
# these are directly proportional

plot(x$ycoord~x$X.y)
#I think all of the .y and .x variables are related to the x and y coordinate fields (lat and lon)

table(x$OBJECTID.x,x$X.y)
table(x$OBJECTID.y,x$X.x)[1:5,1:5]
length(B <- which(x$OBJECTID.x %in% x$X.y==FALSE))
length(C <- which(x$OBJECTID.y %in% x$X.x==FALSE))
#there seems to be some relationship between objectid.x and x.y and vv

summary(x)
length(A <- which(x$X.x %in% x$F1== FALSE))
dim(x)
x$X.x[A]
x$F1[A]
#X.x and F1 seems similar also, except for six observations
C %in% A

# Point reference data
# Date_

#Unknown
#SampleType (S; S,I; W,S,I)
#Sedtype (roots; sed)

x$LabSampleNo
# I-> invertebre
#I believe LabSampleNo related to X.x, but with added what sample was taken on (S,I,W,E)

# Each of the .Water .Bugs and .Sed (and .Root) have the following columns:
# Eval (NN,NT,TS)
# Se.
#.TN
#.NT
#.TS

######## from Jay for Vanessa #########
# eval refers to whether the desired speciman was AVAILABLE to be sampled, not whether measurements were actually
# taken at that point

# TS,TN, NT  == evaluated
#TS - target and sampled
#TN - target and not sampled
#NT - not target and not sampled

# NN, NE, NA == not evaluated

## Appendix B -- R code preparation taken from Jay ##
SE.Data <- read.csv("bnl.se.data.all.clean.csv", header = T)
SE.Data <- subset(SE.Data, stratum != "Unit 4A-T")
SE.Data <- SE.Data[-((nrow(SE.Data) - 3):nrow(SE.Data)), ]
SE.Data$SeSed <- as.numeric(ifelse(SE.Data$SeSed == "noData", NA, as.character(SE.Data$SeSed)))
SE.Data$SeWater <- as.numeric(ifelse(SE.Data$SeWater == "noData", NA, as.character(SE.Data$SeWater)))
SE.Data$SeBugs <- as.numeric(ifelse(SE.Data$SeBugs == "noData", NA, as.character(SE.Data$SeBugs)))
SE.Data$SeRoots <- as.numeric(ifelse(SE.Data$SeRoots == "noData", NA, as.character(SE.Data$SeRoots)))
SE.Data <- subset(SE.Data, EvalSed != "NA")
####################################### SEDIMENT ESTIMATES ############
FixedSeDataAnalysis <- read.csv("SEDataForAnalysisSpSurveyFinal.csv")
FixedSeDataAnalysis <- subset(FixedSeDataAnalysis, EvalSed == "TS")
SedSites <- data.frame(siteID = FixedSeDataAnalysis$siteID, Active = rep(TRUE,
length(FixedSeDataAnalysis$SeSed)))
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = non.area.adj.wgt.sed, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
SedData <- with(FixedSeDataAnalysis, data.frame(siteID = siteID, SeSed))
subpop <- FixedSeDataAnalysis[, c(2, 9)]
colnames(subpop)[2] <- "ManUnit"
sizes <- with(FixedSeDataAnalysis, tapply(init.area, stratum.x, mean))
SedPopSizes.Adj <- list(ManUnit = as.list(sizes))
# using non area adjusted weights
Sed.Est.Local.init <- cont.analysis(sites = SedSites, subpop = subpop,
design = AreaAdj.Design, data.cont = SedData, total = T, popsize = SedPopSizes.Adj,
popcorrect = T, pcfsize = sum(sizes), vartype = "local", conf = 90,
pctval = c(0))
# adjusting using estimated area
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = adj.wgt.sed, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
Sed.Est.Local.adj <- cont.analysis(sites = SedSites, subpop = subpop, design = AreaAdj.Design,
data.cont = SedData, total = T, popsize = SedPopSizes.Adj, popcorrect = T,
pcfsize = sum(sizes), vartype = "local", conf = 90, pctval = c(0))
# adjusting using real area
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = real.wgt.sed, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
Sed.Est.Local.real <- cont.analysis(sites = SedSites, subpop = subpop,
design = AreaAdj.Design, data.cont = SedData, total = T, popsize = SedPopSizes.Adj,
popcorrect = T, pcfsize = sum(sizes), vartype = "local", conf = 90,
pctval = c(0))
# ignoring unequal wts and assuming JUST A SRS within each stratum
24
# [WRONG!!]
Simple.SRS.sed <- ddply(FixedSeDataAnalysis, ~stratum.x, summarize, SampleSize = length(SeSed),
Mean = mean(SeSed), SE = sd(SeSed)/sqrt(length(SeSed)))
maxes.sed <- as.matrix(with(FixedSeDataAnalysis, tapply(SeSed, stratum.x,
max)))
####################################### ROOT ESTIMATES ############
FixedSeDataAnalysis <- read.csv("SEDataForAnalysisSpSurveyFinal.csv")
FixedSeDataAnalysis <- subset(FixedSeDataAnalysis, EvalRoot == "TS")
RootSites <- data.frame(siteID = FixedSeDataAnalysis$siteID, Active = rep(TRUE,
length(FixedSeDataAnalysis$SeRoots)))
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = non.area.adj.wgt.root, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
RootData <- with(FixedSeDataAnalysis, data.frame(siteID = siteID, SeRoots))
subpop <- FixedSeDataAnalysis[, c(2, 9)]
colnames(subpop)[2] <- "ManUnit"
sizes <- with(FixedSeDataAnalysis, tapply(init.area, stratum.x, mean))
sizes <- sizes[c(2, 3, 4, 5, 7)]
RootPopSizes.Adj <- list(ManUnit = as.list(sizes))
# using non area adjusted weights
Root.Est.Local.init <- cont.analysis(sites = RootSites, subpop = subpop,
design = AreaAdj.Design, data.cont = RootData, total = T, popsize = RootPopSizes.Adj,
popcorrect = T, pcfsize = sum(sizes), vartype = "local", conf = 90,
pctval = c(0))
# adjusting using estimated area
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = adj.wgt.root, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
Root.Est.Local.adj <- cont.analysis(sites = RootSites, subpop = subpop,
design = AreaAdj.Design, data.cont = RootData, total = T, popsize = RootPopSizes.Adj,
popcorrect = T, pcfsize = sum(sizes), vartype = "local", conf = 90,
pctval = c(0))
# adjusting using real area
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = real.wgt.root, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
Root.Est.Local.real <- cont.analysis(sites = RootSites, subpop = subpop,
design = AreaAdj.Design, data.cont = RootData, total = T, popsize = RootPopSizes.Adj,
popcorrect = T, pcfsize = sum(sizes), vartype = "local", conf = 90,
pctval = c(0))
# ignoring unequal wts and assuming JUST A SRS within each stratum
# [WRONG!!]
Simple.SRS.root <- ddply(FixedSeDataAnalysis, ~stratum.x, summarize, SampleSize = length(SeRoots),
Mean = mean(SeRoots), SE = sd(SeRoots)/sqrt(length(SeRoots)))
maxes.root <- as.matrix(with(FixedSeDataAnalysis, tapply(SeRoots, stratum.x,
max)))
25



####################################### WATER ESTIMATES ############
FixedSeDataAnalysis <- read.csv("SEDataForAnalysisSpSurveyFinal.csv")
FixedSeDataAnalysis <- subset(FixedSeDataAnalysis, EvalWater == "TS")
WaterSites <- data.frame(siteID = FixedSeDataAnalysis$siteID, Active = rep(TRUE,
length(FixedSeDataAnalysis$SeWater)))
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = non.area.adj.wgt.water, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
WaterData <- with(FixedSeDataAnalysis, data.frame(siteID = siteID, SeWater))
subpop <- FixedSeDataAnalysis[, c(2, 9)]
colnames(subpop)[2] <- "ManUnit"
sizes <- with(FixedSeDataAnalysis, tapply(init.area, stratum.x, mean))
sizes <- sizes[c(1, 3, 7, 8, 9)]
WaterPopSizes.Adj <- list(ManUnit = as.list(sizes))
# using non area adjusted weights
Water.Est.Local.init <- cont.analysis(sites = WaterSites, subpop = subpop,
design = AreaAdj.Design, data.cont = WaterData, total = T, popsize = WaterPopSizes.Adj,
popcorrect = T, pcfsize = sum(sizes), vartype = "local", conf = 90,
pctval = c(0))
# adjusting using estimated area
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = adj.wgt.water, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
Water.Est.Local.adj <- cont.analysis(sites = WaterSites, subpop = subpop,
design = AreaAdj.Design, data.cont = WaterData, total = T, popsize = WaterPopSizes.Adj,
popcorrect = T, pcfsize = sum(sizes), vartype = "local", conf = 90,
pctval = c(0))
# ignoring unequal wts and assuming JUST A SRS within each stratum
# [WRONG!!]
Simple.SRS.water <- ddply(FixedSeDataAnalysis, ~stratum.x, summarize, SampleSize = length(SeWater),
Mean = mean(SeWater), SE = sd(SeWater)/sqrt(length(SeWater)))
maxes.water <- as.matrix(with(FixedSeDataAnalysis, tapply(SeWater, stratum.x,
max)))
####################################### BUGS ESTIMATES #############
FixedSeDataAnalysis <- read.csv("SEDataForAnalysisSpSurveyFinal.csv")
FixedSeDataAnalysis <- subset(FixedSeDataAnalysis, EvalBugs == "TS")
BugSites <- data.frame(siteID = FixedSeDataAnalysis$siteID, Active = rep(TRUE,
length(FixedSeDataAnalysis$SeBugs)))
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = non.area.adj.wgt.bugs, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
BugData <- with(FixedSeDataAnalysis, data.frame(siteID = siteID, SeBugs))
subpop <- FixedSeDataAnalysis[, c(2, 9)]
colnames(subpop)[2] <- "ManUnit"
26
sizes <- with(FixedSeDataAnalysis, tapply(init.area, stratum.x, mean))
sizes <- sizes[c(1, 3, 7, 8, 9)]
BugPopSizes.Adj <- list(ManUnit = as.list(sizes))
# using non area adjusted weights
Bug.Est.Local.init <- cont.analysis(sites = BugSites, subpop = subpop,
design = AreaAdj.Design, data.cont = BugData, total = T, popsize = BugPopSizes.Adj,
popcorrect = T, pcfsize = sum(sizes), vartype = "local", conf = 90,
pctval = c(0))
# adjusting using estimated area
AreaAdj.Design <- with(FixedSeDataAnalysis, data.frame(siteID = siteID,
wgt = adj.wgt.bugs, xcoord = xcoord, ycoord = ycoord, support = 0.008107))
Bug.Est.Local.adj <- cont.analysis(sites = BugSites, subpop = subpop, design = AreaAdj.Design,
data.cont = BugData, total = T, popsize = BugPopSizes.Adj, popcorrect = T,
pcfsize = sum(sizes), vartype = "local", conf = 90, pctval = c(0))
# ignoring unequal wts and assuming JUST A SRS within each stratum
# [WRONG!!]
Simple.SRS.bug <- ddply(FixedSeDataAnalysis, ~stratum.x, summarize, SampleSize = length(SeBugs),
Mean = mean(SeBugs), SE = sd(SeBugs)/sqrt(length(SeBugs)))
maxes.bug <- as.matrix(with(FixedSeDataAnalysis, tapply(SeBugs, stratum.x,
max)))
27



y <- read.csv("bnl_SeEggs_results_2013.csv")
