setwd("~/Documents/toshiba122216/mack_hub/writing_project/data")
require(geoR)
#http://leg.ufpr.br/geoR/geoRdoc/geoRintro.html


se <- read.csv("SEDataForAnalysisSpSurveyFinal.csv")

head(se)

se.cols <- c("siteID", "xcoord", "ycoord", "mdcaty.x", "Date_", "Acres", "Hectares", 
             "SeSed", "SeRoots", "SeWater", "SeBugs")

se_sub <- se[,se.cols]
hist(se_sub$SeSed, breaks = c(seq(.2,21.21, by = 0.01)))
hist(log(se_sub$SeSed))

plot(se_sub$xcoord, se_sub$ycoord)

# fill in empty locations using bivariate linear interpolation
# it works... i don't know what spatial interpolation is
require(akima)
se_all <- se_sub[-c(which(is.na(se_sub$SeSed))),]
int.se <- interp(x = se_all$xcoord, y = se_all$ycoord, 
                 se_all$SeSed)
contour(sort(int.se$x), sort(int.se$y), int.se$z)

persp(int.se, xlim = range(se$xcoord), ylim = range(se$ycoord))

box.xbin <- c(seq(min(se_sub$xcoord) - 1, max(se_sub$xcoord), by = 
                   100000))
x.cut <- cut(se_sub$xcoord, box.xbin)
boxplot(se_sub$SeSed ~ x.cut)
