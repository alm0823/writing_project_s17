setwd("C:/Users/Andrea Mack/Desktop/mack_hub/writing_project/data")
setwd("~/Documents/toshiba122216/mack_hub/writing_project/data")
require(geoR)
#http://leg.ufpr.br/geoR/geoRdoc/geoRintro.html


se <- read.csv("SEDataForAnalysisSpSurveyFinal.csv")
head(se)

#The fields you are likely interested in for now are the xcoord, ycoord 
#(locations for the sample points); SiteID the unique id for a survey location; 

#mdcaty.x = refers to the water zones within a management unit: 
#flooded, saturated, intermittent; stratum.x = the unit name on the refuge so the different management units;

#Acres and Hectares are the size of the areal frame within each of those water zones within a unit; 
#Date_ = survey date when the data was collected;
#SeSed is selenium concentration within a sediment sample;
#SeRoot is selenium concentration within a root sample; 
#SeWater is selenium concentration within a water sample; and 
#SeBugs is selenium concentration within a bug sample.



se.cols <- c("siteID", "xcoord", "ycoord", "mdcaty.x", "Date_", "Acres", "Hectares", 
             "SeSed", "SeRoots", "SeWater", "SeBugs")

se_sub <- se[,se.cols]

#rep.data.action = "first" keeps only first coordinate pair if same coordinate
#pair represented more than once in dataframe
se_geo <- as.geodata(obj = se_sub, coords.col = c(2,3), 
                     data.col = c(8:11), covar.col = c(1,4,5,6,7), 
                     na.action = "none")

summary(se_geo)
#diff number of NAs in each group

plot(se_geo, col.data = 1)
#won't work bc of all the NAs

geo_all <- as.geodata(obj = se_sub, coords.col = c(2,3), 
                     data.col = c(8,9,10,11), 
                     covar.col = c(1,4,5,6,7), 
                     na.action = "none")
#I don't know what the colors mean... 
#documentation says they are the four "Quartiles"
#do they mean the four quantiles?
#quantiles of what?
#I am thinking of the coordinates?

## need to figure out how to do this with all data
plot(geo_all, na.rm = TRUE)
#showing higher se concentration 
#where the authors suggested, units 1 and 2 = low x coord and high y coord

sed_geo <- as.geodata(obj = se_sub, coords.col = c(2,3), 
                      data.col = c(8), 
                      covar.col = c(1,4,5,6,7), 
                      na.action = "ifdata")

plot(sed_geo)
points(sed_geo, pt.divide = "rank.prop")
#shows it more clearly

## continuing from
## http://www.leg.ufpr.br/geoR/geoRdoc/vignette/geoRintro/geoRintrose3.html
## eda plots


# variogram with classical estimator
# switch to modulus estimator by estimator.type = "modulus"
sed_cloud <- variog(sed_geo, option = "cloud")
sed_bin <- variog(sed_geo, uvec = "default")

par(mfrow=c(1,2))
plot(sed_cloud, main = "Cloud Variogram")
plot(sed_bin, main = "Binned Variogram")

# create a "binned" cloud variogram
sed_cloudbin <- variog(sed_geo, uvec = "default",
                       bin.cloud = TRUE)
plot(sed_cloudbin, bin.cloud = T, main = 
       "Classical Estimator, Binned Cloud Variogram")

## directional variograms, angles = 0,45,90,135
## tolerance at default of 22.5 degrees

## resource on interpretting variogram and quality measures

par(mfrow=c(1,1))
sed_vario4 <- variog4(sed_geo)
plot(sed_vario4)

## parameter estimation
## theoretical vs. empirical

plot(sed_bin)
lines.variomodel(cov.model = "exp", cov.pars = c(1,2),
                 nugget = 4)

sed_smooth <- variog(sed_geo, option = "smooth")#, ksmooth = "normal")
lines(sed_smooth, type = "l", lty = 2)
legend(0.4,0.3,c("empirical", "exponential model", 
                 "smoothed"), lty = c(1,1,2), lwd = c(1,3,1))


plot(variog(sed_geo, max.dist = 1)) 
lines.variomodel(cov.model = "exp", cov.pars = c(1, 
                                                   0.3), nug = 0, max.dist = 1) 
lines.variomodel(cov.model = "mat", cov.pars = c(0.85, 
                                                     0.2), nug = 0.1, kappa = 1, max.dist = 1, 
                     lty = 2) 
lines.variomodel(cov.model = "sph", cov.pars = c(0.8, 
                                                      0.8), nug = 0.1, max.dist = 1, lwd = 2)
#ml or reml parameter estimation for
#gaussian random fields
sed_ml <- likfit(sed_geo, ini.cov.pars = c(1,1))
summary(sed_ml)

# no nugget
options(geoR.messages = FALSE)
sed_ml <- likfit(sed_geo, ini = c(1,0.5),
                  fix.nugget = T)
sed_reml <- likfit(sed_geo, ini = c(1,0.5),
                   fix.nugget = T,
                   method = "RML")

sed_ols <- variofit(sed_bin, ini = c(1,0.5), 
                    fix.nugget = T, weights = "equal")

# fixed nugget
sed_ml.fn <- likfit(sed_geo,
                   ini=c(1,0.5), fix.nugget = T,
                  nugget = 0.15)
sed_reml.fn <- likfit(sed_geo, ini = c(1,0.5), 
                     fix.nugget = T,
                     nugget = 0.15,
                     weights = "equal")
sed_ols.fn <- variofit(sed_bin, ini = c(1,0.5),
                      fix.nugget = T,
                      nugget = 1, weights = "equal")
sed_wls.fn <- variofit(sed_bin, ini = c(1,2),
                       fix.nugget = T, nugget = 1)


require(akima)
contour(se$xcoord, se$ycoord, se$EvalSed)

#### Going to Jim's book ####
require(nlme)

sed_complete <- data.frame(se[which(is.na(se$SeSed) == FALSE), ])
sed.fit1 <- gls(SeSed ~ mdcaty.x -1,#no intercept
                data = sed_complete, na.action = "na.omit")#fit ols

plot(Variogram(sed.fit1, form = ~sed_complete$xcoord + sed_complete$ycoord))

#use random picks for range and nugget
sed.spher <- update(sed.fit1, corr = corSpher(c(50,0.1), form = ~ xcoord +
                                                ycoord, nugget = TRUE))


sed.ratio <- update(sed.fit1, corr = corRatio(c(50,0.1), form = ~ xcoord +
                                                ycoord, nugget = TRUE))

sed.lin <- update(sed.fit1, corr = corLin(c(50,0.1), form = ~ xcoord +
                                                ycoord, nugget = TRUE))


sed.exp <- update(sed.fit1, corr = corExp(c(50,0.1), form = ~ xcoord +
                                                ycoord, nugget = TRUE))

sed.gauss <- update(sed.fit1, corr = corGaus(c(50,0.1), form = ~ xcoord +
                                                ycoord, nugget = TRUE))

anova(sed.fit1, sed.spher, sed.ratio, sed.lin, sed.exp, sed.gauss, test = FALSE)


sed.fit1 <- gls(SeSed ~ mdcaty.x, data = sed_complete,#, corr = corSpher)#,form = ~ xcoord + ycoord, nugget = FALSE,
                na.action = "na.omit")#fit ols

########## section 4.2 kattegat salinity data #############
data("kattegat") #geodata object
# NAs are not allowed
summary(kattegat)
str(kattegat)

# from help examples
plot(c(550,770),c(6150,6420),type="n",xlab="X UTM",ylab="Y UTM")
points(kattegat, add=TRUE)
lapply(kattegat$dk, lines, lwd=2)

## hmm can't figure this out for se data 
minx_se <- min(se$xcoord)
maxx_se <- max(se$xcoord)
miny_se <- min(se$ycoord)
maxy_se <- max(se$ycoord)
plot(c(minx_se, maxx_se),c(minx_se,maxx_se),type="n",xlab="X Coord",ylab="Y Coord")
points(se_geo, add = TRUE)
lapply(se_geo$dk, lines, lwd = 2)

##############################################################
################### 4.2 trying to fit model ##################
##############################################################

# prior on phi
set.seed(53216)

phi <- c(seq(10, 100))

# find poster on phi -- discrete unif 10,100
# first need V,R,S

# need F == X matrix
# need R == correlation matrix, using exponential

names(kattegat)

require(stats)
katte_coords <- data.frame(cbind(kattegat$coords[,1], kattegat$coords[,2]))
fn_num <- function(x){
  as.numeric(as.character(x))
}

katte_coords <- apply(katte_coords, 2, fn_num)

# pairwise euclidean distances
u <- as.matrix(dist(katte_coords, method = "euclidean", p = 2))

fn_pu <- function(x,y){
  exp(-x/y)
}

R <- data.frame(matrix(vector(),70,70))
for(i in 1:70){
  for(j in 1:70){
    R[i,j] <- fn_pu(u[i,j],phi)
  }
}

# I belive x1 and x2 in the model are the lat and long coords
# but don' tknow how to check

F <- data.frame(cbind(c(rep(1,70)), katte_coords))
colnames(F) <- c("intercept", "x1", "x2")

# not sure if setting mb = 0 is correct.... check somehow
mb <- 0

# note that they set V_b = 0
R <- apply(R, 2, fn_num)
F <- apply(F, 2, fn_num)

Ft <- t(F)
dim(Ft)
R_inv <- solve(R)
dim(R_inv)



V_btilde <- solve(Ft %*% solve(R) %*% F)
dim(V_btilde)

# define the responses
y <- kattegat$data

Btilde <- V_btilde %*%(Ft%*% solve(R) %*% y)
dim(Btilde)
Btilde # not equal to authors Bhat p. 60

# let's say m_b is 0
n <- 70

S2 <- ((t(y) %*% R_inv %*% y) - (t(Btilde) %*% solve(V_btilde) %*% Btilde))/n
dim(S2)

# figure out whether n_sigma = -p or 0 here, to calc S2, we set = 0
# but discussion below 2.9 says -p

## I'm setting to 0 here

n_sigma <- 0
constant <- (det(V_btilde))^(1/2) * (det(R))^(-1/2) %*% S2^(-(n+n_sigma)/2)

phi.giveny <- phi*constant

## diffuse prior
## only need 1/sigma2

sigma.given. <- 1/rinvchisq(1,n_sigma + n, S2)


