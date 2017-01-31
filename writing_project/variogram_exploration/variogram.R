# simulating variograms reference
# http://www.stat.purdue.edu/~zhanghao/ASS/homework/Project%201.pdf

require(geoR)
require(MASS) #mvrnorm
data("ca20")
ca <- ca20


# empirical variograms
ca_vc <- variog(ca, coords = ca$coords, data = ca$data, estimator.type = c("classical"))
ca_vm <- variog(ca, coords = ca$coords, data = ca$data, estimator.type = "modulus")
par(mfrow=c(1,2))
plot(ca_vc, main = "Classical Emprical")
# robust against outliers (Cressie's)
plot(ca_vm, main = "Modulus/Robust Empirical")

# kenny - issue with variogram is that use many points over again
# theoretically if outlier, showing up many times, for each pairwise distance
# so gotway suggests weighted gls to adjust, but still problem

# simulating from different correlation functions
linear_cor <- function(t,tau2,sigma2){
  I <- diag(1, nrow = dim(as.matrix(t))[1], ncol = dim(as.matrix(t))[2])
  nonspat_var <- tau2*I
  ifelse(t > 0,
  spat_var <- sigma2*as.matrix(t), spat_var <- 0)
  out <- nonspat_var + spat_var
  return(out)
}

linearV <- linear_cor(t = dist(ca$data), tau2 = 9, sigma2 = 0.1)

set.seed(25)
linearZ <- mvrnorm(mu = rep(0,dim(linearV)[1]), Sigma = linearV)

spherical_cor <- function(t,tau2,sigma2,phi){
  I <- diag(1, nrow = dim(as.matrix(t))[1], ncol = 
              dim(as.matrix(t))[2])
  nonspat_var <- tau2*I
spat_var <- ifelse(t<0, 0, ifelse(t >= (1/phi),
  sigma2,sigma2*(1.5*phi*as.matrix(t) - 
              0.5*(phi*as.matrix(t))^3)))
spat_var <- as.matrix(spat_var)
  out <- nonspat_var + spat_var
  return(out)
}

spherV <- spherical_cor(t = dist(ca$data), tau2 = 9, sigma2 = 0.1, phi = 2)

set.seed(25)
spherZ <- mvrnorm(mu = rep(0,dim(spherV)[1]), Sigma = spherV)
{\bf sphere not working}

gaussian_cor <- function(t,tau2,sigma2,phi){
  I <- diag(1, nrow = dim(as.matrix(t))[1], ncol = dim(as.matrix(t))[2])
  nonspat_var <- tau2*I
  spat_var <- ifelse(as.matrix(t) > 0,
                     sigma2*(1-exp(-phi^2*as.matrix(t)^2)),0)
  spat_var <- as.matrix(spat_var)
  out <- nonspat_var + spat_var
  return(out)
}

gaussV <- gaussian_cor(t = dist(ca$data), tau2 = 9, sigma2 = 0.1, phi = 2)

set.seed(25)
gaussZ <- mvrnorm(mu = rep(0,dim(gaussV)[1]), Sigma = gaussV)

exponential_cor <- function(t,tau2,sigma2,phi){
  I <- diag(1, nrow = dim(as.matrix(t))[1], ncol = dim(as.matrix(t))[2])
  nonspat_var <- tau2*I
  spat_var <- ifelse(t > 0,
                     sigma2*(1-exp(-phi*as.matrix(t))), 0)
  spat_var <- as.matrix(spat_var)
  out <- nonspat_var + spat_var
  return(out)
}



expV <- exponential_cor(t = dist(ca$data), tau2 = 9, sigma2 = 0.1, phi = 2)

set.seed(25)
expZ <- mvrnorm(mu = rep(0,dim(V)[1]), Sigma = expV)

##### Matern is not yet working

matern_cor <- function(t,tau2,sigma2,phi,nu){
  I <- diag(1, nrow = dim(as.matrix(t))[1], ncol = dim(as.matrix(t))[2])
  nonspat_var <- tau2*I
  coeff <- ((2*sqrt(nu)*as.matrix(t)*phi)^nu)/((2^(nu-1))*factorial(nu-1))
  bessel_idk <- besselI(2*sqrt(nu)*as.matrix(t)*phi,nu=nu)
  
  spat_var <- ifelse(t>0,
                     as.matrix(sigma2*(1-(coeff%*%bessel_idk)),
                               nrow = dim(as.matrix(t))[1], ncol = dim(as.matrix(t))[2]), 
                     as.matrix(0,nrow = dim(as.matrix(t))[1], ncol = dim(as.matrix(t))[2]))
  spat_var <- as.matrix(spat_var)
  out <- nonspat_var + spat_var
  return(out)
}

maternV <- matern_cor(t = dist(ca$data), tau2 = 9, sigma2 = 0.1, phi = 2, nu = 5)

set.seed(25)
maternZ <- mvrnorm(mu = rep(0,dim(V)[1]), Sigma = maternV)

## No cov
I <- diag(1, nrow = dim(as.matrix(dist(ca$data)))[1], ncol = dim(as.matrix(dist(ca$data)))[2])
Z <- mvrnorm(mu = rep(0,dim(as.matrix(dist(ca$data)))[1]), Sigma = 9*I)

## plot
par(mfrow=c(1,2))
plot(variog(data = linearZ, coords = ca$coords), main = "S2 = 0.1 Linear Cov")
plot(variog(data = spherZ, coords = ca$coords), main = "S2 = 0.1 Spherical Cov")
plot(variog(data = gaussZ, coords = ca$coords), main = "S2 = 0.1 Gaussian Cov")
plot(variog(data = expZ, coords = ca$coords), main = "S2 = 0.1 Exponential Cov")
plot(variog(data = maternZ, coords = ca$coords), main = "S2 = 0.1 Matern Cov")
plot(variog(data = Z, coords = ca$coords), main = "S2 = 0 No Cov")


# fractional sigma2, doesn't have much affect

# change sigma2
sigma2_vec <- c(0.9,1,2,8)

n <- length(ca$data)

linearV <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
linearV[,,i] <- linear_cor(t = dist(ca$data), tau2 = 9, sigma2 = sigma2_vec[i])
}

set.seed(25)
linearZ <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
linearZ[,,i] <- mvrnorm(mu = rep(0,dim(linearV)[1]), Sigma = linearV[,,i])
}

spherV <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
  spherV[,,i] <- spher_cor(t = dist(ca$data), tau2 = 9, sigma2 = sigma2_vec[i])
}

set.seed(25)
spherZ <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
  spherZ[,,i] <- mvrnorm(mu = rep(0,dim(spherV)[1]), Sigma = spherV[,,i])
}

gaussV <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
  gaussV[,,i] <- gauss_cor(t = dist(ca$data), tau2 = 9, sigma2 = sigma2_vec[i])
}

set.seed(25)
gaussZ <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
  gaussZ[,,i] <- mvrnorm(mu = rep(0,dim(gaussV)[1]), Sigma = gaussV[,,i])
}

expV <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
  expV[,,i] <- linear_cor(t = dist(ca$data), tau2 = 9, sigma2 = sigma2_vec[i])
}

set.seed(25)
expZ <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
  expZ[,,i] <- mvrnorm(mu = rep(0,dim(expV)[1]), Sigma = expV[,,i])
}

maternV <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
  maternV[,,i] <- linear_cor(t = dist(ca$data), tau2 = 9, sigma2 = sigma2_vec[i])
}

set.seed(25)
maternZ <- array(0,c(n,n,length(sigma2_vec)))
for(i in 1:length(sigma2_vec)){
  maternZ[,,i] <- mvrnorm(mu = rep(0,dim(maternV)[1]), Sigma = linearV[,,i])
}


for(i in 1:length(sigma2_vec)){
plot(variog(data = linearZ[,,i], coords = ca$coords),  
     main=substitute(paste("S2 =", a, "Linear Cov"), list(a=sigma2_vec[i])))
plot(variog(data = spherZ, coords = ca$coords),  
     main=substitute(paste("S2 =", a, "Spherical Cov"), list(a=sigma2_vec[i])))
plot(variog(data = gaussZ, coords = ca$coords),  
     main=substitute(paste("S2 =", a, "Gauss Cov"), list(a=sigma2_vec[i])))
plot(variog(data = expZ, coords = ca$coords),  
     main=substitute(paste("S2 =", a, "Exp Cov"), list(a=sigma2_vec[i])))
plot(variog(data = maternZ, coords = ca$coords),  
     main=substitute(paste("S2 =", a, "Matern Cov"), list(a=sigma2_vec[i])))
}





