

## 1.  Set up

setwd("/Users/andrew/AndrewFiles/research/rhat")
library("rstan")
library("rstanarm")
library("arm")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## 2.  Rescalilng

z_scale <- function(a){
  n <- length(a)
  r <- rank(a)
  u <- (2*r - 1)/(2*n)
  z <- qnorm(u)
  if (is.null(dim(a))) # assume it's a vector
    return(z)
  else # assume it's an array or matrix
    return(array(z, dim(a)))
}

## 3.  Rhats

rhats <- function(y){
  scales <- array(NA, c(M,6), dimnames=list(NULL, c("within_sd", "loo_sd", "e_hat", "within_scaled_sd", "loo_scaled_sd", "scaled_e_hat")))
  M <- ncol(y)
  y_scaled <- z_scale(y)
  for (m in 1:M){
    scales[m,1] <- sd(y[,m])
    scales[m,2] <- sd(as.vector(y[,-m]))
    scales[m,3] <- scales[m,1]/scales[m,2]
    scales[m,4] <- sd(y_scaled[,m])
    scales[m,5] <- sd(as.vector(y_scaled[,-m]))
    scales[m,6] <- scales[m,4]/scales[m,5]
  }
  total_sd <- sd(as.vector(y))
  rhat <- total_sd/sqrt(mean(scales[,1]^2))
  total_scaled_sd <- sd(as.vector(y_scaled))
  rhat_scaled <- total_scaled_sd/sqrt(mean(scales[,4]^2))
  return(list(rhat=rhat, rhat_scaled=rhat_scaled, scales=scales))
}

## 4.  pfround for list

lfround <- function(a, digits=0){
  if (is.list(a))
    return(lapply(a, function(x) fround(x, digits)))
  else
    fround(a, digits)
}

plfround <- function(a, digits) {
  print(lfround(a, digits), quote=FALSE)
}


## 5a.  Simulation:  normal with locations 0,0,0,0 and scales 1,1,1,1

N <- 1e5
M <- 4
y <- array(NA, c(N,M))
for (m in 1:M){
  y[,m] <- rnorm(N, 0, 1)
}
plfround(rhats(y), 2)

## 5b.  Simulation:  normal with locations 0,0,0,0 and scales 1,1,1,2

for (m in 1:M){
  y[,m] <- rnorm(N, 0, if (m==4) 2 else 1)
}
plfround(rhats(y), 2)

## 5c.  Simulation:  cauchy with locations 0,0,0,0 and scales 1,1,1,1

for (m in 1:M){
  y[,m] <- rcauchy(N, 0, 1)
}
plfround(rhats(y), 2)

## 5d.  Simulation:  cauchy with locations 0,0,0,0 and scales 1,1,1,2

for (m in 1:M){
  y[,m] <- rcauchy(N, 0, if (m==4) 10 else 1)
}
plfround(rhats(y), 2)

## 5e.  Simulation:  cauchy with locations 0,0,0,1 and scales 1,1,1,1

for (m in 1:M){
  y[,m] <- rcauchy(N, if (m==4) 1 else 0, 1)
}
plfround(rhats(y), 2)


