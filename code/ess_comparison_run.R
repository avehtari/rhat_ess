library("rprojroot")
root<-has_file(".root")$make_fix_file()
library(tidyverse)
library(cmdstanr)
library(stableGR)
library(coda)
library(mcmcse)
library(posterior)

#' Simulate 4 chains from AR process, with varying length and phi
Ns=c(100, 1000, 10000)
phis=c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)
mumean=array(dim=c(length(Ns),length(phis),100000))
muvar=array(dim=c(length(Ns),length(phis),100000))
ess1=array(dim=c(length(Ns),length(phis),100000))
ess2=array(dim=c(length(Ns),length(phis),100000))
ess3=array(dim=c(length(Ns),length(phis),100000))
ess4=array(dim=c(length(Ns),length(phis),100000))
for (n in 1:length(Ns)) {
  N=Ns[n];
  print(N);
  for (j in 1:length(phis)) {
    print(j)
    phi=phis[j];
    for (i in 1:100000) {
      set.seed(i)
      x1 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
      x2 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
      x3 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
      x4 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
      mumean[n,j,i]=mean(c(x1,x2,x3,x4))
      muvar[n,j,i]=var(c(x1,x2,x3,x4))
      ess1[n,j,i]=stableGR::n.eff(list(matrix(x1,ncol=1),matrix(x2,ncol=1),matrix(x3,ncol=1),matrix(x4,ncol=1)))$n.eff/4
      ess2[n,j,i]=posterior::ess_basic(matrix(c(x1,x2,x3,x4),ncol=4))
      ess3[n,j,i]=coda::effectiveSize(mcmc.list(mcmc(x1),mcmc(x2),mcmc(x3),mcmc(x4)))
      capture.output(ess4[n,j,i]<-ess(c(x1,x2,x3,x4)))
    }
  }
}
for (j in 1:length(phis)) {
  (tess[j]<-1/var(mumean[3,j,]))
}

df=tibble()
i=0;
for (n in 1:length(Ns)) {
  N=Ns[n];
  for (j in 1:length(phis)) {
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='stableGR::n.eff'
    df[i,'ratio']=mean(ess1[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess1[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess1[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((0-mumean[n,j,])/sqrt(muvar[n,j,]/ess1[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='posterior::ess_basic'
    df[i,'ratio']=mean(ess2[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess2[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess2[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((0-mumean[n,j,])/sqrt(muvar[n,j,]/ess2[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='coda::effectiveSize'
    df[i,'ratio']=mean(ess3[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess3[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess3[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((0-mumean[n,j,])/sqrt(muvar[n,j,]/ess3[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='mcmcse::ess'
    df[i,'ratio']=mean(ess4[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess4[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess4[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((0-mumean[n,j,])/sqrt(muvar[n,j,]/ess4[n,j,]))^2))
  }
}
df$S = factor(df$N, labels = c('S = 4 x 100','S = 4 x 1000','S = 4 x 10000'))
save(file='ess_ar_df.Rdata','df','tess')

#' Simulate 4 chains from AR process, with varying length and
#' phi. This time square the values.
mumean=array(dim=c(length(Ns),length(phis),100000))
muvar=array(dim=c(length(Ns),length(phis),100000))
ess1=array(dim=c(length(Ns),length(phis),100000))
ess2=array(dim=c(length(Ns),length(phis),100000))
ess3=array(dim=c(length(Ns),length(phis),100000))
ess4=array(dim=c(length(Ns),length(phis),100000))
for (n in 1:length(Ns)) {
  N=Ns[n];
  print(N);
  for (j in 1:length(phis)) {
    print(j)
    phi=phis[j];
    for (i in 1:100000) {
      set.seed(i)
      x1 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))^2;
      x2 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))^2;
      x3 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))^2;
      x4 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))^2;
      mumean[n,j,i]=mean(c(x1,x2,x3,x4))
      muvar[n,j,i]=var(c(x1,x2,x3,x4))
      ess1[n,j,i]=stableGR::n.eff(list(matrix(x1,ncol=1),matrix(x2,ncol=1),matrix(x3,ncol=1),matrix(x4,ncol=1)))$n.eff/4
      ess2[n,j,i]=posterior::ess_basic(matrix(c(x1,x2,x3,x4),ncol=4))
      ess3[n,j,i]=coda::effectiveSize(mcmc.list(mcmc(x1),mcmc(x2),mcmc(x3),mcmc(x4)))
      capture.output(ess4[n,j,i]<-ess(c(x1,x2,x3,x4)))
    }
  }
}
for (j in 1:length(phis)) {
  (tess[j]<-2/var(mumean[3,j,]))
}

df=tibble()
i=0;
for (n in 1:length(Ns)) {
  N=Ns[n];
  for (j in 1:length(phis)) {
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='stableGR::n.eff'
    df[i,'ratio']=mean(ess1[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess1[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess1[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((1-mumean[n,j,])/sqrt(muvar[n,j,]/ess1[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='posterior::ess_basic'
    df[i,'ratio']=mean(ess2[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess2[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess2[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((1-mumean[n,j,])/sqrt(muvar[n,j,]/ess2[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='coda::effectiveSize'
    df[i,'ratio']=mean(ess3[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess3[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess3[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((1-mumean[n,j,])/sqrt(muvar[n,j,]/ess3[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='mcmcse::ess'
    df[i,'ratio']=mean(ess4[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess4[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess4[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((1-mumean[n,j,])/sqrt(muvar[n,j,]/ess4[n,j,]))^2))
  }
}
df$S = factor(df$N, labels = c('S = 4 x 100','S = 4 x 1000','S = 4 x 10000'))
save(file='ess_ar2_df.Rdata','df','tess')

#' Simulate 4 chains from AR process, with varying length and
#' phi. This time estimate p(x<0)
mumean=array(dim=c(length(Ns),length(phis),100000))
muvar=array(dim=c(length(Ns),length(phis),100000))
ess1=array(dim=c(length(Ns),length(phis),100000))
ess2=array(dim=c(length(Ns),length(phis),100000))
ess3=array(dim=c(length(Ns),length(phis),100000))
ess4=array(dim=c(length(Ns),length(phis),100000))
for (n in 1:length(Ns)) {
  N=Ns[n];
  print(N);
  for (j in 1:length(phis)) {
    print(j)
    phi=phis[j];
    for (i in 1:100000) {
      set.seed(i)
      x1 = (arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))<=0)+0;
      x2 = (arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))<=0)+0;
      x3 = (arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))<=0)+0;
      x4 = (arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))<=0)+0;
      mumean[n,j,i]=mean(c(x1,x2,x3,x4))
      muvar[n,j,i]=var(c(x1,x2,x3,x4))
      ## capture.output(ess1[n,j,i]<-stableGR::n.eff(list(matrix(x1,ncol=1),matrix(x2,ncol=1),matrix(x3,ncol=1),matrix(x4,ncol=1)))$n.eff/4)
      ## ess2[n,j,i]=posterior::ess_basic(matrix(c(x1,x2,x3,x4),ncol=4))
      ## capture.output(ess3[n,j,i]<-coda::effectiveSize(mcmc.list(mcmc(x1),mcmc(x2),mcmc(x3),mcmc(x4))))
      capture.output(ess4[n,j,i]<-ess(c(x1,x2,x3,x4)))
    }
  }
}

for (j in 1:length(phis)) {
  (tess[j]<-1/4/var(mumean[3,j,]))
}

df=tibble()
i=0;
for (n in 1:length(Ns)) {
  N=Ns[n];
  for (j in 1:length(phis)) {
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='stableGR::n.eff'
    df[i,'ratio']=mean(ess1[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess1[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess1[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((0.5-mumean[n,j,])/sqrt(muvar[n,j,]/ess1[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='posterior::ess_basic'
    df[i,'ratio']=mean(ess2[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess2[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess2[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((0.5-mumean[n,j,])/sqrt(muvar[n,j,]/ess2[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='coda::effectiveSize'
    df[i,'ratio']=mean(ess3[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess3[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess3[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((0.5-mumean[n,j,])/sqrt(muvar[n,j,]/ess3[n,j,]))^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='mcmcse::ess'
    df[i,'ratio']=mean(ess4[n,j,]*10^(3-n)/tess[j])
    df[i,'sd']=sd(ess4[n,j,]*10^(3-n)/tess[j])
    df[i,'rmse']=sqrt(mean((1-ess4[n,j,]*10^(3-n)/tess[j])^2))
    df[i,'z']=sqrt(mean(abs((0.5-mumean[n,j,])/sqrt(muvar[n,j,]/ess4[n,j,]))^2))
  }
}
df$S = factor(df$N, labels = c('S = 4 x 100','S = 4 x 1000','S = 4 x 10000'))
save(file='ess_arp_df.Rdata','df','tess')


#' Simulate 4 chains from AR process, with varying length and
#' phi. This time estimate median.
mumean=array(dim=c(length(Ns),length(phis),100000))
muvar=array(dim=c(length(Ns),length(phis),100000))
mcse1=array(dim=c(length(Ns),length(phis),100000))
mcse2=array(dim=c(length(Ns),length(phis),100000))
mcse3=array(dim=c(length(Ns),length(phis),100000))
mcse4=array(dim=c(length(Ns),length(phis),100000))
for (n in 1:length(Ns)) {
  N=Ns[n];
  print(N);
  for (j in 1:length(phis)) {
    print(j)
    phi=phis[j];
    for (i in 1:100000) {
      set.seed(i)
      x1 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
      x2 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
      x3 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
      x4 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
      mumean[n,j,i]=median(c(x1,x2,x3,x4))
      muvar[n,j,i]=var(c(x1,x2,x3,x4))
      mcse2[n,j,i]=posterior::mcse_quantile(matrix(c(x1,x2,x3,x4),ncol=4),0.5)
      capture.output(mcse4[n,j,i]<-mcse.q(c(x1,x2,x3,x4),0.5)$se)
    }
  }
}

for (j in 1:length(phis)) {
  (tess[j]<-(pi/2)/var(mumean[3,j,]))
}

df=tibble()
i=0;
for (n in 1:length(Ns)) {
  N=Ns[n];
  for (j in 1:length(phis)) {
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='posterior::mcse_quantile'
    df[i,'z']=sqrt(mean(abs((0-mumean[n,j,])/mcse2[n,j,])^2))
    i=i+1
    df[i,'phi']=phis[j]
    df[i,'N']=Ns[n]
    df[i,'method']='mcmcse::mcse.q'
    df[i,'z']=sqrt(mean(abs((0-mumean[n,j,])/mcse4[n,j,])^2))
  }
}
df$S = factor(df$N, labels = c('S = 4 x 100','S = 4 x 1000','S = 4 x 10000'))
save(file='ess_arq_df.Rdata','df','tess')


#' Sample from a bimodal distribution with Stan. Vary the number of
#' chains from 1 to 16.
N=20
set.seed(747)
y=c(rnorm(N/2, mean=-5, sd=1),rnorm(N/2, mean=5, sd=1));
data_tt <-list(N = N, y = y)
#' Student's t model
code_tt <- root("ESScomparison", "student.stan")
#' Sample
#+ message=TRUE, error=TRUE, warning=TRUE
mod_tt <- cmdstan_model(stan_file = code_tt)

chainss=c(1,2,4,8,16)
mumean= muvar= ess1= ess2= ess3= array(dim=c(length(chainss), 100))
df2 = tibble()
ii=0
for (j in 1:length(chainss)) {
  print(j)
  for (i in 1:100) {
    capture.output(fit_tt <- mod_tt$sample(data = data_tt, seed = i,
                                           refresh = 0, chains = chainss[j]))
    draws <- as_draws_rvars(fit_tt$draws(variables="mu"))
    ii=ii+1
    df2[ii,'method']=4;
    df2[ii,'chains']=chainss[j]
    df2[ii,'ess']=posterior::ess_basic(draws$mu)
    ii=ii+1
     draws <- as_draws_array(draws)
    df2[ii,'method']=1;
    df2[ii,'chains']=chainss[j]
    df2[ii,'ess']=stableGR::n.eff(lapply(1:chainss[j],function(c) { matrix(draws[,c,1]) }))$n.eff/chainss[j]
    ii=ii+1
    df2[ii,'method']=3;
    df2[ii,'chains']=chainss[j]
    df2[ii,'ess']=coda::effectiveSize(mcmc.list(lapply(1:chainss[j],function(c) { mcmc(matrix(draws[,c,1])) })))
    ii=ii+1
    df2[ii,'method']=2;
    df2[ii,'chains']=chainss[j]
    df2[ii,'ess']=mcmcse::ess(as.vector(draws))
  }
}
df2$methodlab=factor(df2$method, labels = c('stableGR','mcmcse','coda','posterior'))
fit_tt <- mod_tt$sample(data = data_tt, seed = 1, refresh = 0)
draws <- fit_tt$draws(variables="mu")
save(file='ess_df2.Rdata','df2','draws')
