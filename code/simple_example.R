library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(bayesplot::theme_default(base_family = "sans"))


source(here::here("monitornew.R"))

set.seed(666)
nreps = 1000
nchains=4;
niter = 1000;
npars=4
sd = matrix(data = 1.0,nrow=nchains, ncol= npars-2);
sd[1,1] = 1/3;
means = c(2,0,0,0)
ar_param = 0.3

Rhat = data.frame(rhat = rep(NA,2*npars*nreps), 
                  version = rep(NA,2*npars*nreps), 
                  par = rep(NA,2*npars*nreps),
                  neff = rep(NA,2*npars*nreps),
                  neff_bulk = rep(NA,2*npars*nreps),
                  neff_tail = rep(NA,2*npars*nreps))

sims = array(data=0.0,dim = c(niter,nchains,npars))
for (rep in 1:nreps) {
  for(chain in 1:nchains) {
    for (par in 1:2) {
      sims[,chain,par] = arima.sim(n=niter, list(ar=ar_param),sd=sqrt(1-ar_param^2))*sd[chain,par]
    }
    sims[,chain,3] = (arima.sim(n=niter, list(ar=ar_param),sd=sqrt(1-ar_param^2)))/arima.sim(n=niter, list(ar=ar_param),sd=sqrt(1-ar_param^2)) + means[chain]
    sims[,chain,4] = arima.sim(n=niter, list(ar=ar_param),sd=sqrt(1-ar_param^2))/arima.sim(n=niter, list(ar=ar_param),sd=sqrt(1-ar_param^2))
  }
  index_start =  2*npars*(rep-1) + 1;
  index_end = 2*npars*rep;
  nothing <- capture.output(monitor_old <- rstan::monitor(sims,warmup=0))
  nothing <- capture.output(monitor_new <- monitor(sims,warmup=0))
  par_names <- c("Finite mean, different variance", "Finite mean, same variance", "Infinite mean, different location","Infinite mean, same location")
  
  Rhat$rhat[index_start:index_end] = c(monitor_old[,10],monitor_new$Rhat)
  Rhat$version[index_start:index_end] = c(rep("old",npars),rep("new",npars))
  Rhat$par[index_start:index_end] = c(par_names,par_names)
  Rhat$neff[index_start:index_end] = c(monitor_old[,9],monitor_new$Bulk_ESS)
  Rhat$neff_bulk[index_start:index_end] = c(rep(NA,npars),monitor_new$Bulk_ESS)
  Rhat$neff_tail[index_start:index_end] = c(rep(NA,npars), monitor_new$Tail_ESS)
}

Rhat %>% ggplot(aes(rhat,colour=NULL, fill=version)) + 
  geom_histogram(position="identity", alpha=0.9) + facet_wrap(~par,nrow = 2,ncol = 2) +
  scale_fill_manual(name=paste("Version of Rhat"),
                      breaks=c("new", "old"),
                      labels=c("This paper","Gelman et al. (2013)"),
                      values=c("new" = "#d1e1ec", "old" ="#a25050")) +
  guides(colour=FALSE) + labs(y="", x="Rhat")+ labs(y="", x="Rhat")
ggsave(file="simple_rhat_compare.png",width=12, height = 7.5, units="in")
ggsave(file="../paper/graphics/simple_rhat_compare.png",width=12, height = 7.5, units="in")
Rhat %>% ggplot(aes(neff,colour=version, fill=version)) + geom_histogram() + facet_wrap(~par,nrow = 2,ncol = 2) + theme_minimal() + theme(strip.text = element_text(size=16),axis.title = element_text(size=16), axis.text = element_text(size=16)) 
ggsave(file="simple_neff_compare.png",width=12, height = 7.5, units="in")
Rhat %>% select(neff_bulk,neff_tail,par) %>% gather(`neff_bulk`,`neff_tail`,key="type",value="neff") %>% drop_na() %>%
  ggplot(aes(x=neff,colour=type,fill=type)) + geom_histogram() + facet_wrap(~par,nrow = 2,ncol = 2) + theme_minimal()+ theme(strip.text = element_text(size=16),axis.title = element_text(size=16), axis.text = element_text(size=16))
ggsave(file="simple_bulk_vs_tail.png",width=12, height = 7.5, units="in")

