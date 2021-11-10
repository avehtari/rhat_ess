#' ---
#' title: "Comparison of MCMC effective sample size estimators"
#' author: "Aki Vehtari"
#' date: "First version 2021-11-03. Last modified `r format(Sys.Date())`."
#' output:
#'   html_document:
#'     theme: readable
#'     toc: true
#'     toc_depth: 3
#'     toc_float: true
#'     code_download: true
#' ---
#'
#' ## Post Script Appendix
#' 
#' This notebook is a post Script appendix to the paper <br>
#'
#' - Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter,
#' Paul-Christian Bürkner (2021): Rank-normalization, folding, and
#' localization: An improved $\widehat{R}$ for assessing convergence of MCMC. *Bayesian analysis*, 16(2):667-718. [doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221).
#'
#' The code used to make this notebook is available at [https://github.com/avehtari/rhat_ess](https://github.com/avehtari/rhat_ess).
#'
#' ## Introduction
#' 
#' The MCMC effective sample size (ESS) and Monte Carlo standard error
#' (MCSE) estimated for one chain includes estimation of the
#' correlation between the iterations, for example, using
#' autocorrelation time (or spectral density at frequency zero). As a
#' finite number of MCMC draws are available the autocorrelation
#' estimates for bigger lags are noisier.  In the paper, we wrote
#' about autocorrelation estimators
#'
#' > In our experiments, Geyer's [(1992)] truncation had superior
#' > stability compared to flat-top (Doss et al., 2014) and slug-sail
#' > (Vats and Knudson, 2018) lag window approaches.
#'
#' (Vats and Knudson (2018) has since been published as Vats and
#' Knudson (2021)).  As the main point of the paper was not comparison
#' of autocorrelation estimators, the experimental results were not
#' included in the paper.
#'
#' This notebook includes additional experiments comparing the Geyer's
#' (1992) truncation approach to three other methods for estimating
#' effective sample size (ESS) and corresponding Monte Carlo standard
#' error (MCSE). We also compare the new quantile MCSE method we
#' proposed in the paper (Vehtari et al., 2021) to batch means
#' quantile MCSE method by (Gong and Flegal, 2015). Finally, we
#' compare the behavior of ESS implementations in four R packages in
#' case of multiple chains and well separated multimodal distribution.
#'
#' ### The methods and R packages compared
#'
#' - [`stableGR`](https://cran.r-project.org/package=stableGR) (`n.eff`) computes a weighted sum of autocorrelations using slug-sail lag window (Vats and Knudson, 2018) 
#' - ['mcmcse'](https://cran.r-project.org/package=mcmcse) (`ess` and `mcse.q`) uses by default a batch means approach (Gong and Flegal, 2015) 
#' - [`coda`](https://cran.r-project.org/package=coda) (`effectiveSize`) computes the spectral density at frequency zero by fitting an AR model to the chain (Heidelberger and Welch, 1981)
#' - [`posterior`](https://cran.r-project.org/package=posterior): `ess_basic` computes sum of autocorrelations using Geyer's truncation rule (Geyer, 1992), and `mcse_quantile` uses in addition an inverse approach for quantile MCSE. The `posterior` package also takes into account the the possibility that the chains are not mixing by using also multi-chain $\widehat{R}$ in ESS and MCSE computations.
#'
#' ### tl;dr
#'
#' In case of well mixing chains, `ess_basic` and `mcse_quantile` in
#' the `posterior` package are the most accurate methods. In case of
#' well separated modes, the `posterior` package also provides
#' appropriately small ESS indicating that the chains are not mixing.
#' 
#+ setup, include=FALSE
knitr::opts_chunk$set(message=FALSE, error=FALSE, warning=FALSE, comment=NA, cache=FALSE)
#+ load_packages, echo=FALSE
library(tidyverse)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(khroma)
library(latex2exp)
#'
#' ## Comparison of ESS estimators with AR(1) simulation
#'
#' ### Mean of x
#' 
#' We simulate $M=4$ chains from a known AR(1) process with varying
#' parameter $\phi$, and normal(0,1) marginal distribution.
#' ```
#' x1 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
#' x2 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
#' x3 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
#' x4 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)));
#' ```
#'
#' We vary the length of the chains, $N$, and the AR process parameter $\phi$
#' ```
#' Ns=c(100, 1000, 10000)
#' phis=c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)
#' ```
#'
#' $N=1000$ is the default in Stan, and $N=100$ and $N=10000$ are used
#' to illustrate the behavior with less or more MCMC iterations. The
#' total number of draws is $S=M*N$.
#'
#' The values of $\phi$ have been chosen to reflect typical ESS/$S$
#' values seen from Stan. The chosen $\phi$ values correspond to ESS/$S$
#' values $(1.9, 1.2, 0.8, 0.5, 0.3, 0.2, 0.1)$. The two first $\phi$
#' values being negative produce antithetic chains with ESS>$S$, which
#' is not rare when using Stan's dynamic Hamiltonian Monte Carlo.
#'
#' Simulation is repeated $100\,000$ times for each combination of $N$
#' and $\phi$. For each simulation we compute the empirical mean. The
#' variance of the empirical means and known true variance are
#' used to compute the asymptotic efficiency, that is, true ESS/$S$ when
#' $S\rightarrow \infty$. For each simulation we compute three ESS
#' estimates and use those to compute Monte Carlo standard error
#' (MCSE) which is compared to the actual error, that is, the
#' difference between the empirical mean $\hat{\mu}$ and known true
#' mean $\mu=0$.
#'
#' The following figure summarizes the results.  In the case of
#' estimating some posterior expectation $\mathrm{E}[g(x)]$, the used
#' ESS estimators assume finite variance and MCSE is simply
#' $\mathrm{sd}[g(x)]/\sqrt{\mathrm{ESS}}$. Thus accuracy of ESS is
#' directly related to the accuracy of MCSE. We report root mean
#' square of standardized errors $|\mu -
#' \hat{\mu}|/\widehat{\mathrm{MCSE}}$, which is in case of well calibrated ESS
#' and MCSE should be close to 1. Values less than 1 indicate
#' underestimated ESS and overestimated MCSE, and values larger than 1
#' indicate overestimated ESS and underestimated MCSE.
#' 
#+ plot1, echo=FALSE, fig.dim=c(9,3)
load('ess_ar_df.Rdata')
Ns=c(100, 1000, 10000)
phis=c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)
ggplot(df, aes(x=phi, y=z, 
               group=method,
               color=method)) +
  geom_line(alpha=0.8, size=1) +
  geom_hline(yintercept = 1, linetype='dashed') +
  coord_cartesian(ylim=c(0.91,1.3), clip='off')+
  labs(x='True ESS/S', y=TeX('RMS of $|\\mu -\\hat{\\mu}| / \\widehat{MCSE}$'))+
  guides(linetype = guide_legend(title = "N"),
         color = guide_legend(title = "Package: Method"))+
  scale_x_continuous(breaks=phis, labels = round(tess/40000,1))+
  scale_color_bright(breaks=c("stableGR::n.eff","mcmcse::ess","coda::effectiveSize","posterior::ess_basic"), labels=c('stableGR: Slug-sail','mcmcse: Batch means','coda: AR','posterior: Geyer'))+
  facet_grid(cols=vars(S=S))
#'
#' `posterior: Geyer` is on average clearly closer to 1 than the other
#' methods. `coda: AR` method is overconfident with short chains, but
#' has good behavior with longer chains. `stableGR: Slug-sail` has
#' quite unpredictable behavior being either over- or underconfident,
#' except in case of very long chains when the method is consistently
#' overconfident for antithetic chains and underconfident for
#' "thetic" chains. Not shown in the plot, but slug-sail estimator
#' had overall much higher variance than Geyer's and AR approaches.
#'
#' Although the differences between methods are clear, the effect of
#' using the methods with bigger errors is likely to be small in most
#' actual Bayesian modeling use cases. Except for the shortest chains
#' with the highest autocorrelations, the errors in MCSE are at
#' highest 20\%, which should provide sufficient information whether
#' more draws are needed. With more draws the accuracy of all methods
#' improve.
#'
#' ### Mean of x^2
#'
#' In the above the chains were simulated from AR(1) process, so we would
#' expect that `coda: AR` method would have some advantage. Also the
#' chains have normal(0,1) marginal which makes the empirical means
#' also to behave nicely. As an alternative we test the methods also
#' with the same AR processes, but values are squared so that the
#' chains are non-linear transformation of the AR process and the
#' marginal is skewed. Squared variable is, for example, a natural part
#' of variance computation.
#'```
#' x1 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))^2;
#' x2 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))^2;
#' x3 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))^2;
#' x4 = arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))^2;
#' ```
#'
#' The chosen $\phi$ values correspond to now ESS/$S$ values $(0.8, 1.0,
#' 1.0, 0.8, 0.6, 0.3, 0.1)$. The squaring breaks the antithetic
#' behavior, but for simplicity we show the results for all different
#' values of $\phi$.
#'
#+ plot2, echo=FALSE, fig.dim=c(9,3)
load('ess_ar2_df.Rdata')
ggplot(df, aes(x=phi, y=z, 
               group=method,
               color=method)) +
  geom_line(alpha=0.8, size=1) +
  geom_hline(yintercept = 1, linetype='dashed') +
  coord_cartesian(ylim=c(0.91,1.4), clip='off')+
  labs(x='True ESS/S', y=TeX('RMS of $|\\mu -\\hat{\\mu}| / \\widehat{MCSE}$'))+
  guides(linetype = guide_legend(title = "N"),
         color = guide_legend(title = "Method"))+
  scale_x_continuous(breaks=phis, labels = round(tess/40000,1))+
  scale_color_bright(breaks=c("stableGR::n.eff","mcmcse::ess","coda::effectiveSize","posterior::ess_basic"), labels=c('stableGR: Slug-sail','mcmcse: Batch means','coda: AR','posterior: Geyer'))+
  facet_grid(cols=vars(S=S))
#'
#' Again `posterior: Geyer` is on average clearly closer to 1 than the
#' other methods. All methods have bigger errors for the short chains
#' (S=4 x 100), especially when ESS/$S$ is small, which is natural as
#' the draws are from a skewed distribution and the variability in
#' means and correlations increases. With increasing number of draws
#' the error decreases for all methods.
#' 
#' ### Probability p(x<0)
#'
#' Next we examine the behavior for estimating probabilities. The
#' sequence is now binary with values 0 and 1, computed with an
#' indicator function I(x<0).
#'```
#' x1 = as.double(arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))<0);
#' x2 = as.double(arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))<0);
#' x3 = as.double(arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))<0);
#' x4 = as.double(arima.sim(n = N, list(ar = c(phi)), sd = sqrt((1-phi^2)))<0);
#' ```
#'
#' The chosen $\phi$ values correspond to now ESS/$S$ values $(1.4,
#' 1.1, 0.9, 0.6, 0.4, 0.2, 0.1)$. We see super-efficiency for
#' antithetic chains.
#' 
#+ plot3, echo=FALSE, fig.dim=c(9,3)
load('ess_arp_df.Rdata')
ggplot(df, aes(x=phi, y=z, 
               group=method,
               color=method)) +
  geom_line(alpha=0.8, size=1) +
  geom_hline(yintercept = 1, linetype='dashed') +
  coord_cartesian(ylim=c(0.91,1.4), clip='off')+
  labs(x='True ESS/S', y=TeX('RMS of $|\\mu -\\hat{\\mu}| / \\widehat{MCSE}$'))+
  guides(linetype = guide_legend(title = "N"),
         color = guide_legend(title = "Method"))+
  scale_x_continuous(breaks=phis, labels = round(tess/40000,1))+
  scale_color_bright(breaks=c("stableGR::n.eff","mcmcse::ess","coda::effectiveSize","posterior::ess_basic"), labels=c('stableGR: Slug-sail','mcmcse: Batch means','coda: AR','posterior: Geyer'))+
  facet_grid(cols=vars(S=S))
#'
#' Again `posterior: Geyer` is on average clearly closer to 1 than the
#' other methods. All methods have bigger errors for the short chains
#' (S=4 x 100), especially when ESS/$S$ is small, which is natural as
#' binary observation contain less information than continuous
#' observations.  With increasing number of draws the error decreases
#' for all methods.
#'
#' ### Median of x
#'
#' Quantiles are not expectations, but can be easily estimated from
#' MCMC draws. Computing MCSE requires additional steps. `stableGR`
#' and `coda` packages don't have a function for computing MCSE for
#' quantiles (or at least it is not mentioned in the
#' documentation). `mcmcse` package function `mcse.q` uses batch means
#' approach. `posterior` package function `mcse_quantile` uses the
#' inversion approach presented by Vehtari et al. (2021).
#' 
#+ plot4, echo=FALSE, fig.dim=c(9,3)
load('ess_arq_df.Rdata')
clr<-colour("bright",names=FALSE)(7)
ggplot(df, aes(x=phi, y=z, 
               group=method,
               color=method)) +
  geom_line(alpha=0.8, size=1) +
  geom_hline(yintercept = 1, linetype='dashed') +
  coord_cartesian(ylim=c(0.95,1.25), clip='off')+
  labs(x='True ESS/S', y=TeX('RMS of $|\\mu -\\hat{\\mu}| / \\widehat{MCSE}$'))+
  guides(linetype = guide_legend(title = "N"),
         color = guide_legend(title = "Method"))+
  scale_x_continuous(breaks=phis, labels = round(tess/40000,1))+
  scale_color_manual(breaks=c("mcmcse::mcse.q","posterior::mcse_quantile"), labels=c('mcmcse: Batch means','posterior: Geyer + inv.prob.'), values=c(clr[2],clr[4]))+
  facet_grid(cols=vars(S=S))
#'
#' The plots for `posterior: Geyer` and `mcmcse: batch means` look
#' similar as in the probability estimation which is expected due to
#' the close connection between computing cumulative probability and
#' quantile value.  `posterior: Geyer` line is on average closer to 1.
#' 
#' ## Comparison of ESS estimators for bimodal distribution
#'
#' Above, AR process was used to simulate MCMC sampling from a
#' unimodal distribution. We advocate running many chains not just for
#' the parallelization of the computation, but also to be able to
#' detect possible multimodality. Here we sample from a bimodal
#' distribution using Stan. We vary the number of chains and examine
#' the behavior of the ESS estimators in the same three packages.
#'
#' `posterior` package is the only one of these packages taking into
#' account the the possibility that the chains are not mixing by using
#' also multi-chain \widehat{R} in ESS and MCSE
#' computations. `stableGR` and `coda` accept multiple chains as
#' input, but don't show in ESS estimates indication of between chain
#' mixing problems. `mcmcse` accepts only one chain as input, and
#' multiple chains were simply concatenated to provide a "single"
#' chain as input.
#'
#' The model is a simple Student's $t_4$ distribution with unknown
#' location and scale.
#' ```
#' data {int<lower=0> N; vector[N] y }
#' parameters {real mu;}
#' model {mu ~ student_t(4, 0, 100); y ~ student_t(4, mu, 1);}
#' ```
#' The data is generated from a bimodal distribution.
#' ```
#' N=20; y=c(rnorm(N/2, mean=-5, sd=1),rnorm(N/2, mean=5, sd=1));
#' ```
#'
#' With this data, the posterior is bimodal with well separated
#' modes. The modes don't necessarily have equal posterior
#' mass. Stan's dynamic HMC is not able to jump from one mode to
#' another. The random initialization determines to which mode each
#' chain ends up, and it is not possible to infer the actual relative
#' posterior masses from the chains. The number of chains in each mode
#' is random, and the probability of a chain with random
#' initialization ending to a mode depends on the volume of the
#' attraction area of that mode.  In such case we would like to get a
#' warning that the chains are not mixing.
#'
#' The following plot shows an example of trace plot when one chain
#' is stuck in one mode and three in another mode.
#+ echo=FALSE, fig.dim=c(9,3)
load('ess_df2.Rdata')
mcmc_trace(draws)
#'
#' If we would repeat the sampling with random initializations and four
#' chains, we would observe 0--4 chains stuck in one mode and the rest
#' in the other mode. So there is still a possibility that all the
#' chains are stuck in one mode and we think that the posterior is
#' unimodal. With just one chain, we have probability of 0 recognizing
#' that there are well separated modes. With more chains, it is
#' possible that our random initialization is within the attraction
#' area of just one mode, and then having more chains doesn't help. In
#' general, with increasing number of chains, we have higher
#' probability of seeing more than one mode.
#'
#' How should ESS estimates look like when we have chains that are not
#' mixing? If the modes are well separated so that between-chain
#' variance is much bigger than within-chain variances, each chain is
#' mostly indicating the loaction of the mode and not worth much more
#' than one observation.  The `posterior` package combines Geyer
#' truncated autocorrelation estimates with multi-chain convergence
#' diagnostic $\widehat{R}$. If the chains are not mixing, then the
#' draws within each chain are considered to be more correlated within
#' the chain and the multi-chain ESS estimate gets smaller. If the
#' modes are well separated, multi-chain ESS is close to the number of
#' chains. Arguably even this is an overestimate as increasing the
#' number of non-mixing chains does not provide reliable estimate of
#' the relative masses of the masses.
#'
#' The following plot shows ESS per chain estimates using the four
#' packages and with increasing number of chains. Simulations have
#' been repeated 100 times.
#' 
#+ echo=FALSE, fig.dim=c(9,4)
load('ess_df2.Rdata')
ggplot(df2, aes(x=log2(chains),
                y=ess/chains, color=factor(method))) +
  geom_jitter(width=0.07,height=0,alpha=1,shape=1) +
  geom_hline(yintercept = 1, linetype='dashed') +
  labs(x='Number of chains', y='Estimated ESS per chain')+
  guides(color = "none")+
  scale_x_continuous(breaks=0:4, labels = c(1,2,4,8,16))+
  scale_y_continuous(trans='log10')+
  scale_color_bright(breaks=c("1","2","3","4"), labels=c('stableGR','mcmcse','coda','posterior'))+
  facet_grid(cols=vars(methodlab))
#'
#' With one chain, all packages consider only the within-chain
#' correlation and report similar ESS per chain values (~350). When
#' different chains end up to stuck in different well-separated modes,
#' we would prefer to get an estimate that is indicating poor
#' efficiency. With increasing number of chains, `coda` keeps
#' reporting similar ESS per chain (with less variation) and provides
#' no indication of non-mixing chains.  When at least on chain is
#' stuck in a different mode, `stableGR` reports ESS per chain of
#' approximately 60, which is much less than 350, but we could see
#' similar low values due to high autocorrelation and the total ESS is
#' high enough not to indicate that chains are not mixing
#' well. `mcmcse` does not have support for multiple chains, but
#' combining the chains to one chain, changes the autocorrelation
#' estimate so that when the number of chains increases, eventually
#' ESS per chain decreases close to 1, indicating inefficient
#' sampling.  `posterior` reports ESS per chain close to 1 always when
#' at least one chain is in the other mode, which is a clear
#' indication that something is wrong in the mixing of the chains.
#' 
#' ## References
#' 
#' - Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statistical Science*, 7:473–483.
#' - Gong, L., and Flegal, J. M. (2015). A practical sequential stopping rule for high-dimensional Markov chain Monte Carlo. *Journal of Computational and Graphical Statistics*, 25:684—700.
#' - Heidelberger, P., and Welch, P. D. (1981). A spectral method for confidence interval generation and run length control in simulations. *Communications of the ACM*, 24:233—245.
#' - Vats, D., and Knudson, C. (2018). Revisiting the Gelman-Rubin diagnostic. *arXiv preprint*, arXiv:1812.09384.
#' - Vats, D., and Knudson, C. (2021). Revisiting the Gelman-Rubin diagnostic. *Statistical Science*, 36(4): 518—529.
#' - Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., Bürkner, P.-C. (2021): Rank-normalization, folding, and localization: An improved $\widehat{R}$ for assessing convergence of MCMC. *Bayesian analysis*, 16(2):667—718. [doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221).
#'
