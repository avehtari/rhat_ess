Repo for split-R-hat and effective sample size estimation.

Draft is based on material from Stan manual and BDA3.

Description of effective sample size computation has been fixed
 - autocorrelations are computed with FFT instead of covariogram
 - effective sample size `N_eff > N`

New ideas
 - focus 1. Rhat, 2. n_eff, 3. MCSE
 - separate Rhat diagnostic and n_eff
 - separate n_eff in general and MCSE using n_eff
 - more dicussion on infinite variances and non-finite means
 - rank normalization
 - unequal variances, min \hat[E]?

Implementation todo
- Bayesplot
  + rank plots
  + local relative efficiency plots
    - small intervals
    - quantiles
- Default output?
  + for all quantities
  + Median, Mad
- C++
  + Rank-normalization 
  + Folded
