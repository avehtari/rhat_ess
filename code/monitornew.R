# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2018 Aki Vehtari, Paul Bürkner
#
# RStan is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# RStan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# stan_prob_autocovariance <- function(v) { 
#   .Call("stan_prob_autocovariance", v)
# }

fft_next_good_size <- function(N) {
  # Find the optimal next size for the FFT so that
  # a minimum number of zeros are padded.
  if (N <= 2)
    return(2)
  while (TRUE) {
    m = N
    while ((m %% 2) == 0) m = m / 2
    while ((m %% 3) == 0) m = m / 3
    while ((m %% 5) == 0) m = m / 5
    if (m <= 1)
      return(N)
    N = N + 1
  }
}

autocovariance <- function(y) {
  # Compute autocovariance estimates for every lag for the specified
  # input sequence using a fast Fourier transform approach.
  N <- length(y)
  M <- fft_next_good_size(N)
  Mt2 <- 2 * M
  yc <- y - mean(y)
  yc <- c(yc, rep.int(0, Mt2 - N))
  transform <- fft(yc)
  ac <- fft(Conj(transform) * transform, inverse = TRUE)
  ac <- Re(ac)[1:N] / (N * 2 * seq(N, 1, by = -1))
  ac
}
   
autocorrelation <- function(y) {
  # Compute autocorrelation estimates for every lag for the specified
  # input sequence using a fast Fourier transform approach.
  ac <- autocovariance(y)
  ac <- ac / ac[1]
}

z_scale <- function(x) {
  S <- length(x)
  r <- rank(x, ties.method = 'average')
  z <- qnorm((r - 1 / 2) / S)
  if (!is.null(dim(x))) {
    # output should have the input dimension
    z <- array(z, dim = dim(x), dimnames = dimnames(x))
  }
  z
}

u_scale <- function(x) {
  S <- length(x)
  r <- rank(x, ties.method = 'average')
  u <- (r - 1 / 2) / S
  if (!is.null(dim(x))) {
    # output should have the input dimension
    u <- array(u, dim = dim(x), dimnames = dimnames(x))
  }
  u
}

r_scale <- function(x) {
  S <- length(x)
  r <- rank(x, ties.method = 'average')
  if (!is.null(dim(x))) {
    # output should have the input dimension
    r <- array(r, dim = dim(x), dimnames = dimnames(x))
  }
  r
}

ess_rfun <- function(sims) {
  # Compute the effective sample size for samples of several chains 
  # for one parameter; see the C++ code of function  
  # effective_sample_size in chains.cpp 
  # Args:
  #   sims: a 2D array _without_ warmup samples (# iter * # chains) 
  # Returns:
  #   A single numeric value
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  chains <- ncol(sims)
  n_samples <- nrow(sims)

  acov <- lapply(seq_len(chains), FUN = function(i) autocovariance(sims[, i])) 
  acov <- do.call(cbind, acov)
  chain_mean <- apply(sims, 2, mean)
  mean_var <- mean(acov[1, ]) * n_samples / (n_samples - 1) 
  var_plus <- mean_var * (n_samples - 1) / n_samples
  if (chains > 1) 
    var_plus <- var_plus + var(chain_mean)
  
  # Geyer's initial positive sequence
  rho_hat_t <- rep.int(0, n_samples)
  t <- 0
  rho_hat_even <- 1;
  rho_hat_t[t + 1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
  rho_hat_t[t + 2] <- rho_hat_odd
  t <- 2  
  while (t < nrow(acov) - 1 && !is.nan(rho_hat_even + rho_hat_odd) && 
         (rho_hat_even + rho_hat_odd > 0)) {
    rho_hat_even = 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
    rho_hat_odd = 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_t[t + 1] <- rho_hat_even
      rho_hat_t[t + 2] <- rho_hat_odd
    }
    t <- t + 2
  }
  max_t <- t
  # Geyer's initial monotone sequence
  t <- 2
  while (t <= max_t - 2) {  
    if (rho_hat_t[t + 1] + rho_hat_t[t + 2] >
  	rho_hat_t[t - 1] + rho_hat_t[t]) {
      rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2;
      rho_hat_t[t + 2] = rho_hat_t[t + 1];
    }
    t <- t + 2
  }
  ess <- chains * n_samples
  ess <- ess / (-1 + 2 * sum(rho_hat_t[1:max_t]))
  ess 
} 

rhat_rfun <- function(sims) {
  # Compute the rhat convergence diagnostic for a single parameter
  # For split-rhat, just call this with splitted chains
  # Args:
  #   sims: a 2D array _without_ warmup samples (# iter * # chains) 
  # Returns:
  #   A single numeric value
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  chain_mean <- numeric(chains)
  chain_var <- numeric(chains)
  for (i in seq_len(chains)) {
    chain_mean[i] <- mean(sims[, i])
    chain_var[i] <- var(sims[, i])
  } 
  var_between <- n_samples * var(chain_mean)
  var_within <- mean(chain_var) 
  sqrt((var_between/var_within + n_samples - 1) / n_samples)
} 

quantile_mcse <- function(sims, prob = NULL) {
  # compute Markov-chain SE of quantiles for a single parameter
  # prob must be a single quantile
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  I <- sims < quantile(sims, prob)
  Seff <- ess_rfun(z_scale(split_chains(I)))
  q <- c(0.1586553, 0.8413447, 0.05, 0.95)
  a <- qbeta(q, Seff * prob + 1, Seff * (1 - prob) + 1)
  ssims <- sort(sims)
  S <- length(ssims)
  th1 <- ssims[max(round(a[1] * S), 1)]
  th2 <- ssims[min(round(a[2] * S), S)]
  mcse <- (th2 - th1) / 2
  th1 <- ssims[max(round(a[3] * S), 1)]
  th2 <- ssims[min(round(a[4] * S), S)]
  data.frame(mcse = mcse, q05 = th1, q95 = th2, Seff = Seff)
}

split_chains <- function(sims) {
  # split Markov chains
  # Args:
  #   sims: a 2D array of samples (# iter * # chains) 
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  niter <- dim(sims)[1]
  half <- floor(niter / 2)
  cbind(sims[1:half, ], sims[(half + 1):niter, ])
}

monitornew <- function(sims, warmup = 0, probs = c(0.05, 0.50, 0.95)) { 
  # print the summary for a general simulation results 
  # of 3D array: # iter * # chains * # parameters 
  # Args:
  #   sims: a 3D array described above 
  #   warmup: the number of iterations used for warmup 
  #   probs: probs of summarizing quantiles 
  #   print: print out the results
  # Return: 
  #   A summary matrix of class 'simsummary'
  if (inherits(sims, "stanfit")) {
    chains <- sims@sim$chains
    iter <- sims@sim$iter
    warmup <- sims@sim$warmup
    parnames <- names(sims)
    sims <- as.array(sims)
  } else {
    dim_sims <- dim(sims)
    if (is.null(dim_sims)) {
      dim(sims) <- c(length(sims), 1, 1) 
    } else if (length(dim_sims) == 2) {
      dim(sims) <- c(dim_sims, 1)
    } else if (length(dim_sims) > 3) {
      stop("'sims' has more than 3 dimensions") 
    }
    parnames <- dimnames(sims)[[3]]
    if (is.null(parnames)) {
      parnames <- paste0("V", seq_len(dim(sims)[3]))
    }
    iter <- dim(sims)[1]
    chains <- dim(sims)[2]
    if (warmup > dim(sims)[1]) {
      stop("warmup is larger than the total number of iterations")
    }
    if (warmup >= 1) {
      sims <- sims[-seq_len(warmup), , , drop = FALSE] 
    }
  }
  
  mcse_fun <- function(p, sims) quantile_mcse(sims, p)$mcse
  summary <- vector("list", length(parnames))
  summary <- setNames(summary, parnames)
  for (i in seq_along(summary)) {
    sims_i <- sims[, , i]
    quan <- unname(quantile(sims_i, probs = probs))
    mcse <- sapply(probs, mcse_fun, sims_i)
    
    zsims_split <- z_scale(split_chains(sims_i))
    zsplit_rhat <- rhat_rfun(zsims_split)
    zsplit_ess <- ess_rfun(zsims_split)
    zsplit_ress <- zsplit_ess / prod(dim(sims_i))
    
    sims_folded <- abs(sims_i - median(sims_i))
    zsims_folded_split <- z_scale(split_chains(sims_folded))
    zfsplit_rhat <- rhat_rfun(zsims_folded_split)
    zfsplit_ess <- ess_rfun(zsims_folded_split)
    zfsplit_ress <- zfsplit_ess / prod(dim(sims_i))
    rhat <- max(zsplit_rhat, zfsplit_rhat)

    summary[[i]] <- c(quan, mcse, rhat, zsplit_ress, zfsplit_ress)
  }
  
  summary <- do.call(rbind, summary)
  probs_str <- paste0("Q", probs * 100)
  mcse_str <- paste0("SE_", probs_str)
  colnames(summary) <- c(probs_str, mcse_str, "Rhat", "Bulk_Reff", "Tail_Reff")
  rownames(summary) <- parnames
  structure(
    summary,
    chains = chains,
    iter = iter,
    warmup = warmup,
    class = "simsummary" 
  )
} 

print.simsummary <- function(x, digits = 2, ...) {
  atts <- attributes(x)
  rm_atts <- c("chains", "iter", "warmup")
  attributes(x)[rm_atts] <- NULL
  cat(
    "Inference for the input samples (", atts$chains, 
    " chains: each with iter = ", atts$iter, 
    "; warmup = ", atts$warmup, "):\n\n", sep = ""
  )
  class(x) <- "matrix"
  print(round(x, digits), ...)
  cat(
    "\nFor each parameter, Bulk_Reff and Tail_Reff are crude measures of relative\n",
    "effective sample size for bulk and tail quantities respectively (good mixing\n",
    "Reff > 0.1), and Rhat is the potential scale reduction factor on rank normalized\n",
    "split chains (at convergence, Rhat = 1).\n", sep = ""
  )
  invisible(x)
}

# monitor_simple <- function(sims, ...) {
#   out <- monitornew(sims, warmup = 0, probs = 0.5, print = FALSE, ...)
#   out <- as.data.frame(out)
#   out$par <- seq_len(nrow(out))
#   out
# }

# outdated version of 'monitornew'
# monitornew <- function(sims, warmup = floor(dim(sims)[1] / 2), 
#                        probs = c(0.025, 0.25, 0.50, 0.75, 0.975), 
#                        digits_summary = 1, print = TRUE, ...) { 
#   # print the summary for a general simulation results 
#   # of 3-d array: # iter * # chains * # parameters 
#   # Args:
#   #   sims: a 3-d array described above 
#   #   warmup: the number of iterations used for warmup 
#   #   probs: probs of summarizing quantiles 
#   #   print: print out the results
#   # 
#   # Return: 
#   #   A summary array  
#   dim_sims <- dim(sims)
#   if (is.null(dim_sims))
#       dim(sims) <- c(length(sims), 1, 1)
#   if (length(dim_sims) == 2)
#       dim(sims) <- c(dim_sims, 1)
#   if (length(dim_sims) > 3) 
#     stop("'sims' has more than 3 dimensions")
#   dim_sims <- dim(sims)
#   if (warmup > dim_sims[1])
#     stop("warmup is larger than the total number of iterations")
#   if (is(sims, "stanfit")) {
#     warmup <- 0L
#     sims <- as.array(sims)
#   }
#   dimnames_sims <- dimnames(sims)
#   parnames <- dimnames_sims[[3]]
#   num_par <- dim_sims[3]
#   
#   if (is.null(parnames)) parnames <- paste0("V", 1:num_par)
#   sims_wow <- if (warmup >= 1) apply(sims, c(2, 3), FUN = function(x) x[-(1:warmup)]) else sims 
#   dim_sims <- dim(sims_wow)
#   n_samples <- dim_sims[1]
#   n_chains <- dim_sims[2]
#   half_n <- floor(n_samples / 2)
#   idx_2nd <- n_samples - half_n + 1
#   m <- apply(sims_wow, 3, mean)
#   sd <- sapply(1:num_par, FUN = function(i) sd(as.vector(sims_wow[,,i]))) 
#   quan <- lapply(1:num_par, FUN = function(i) quantile(sims_wow[,,i], probs = probs))
#   probs_str <- names(quan[[1]])
#   quan <- do.call(rbind, quan)
#   rhat <- sapply(1:num_par, FUN = function(i) rhat_rfun(sims_wow[,,i]))
#   split_rhat <- sapply(1:num_par, FUN = function(i)
#     rhat_rfun(cbind(sims_wow[1:half_n,,i],sims_wow[idx_2nd:n_samples,,i])))
#   ess <- sapply(1:num_par, FUN = function(i) ess_rfun(sims_wow[,,i]))
#   ress <-ess/n_samples/n_chains
#   split_ess <- sapply(1:num_par, FUN = function(i)
#     ess_rfun(cbind(sims_wow[1:half_n,,i],sims_wow[idx_2nd:n_samples,,i])))
#   sem <- sd / sqrt(ess)
#   zrhat <- sapply(1:num_par, FUN = function(i) rhat_rfun(z_scale(sims_wow[,,i])))
#   zsplit_rhat <- sapply(1:num_par, FUN = function(i)
#       rhat_rfun(z_scale(cbind(sims_wow[1:half_n,,i],sims_wow[idx_2nd:n_samples,,i]))))
#   sims_centered <- sweep(sims_wow,3,apply(sims_wow,3,median))
#   sims_folded <- abs(sims_centered)
#   sims_med <- (sims_centered<0)*1
#   sims_mad <- (sweep(sims_folded,3,apply(sims_folded,3,median))<0)*1
#   zfsplit_rhat <- sapply(1:num_par, FUN = function(i)
#       rhat_rfun(z_scale(cbind(sims_folded[1:half_n,,i],sims_folded[idx_2nd:n_samples,,i]))))
#   ## zfsplit_rhat2 <- sapply(1:num_par, FUN = function(i) {
#   ##     sims_split <- cbind(sims_wow[1:half_n,,i],sims_wow[idx_2nd:n_samples,,i])
#   ##     sims_folded <- abs(sweep(sims_split,2,apply(sims_split,2,median)))
#   ##     rhat_rfun(z_scale(sims_folded))
#   ## })
#   zess <- sapply(1:num_par, FUN = function(i) ess_rfun(z_scale(sims_wow[,,i])))
#   zsplit_ess <- sapply(1:num_par, FUN = function(i)
#     ess_rfun(z_scale(cbind(sims_wow[1:half_n,,i],sims_wow[idx_2nd:n_samples,,i]))))
#   zsplit_ress <- zsplit_ess/n_samples/n_chains
#   zfsplit_ess <- sapply(1:num_par, FUN = function(i)
#       ess_rfun(z_scale(cbind(sims_folded[1:half_n,,i],sims_folded[idx_2nd:n_samples,,i]))))
#   zfsplit_ress <- zfsplit_ess/n_samples/n_chains
#   medsplit_ess <- sapply(1:num_par, FUN = function(i)
#       ess_rfun(z_scale(cbind(sims_med[1:half_n,,i],sims_med[idx_2nd:n_samples,,i]))))
#   medsplit_ress <- medsplit_ess/n_samples/n_chains
#   madsplit_ess <- sapply(1:num_par, FUN = function(i)
#       ess_rfun(z_scale(cbind(sims_mad[1:half_n,,i],sims_mad[idx_2nd:n_samples,,i]))))
#   madsplit_ress <- madsplit_ess/n_samples/n_chains
#   sem <- sd / sqrt(ess)
#   
#   summary <- cbind(m, sem, sd, quan, ess, ress, split_ess, zess, zsplit_ess, zsplit_ress, rhat, split_rhat, zrhat, zsplit_rhat, zfsplit_rhat, zfsplit_ess, zfsplit_ress, medsplit_ess, medsplit_ress, madsplit_ess, madsplit_ress)
#   colnames(summary) <- c("mean", "se_mean", "sd", probs_str, "neff", "reff", "sneff", "zneff", "zsneff", "zsreff", "Rhat", "sRhat", "zRhat", "zsRhat", "zfsRhat", "zfsneff", "zfsreff", "medsneff", "medsreff", "madsneff", "madsreff")
#   rownames(summary) <- parnames 
#   if (print) {
#     cat("Inference for the input samples (")
#     cat(dim_sims[2], " chains: each with iter=", dim_sims[1], "; warmup=", warmup, "):\n\n", sep = "")
#     # round n_eff to integers
#     summary[, 'neff'] <- round(summary[, 'neff'], 0)
#     summary[, 'sneff'] <- round(summary[, 'sneff'], 0)
#     summary[, 'zneff'] <- round(summary[, 'zneff'], 0)
#     summary[, 'zsneff'] <- round(summary[, 'zsneff'], 0)
#     print(round(summary, digits_summary), ...)
#  
#     cat("\nFor each parameter, n_eff is a crude measure of effective sample size,\n", 
#         "and Rhat is the potential scale reduction factor on split chains (at \n",
#         "convergence, Rhat=1).\n", sep = '')
#   } 
#   invisible(summary) 
# } 
