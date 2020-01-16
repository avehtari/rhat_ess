# Copyright (C) 2018 Aki Vehtari, Paul Bürkner

plot_ranknorm <- function(theta, n, m = 1, interval = FALSE) {
  df <- data.frame(theta = theta) %>%
    mutate(
      gid = gl(m, n),
      r = r_scale(theta),
      u = u_scale(theta),
      z = z_scale(theta)
    )
  size <- 1.1
  alpha <- if (m > 1) 0.5 else 1
  blue <- color_scheme_get(scheme = 'blue', i = 4)[[1]]
  p2 <- ggplot(df, aes(x = theta, grp = gid)) +
    stat_ecdf(color = blue, size = size, alpha = alpha, pad = FALSE) +
    labs(x = 'theta', y = 'ECDF')
  p3 <- ggplot(df, aes(x = r / (n * m), grp = gid)) +
    stat_ecdf(color = blue, size = size, alpha = alpha, pad = FALSE) +
    labs(x = 'Scaled rank', y = 'ECDF')
  
  if (interval) {
    df <- df %>% mutate(
      psd = sqrt((r + 1) * (n - r + 1) / (n + 2)^2 / (n + 3)),
      pm2sd = (r + 1) / (n + 2) - 1.96 * psd,
      pp2sd = (r + 1) / (n + 2) + 1.96 * psd,
      p975 = qbeta(0.975, r + 1, n - r + 1),
      p025 = qbeta(0.025, r + 1, n - r + 1)
    )
    p2b <- p2 + 
      geom_line(data = df, aes(y = p025), color = blue) +
      geom_line(data = df, aes(y = p975), color = blue) +
      ylab('ECDF + beta 95% interval')
    p3b <- p3 + 
      geom_line(data = df, aes(y = p025), color = blue) +
      geom_line(data = df, aes(y = p975), color = blue) + 
      ylab('ECDF + beta 95% interval')
    p3c <- p3 + 
      geom_line(data = df, aes(y = pm2sd), color = blue) +
      geom_line(data = df, aes(y = pp2sd), color = blue) +
      ylab('ECDF + normal approximated 95% interval')
    out <- gridExtra::grid.arrange(p2b, p3b, p3c, nrow = 1)
  } else {
    p1 <- bayesplot::mcmc_hist(as.data.frame(theta)) +
      xlab('theta')
    p4 <- ggplot(df, aes(x = z, grp = gid)) +
      stat_ecdf(color = blue, size = size, alpha = alpha, pad = FALSE) +
      labs(x = 'z', y = 'ECDF')
    out <- gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 1)
  }
  invisible(out)
}

mcmc_hist_r_scale <- function(x, nbreaks = 50, ...) {
  max <- prod(dim(x)[1:2])
  bayesplot::mcmc_hist(
    r_scale(x), 
    breaks = seq(0, max, by = max / nbreaks) + 0.5,
    ...
  ) +
    theme(axis.line.y = element_blank())
}

plot_rhat <- function(res) {
  res$par <- rownames(res)
  p1 <- ggplot(res, aes(x = par, y = sRhat)) + 
    geom_point() + 
    ggtitle('Classic split-Rhat') + 
    geom_hline(yintercept = 1.005, linetype = 'dashed') + 
    geom_hline(yintercept = 1) + 
    ylim(c(.99 ,1.26)) +
    xlab("Parameters") +
    theme(axis.text = element_blank())

  p2 <- ggplot(res, aes(x = par, y = zsRhat)) + 
    geom_point() + 
    ggtitle('Rank normalized split-Rhat') + 
    geom_hline(yintercept = 1.005, linetype = 'dashed') + 
    geom_hline(yintercept = 1) + 
    ylim(c(.99,1.26)) +
    xlab("Parameters") +
    theme(axis.text = element_blank())
  
  p3 <- ggplot(res, aes(x = par, y = zfsRhat)) + 
    geom_point() + 
    ggtitle('Folded rank norm. split-Rhat') + 
    geom_hline(yintercept = 1.005, linetype = 'dashed') + 
    geom_hline(yintercept = 1) + 
    ylim(c(.99,1.26)) +
    xlab("Parameters") +
    theme(axis.text = element_blank())
  
  gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
}

plot_ess <- function(res) {
  att <- attributes(res)
  max_seff <- with(res,
    max(c(seff, zsseff, zfsseff, medsseff, madsseff), na.rm = TRUE)
  )
  S <- att$iter * att$chains
  if (!length(S)) S <- max_seff
  ymax <- round(max(S, max_seff * 1.15))
  ylimits <- c(0, ymax)
  res$par <- rownames(res)
  
  p1 <- ggplot(res, aes(x = par, y = seff)) + 
    geom_point() + 
    ggtitle('Classic ESS') + 
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 400, linetype = 'dashed') + 
    scale_y_continuous(limits = ylimits) +
    xlab("Parameters") +
    theme(axis.text = element_blank())
  
  p2 <- ggplot(res, aes(x = par, y = zsseff)) + 
    geom_point() + 
    ggtitle('Bulk-ESS') + 
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = 400, linetype = 'dashed') + 
    scale_y_continuous(limits = ylimits) +
    xlab("Parameters") +
    theme(axis.text = element_blank())
  
  p3 <- ggplot(res, aes(x = par, y = tailseff)) + 
    geom_point() + 
    ggtitle('Tail-ESS') + 
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = 400, linetype = 'dashed') + 
    scale_y_continuous(limits = ylimits) +
    xlab("Parameters") +
    theme(axis.text = element_blank())
  
  p4 <- ggplot(res, aes(x = par, y = medsseff)) + 
    geom_point() + 
    ggtitle('Median-ESS') + 
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = 400, linetype = 'dashed') + 
    scale_y_continuous(limits = ylimits) +
    xlab("Parameters") +
    theme(axis.text = element_blank())
  
  p5 <- ggplot(res, aes(x = par, y = madsseff)) + 
    geom_point() + 
    ggtitle('MAD-ESS') + 
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = 400, linetype = 'dashed') + 
    scale_y_continuous(limits = ylimits) +
    xlab("Parameters") +
    theme(axis.text = element_blank())
  
  blank <- grid::grid.rect(gp = grid::gpar(col = "white"), draw = FALSE)
  gridExtra::grid.arrange(p1, p2, p3, blank, p4, p5, nrow = 2)
}

plot_local_ess <- function(fit, par, nalpha = 20, rank = TRUE) {
  if (length(par) != 1L) {
    stop("'par' should be of length 1.")
  }
  if (inherits(fit, "stanfit")) {
    if (!is.character(par)) {
      par <- names(fit)[par]
    }
    sims <- as.array(fit, pars = par)[, , 1]
    params <- as.data.frame(fit, pars = par)
    sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
    divergent <- do.call(rbind, sampler_params)[, 'divergent__']
    max_depth <- attr(fit@sim$samples[[1]], "args")$control$max_treedepth
    treedepths <- do.call(rbind, sampler_params)[, 'treedepth__']
    params$divergent <- divergent
    params$max_depth <- (treedepths == max_depth) * 1
    params$urank <- u_scale(params[, par])
    params$value <- params[, par]
  } else {
    if (!is.character(par)) {
      par <- dimnames(fit)[[3]][par]
    }
    sims <- fit[, , par]
    params <- data.frame(value = as.vector(sims))
    params$divergent <- 0
    params$max_depth <- 0
    params$urank <- u_scale(params$value)
  }   
  
  # compute local Seff
  delta <- 1 / nalpha
  alphas <- seq(0, 1 - delta, by = delta)
  zsseffs <- rep(NA, length(alphas))
  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    I <- sims > quantile(sims, alpha) & sims <= quantile(sims, alpha + delta)
    zsseffs[i] <- ess_rfun(split_chains(I))
  }
  S <- prod(dim(I))
  
  # create the plot
  df <- data.frame(
    quantile = seq(0, 1, by = delta),
    value = quantile(params$value, seq(0, 1, by = delta)),
    zsseff = c(zsseffs, zsseffs[nalpha])
  )
  ymax <- max(S, round(max(zsseffs, na.rm = TRUE) * 1.15, 1))
  xname <- if (rank) "quantile" else "value"
  xrug <- if (rank) "urank" else "value"
  out <- ggplot(data = df, aes_string(x = xname, y = "zsseff")) +
    geom_step() + 
    geom_hline(yintercept = c(0, 1)) + 
    geom_hline(yintercept = 400, linetype = 'dashed') + 
    scale_y_continuous(
      breaks = seq(0, ymax, by = round(0.25*S)), 
      limits = c(0, ymax)
    ) +
    geom_rug(
      data = params[params$divergent == 1, ], 
      aes_string(x = xrug, y = NULL), sides = "b", color = "red"
    ) +
    geom_rug(
      data = params[params$max_depth == 1, ], 
      aes_string(x = xrug, y = NULL), sides = "b", color = "orange"
    ) +
    ylab('ESS for small intervals')
  if (rank) {
    out <- out +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
      xlab('Quantile')
  } else {
    out <- out + xlab(par)
  }
  out
}

plot_quantile_ess <- function(fit, par, nalpha = 20, rank = TRUE) {
  if (length(par) != 1L) {
    stop("'par' should be of length 1.")
  }
  if (inherits(fit, "stanfit")) {
    if (!is.character(par)) {
      par <- names(fit)[par]
    }
    sims <- as.array(fit, pars = par)[, , 1]
    params <- as.data.frame(fit, pars = par)
    sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
    divergent <- do.call(rbind, sampler_params)[, 'divergent__']
    max_depth <- attr(fit@sim$samples[[1]], "args")$control$max_treedepth
    treedepths <- do.call(rbind, sampler_params)[, 'treedepth__']
    params$divergent <- divergent
    params$max_depth <- (treedepths == max_depth) * 1
    params$urank <- u_scale(params[, par])
    params$value <- params[, par]
  } else {
    if (!is.character(par)) {
      par <- dimnames(fit)[[3]][par]
    }
    sims <- fit[, , par]
    params <- data.frame(value = as.vector(sims))
    params$divergent <- 0
    params$max_depth <- 0
    params$urank <- u_scale(params$value)
  }
  
  # compute quantile Seff
  delta <- 1 / nalpha
  alphas <- seq(delta, 1 - delta, by = delta)
  zsseffs <- rep(NA, length(alphas))
  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    I <- sims <= quantile(sims, alpha)
    zsseffs[i] <- ess_rfun(split_chains(I))
  }
  S <- prod(dim(I))
  
  # create the plot
  df <- data.frame(
    quantile = seq(delta, 1 - delta, by = delta), 
    value = quantile(params$value, seq(delta, 1 - delta, by = delta)),
    zsseff = zsseffs
  )
  ymax <- max(S, round(max(zsseffs, na.rm = TRUE) * 1.15, 1))
  xname <- if (rank) "quantile" else "value"
  xrug <- if (rank) "urank" else "value"
  out <- ggplot(data = df, aes_string(x = xname, y = "zsseff")) +
    geom_point() + 
    geom_hline(yintercept = c(0, 1)) + 
    geom_hline(yintercept = 400, linetype = 'dashed') + 
    scale_y_continuous(
      breaks = seq(0, ymax, by = round(0.25*S)), 
      limits = c(0, ymax)
    ) +
    geom_rug(
      data = params[params$divergent == 1, ], 
      aes_string(x = xrug, y = NULL), sides = "b", color = "red"
    ) +
    geom_rug(
      data = params[params$max_depth == 1, ], 
      aes_string(x = xrug, y = NULL), sides = "b", color = "orange"
    ) +
    ylab("ESS for quantiles")
  if (rank) {
    out <- out +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
      xlab('Quantile')
  } else {
    out <- out + xlab(par)
  }
  out
}

plot_change_ess <- function(fit, par, breaks = seq(0.1, 1, 0.05), 
                            yaxis = c("absolute", "relative")) {
  if (length(par) != 1L) {
    stop("'par' should be of length 1.")
  }
  yaxis <- match.arg(yaxis)
  if (inherits(fit, "stanfit")) {
    if (!is.character(par)) {
      par <- names(fit)[par]
    }
    sims <- as.array(fit, pars = par)[, , 1]
  } else {
    if (!is.character(par)) {
      par <- dimnames(fit)[[3]][par]
    }
    sims <- fit[, , par]
  }
  
  iter_breaks <- round(breaks * NROW(sims))
  nbreaks <- length(iter_breaks)
  bulk_seff <- tail_seff <- bulk_reff <- 
    tail_reff <- rep(NA, length(nbreaks))
  for (i in seq_along(iter_breaks)) {
    sims_i <- sims[seq_len(iter_breaks[i]), ]
    nsamples <- prod(dim(sims_i))
    bulk_seff[i] <- ess_rfun(z_scale(split_chains(sims_i)))
    tail_seff[i] <- ess_tail(sims_i)
    bulk_reff[i] <- bulk_seff[i] / nsamples
    tail_reff[i] <- tail_seff[i] / nsamples
  }
  df <- data.frame(
    breaks = breaks,
    ndraws = iter_breaks * NCOL(sims),
    seff = c(bulk_seff, tail_seff),
    reff = c(bulk_reff, tail_reff), 
    type = rep(c("bulk", "tail"), each = nbreaks)
  )
  blues <- bayesplot::color_scheme_get(scheme = "blue", i = c(4, 2))
  blues <- unname(unlist(blues))
  if (yaxis == "absolute") {
    out <- ggplot(df, aes(ndraws, seff, color = type)) +
      ylab("ESS") +
      geom_hline(yintercept = 0, linetype = 1) +
      geom_hline(yintercept = 400, linetype = 2)
  } else if (yaxis == "relative") {
    out <- ggplot(df, aes(ndraws, reff, color = type)) +
      ylab("Relative efficiency") +
      geom_hline(yintercept = 0, linetype = 2)
  }
  out +  
    geom_line() +
    geom_point() +
    xlab("Total number of draws") +
    scale_colour_manual(values = blues)
}
