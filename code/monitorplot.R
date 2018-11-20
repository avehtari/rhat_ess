# Copyright (C) 2018 Aki Vehtari, Paul Bürkner

plotranknorm <- function(theta, n, m = 1, interval = FALSE) {
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
    out <- grid.arrange(p2b, p3b, p3c, nrow = 1)
  } else {
    p1 <- mcmc_hist(as.data.frame(theta)) +
      xlab('theta')
    p4 <- ggplot(df, aes(x = z, grp = gid)) +
      stat_ecdf(color = blue, size = size, alpha = alpha, pad = FALSE) +
      labs(x = 'z', y = 'ECDF')
    out <- grid.arrange(p1, p2, p3, p4, nrow = 1)
  }
  invisible(out)
}

mcmc_hist_r_scale <- function(x, nbreaks = 50) {
  max <- prod(dim(x)[1:2])
  mcmc_hist(
    r_scale(x), 
    breaks = seq(0, max, by = max / nbreaks) + 0.5
  )
}

plot_rhat <- function(res) {
  res$par <- rownames(res)
  p1 <- ggplot(res, aes(x = par, y = sRhat)) + 
    geom_point() + 
    ggtitle('Classic split-Rhat') + 
    geom_hline(yintercept = 1.005, linetype = 'dashed') + 
    geom_hline(yintercept = 1) + 
    ylim(c(.99 ,1.26))

  p2 <- ggplot(res, aes(x = par, y = zsRhat)) + 
    geom_point() + 
    ggtitle('Rank normalized split-Rhat') + 
    geom_hline(yintercept = 1.005, linetype = 'dashed') + 
    geom_hline(yintercept = 1) + 
    ylim(c(.99,1.26))
  
  p3 <- ggplot(res, aes(x = par, y = zfsRhat)) + 
    geom_point() + 
    ggtitle('Folded rank norm. split-Rhat') + 
    geom_hline(yintercept = 1.005, linetype = 'dashed') + 
    geom_hline(yintercept = 1) + 
    ylim(c(.99,1.26))
  
  grid.arrange(p1, p2, p3, nrow = 1)
}

plot_reff <- function(res) {
  ymax <- 2
  ybreaks <- seq(0, ymax, by = 0.25)
  ylimits <- c(0, ymax)
  res$par <- rownames(res)
  
  p1 <- ggplot(res, aes(x = par, y = reff)) + 
    geom_point() + 
    ggtitle('Classic Reff') + 
    geom_hline(yintercept = c(0,1)) +
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits) 
  
  p2 <- ggplot(res, aes(x = par, y = zsreff)) + 
    geom_point() + 
    ggtitle('New Bulk-Reff') + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits)
  
  p3 <- ggplot(res, aes(x = par, y = zfsreff)) + 
    geom_point() + 
    ggtitle('New Tail-Reff') + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits)
  
  p4 <- ggplot(res, aes(x = par, y = medsreff)) + 
    geom_point() + 
    ggtitle('Median-Reff') + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits)
  
  p5 <- ggplot(res, aes(x = par, y = madsreff)) + 
    geom_point() + 
    ggtitle('MAD-Reff') + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits)
  
  blank <- grid::grid.rect(gp = grid::gpar(col = "white"), draw = FALSE)
  grid.arrange(p1, p2, p3, blank, p4, p5, nrow = 2)
}

plot_local_reff <- function(fit, par, nalpha = 20, rank = TRUE) {
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
  
  # compute local Reff
  delta <- 1 / nalpha
  alphas <- seq(0, 1 - delta, by = delta)
  zsreffs <- rep(NA, length(alphas))
  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    I <- sims > quantile(sims, alpha) & sims <= quantile(sims, alpha + delta)
    zsreffs[i] <- ess_rfun(z_scale(split_chains(I))) / prod(dim(I))
  }
  
  # create the plot
  df <- data.frame(
    quantile = seq(0, 1, by = delta),
    value = quantile(params$value, seq(0, 1, by = delta)),
    zsreff = c(zsreffs, zsreffs[nalpha])
  )
  ymax <- max(1, round(max(zsreffs, na.rm = TRUE) + 0.15, 1))
  xname <- if (rank) "quantile" else "value"
  xrug <- if (rank) "urank" else "value"
  out <- ggplot(data = df, aes_string(x = xname, y = "zsreff")) +
    geom_step() + 
    geom_hline(yintercept = c(0, 1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(
      breaks = seq(0, ymax, by = 0.25), 
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
    ylab('Reff of small intervals')
  if (rank) {
    out <- out +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
      xlab('Quantile')
  } else {
    out <- out + xlab(par)
  }
  out
}

plot_quantile_reff <- function(fit, par, nalpha = 20, rank = TRUE) {
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
  
  # compute quantile Reff
  delta <- 1 / nalpha
  alphas <- seq(delta, 1 - delta, by = delta)
  zsreffs <- rep(NA, length(alphas))
  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    I <- sims <= quantile(sims, alpha)
    zsreffs[i] <- ess_rfun(z_scale(split_chains(I))) / prod(dim(I))
  }
  
  # create the plot
  df <- data.frame(
    quantile = seq(delta, 1 - delta, by = delta), 
    value = quantile(params$value, seq(delta, 1 - delta, by = delta)),
    zsreff = zsreffs
  )
  ymax <- max(1, round(max(zsreffs, na.rm = TRUE) + 0.15, 1))
  xname <- if (rank) "quantile" else "value"
  xrug <- if (rank) "urank" else "value"
  out <- ggplot(data = df, aes_string(x = xname, y = "zsreff")) +
    geom_point() + 
    geom_hline(yintercept = c(0, 1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(
      breaks = seq(0, ymax, by = 0.25), 
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
    ylab("Reff of quantiles")
  if (rank) {
    out <- out +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + 
      xlab('Quantile')
  } else {
    out <- out + xlab(par)
  }
  out
}
