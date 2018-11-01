# Copyright (C) 2018 Aki Vehtari, Paul Bürkner

plotranknorm <- function(theta, n, m = 1, interval = NULL) {
    df <- data.frame(theta = theta) %>%
        mutate(gid = gl(m, n),
               r = r_scale(theta),
               u = u_scale(theta),
               z = z_scale(theta))
    blue <- color_scheme_get(scheme ='blue', i = 4)[[1]]
    p2 <- ggplot(df, aes(x=theta, grp=gid)) +
        stat_ecdf(color=blue, alpha=0.5, pad = FALSE) +
        labs(x=TeX('$\\theta$'), y='ECDF')
    p3 <- ggplot(df, aes(x=r/(n*m), grp=gid)) +
        stat_ecdf(color=blue, alpha=0.5, pad = FALSE) +
        labs(x='Scaled rank', y='ECDF')
    if (is.null(interval)) {
      p1 <- mcmc_hist(as.data.frame(theta)) +
         xlab(TeX('$\\theta$'))
      p4 <- ggplot(df, aes(x=z, grp=gid)) +
        stat_ecdf(color=blue, alpha=0.5, pad = FALSE) +
        labs(x='z', y='ECDF')
      grid.arrange(p1, p2, p3, p4, nrow = 1)
    } else {
    	df <- mutate(df,
          psd = sqrt((r+1)*(n-r+1)/(n+2)^2/(n+3)), #sqrt(n*u*(1-u))/n,
          pm2sd = (r+1)/(n+2)-1.96*psd,
          pp2sd = (r+1)/(n+2)+1.96*psd,
          p975 = qbeta(0.975,r+1,n-r+1),
          p025 = qbeta(0.025,r+1,n-r+1))
        p2b <- p2 + geom_linerange(data=df, aes(ymin=p025, ymax=p975), color=blue) +
          ylab('ECDF + beta 95% interval')
        p3b <- p3 + geom_linerange(data=df, aes(ymin=p025, ymax=p975), color=blue) + 
          ylab('ECDF + beta 95% interval')
        p3c <- p3 + geom_linerange(data=df, aes(ymin=pm2sd, ymax=pp2sd), color=blue) +
          ylab('ECDF + normal approximated 95% interval')
        grid.arrange(p2b, p3b, p3c, nrow = 1)
    }
}

mcmc_hist_r_scale <- function(x, nbreaks = 50) {
  max <- prod(dim(x)[1:2])
  mcmc_hist(
    r_scale(x), 
    breaks = seq(0, max, by = max / nbreaks) + 0.5
  )
}

plot_rhat <- function(res) {
  p1 <- ggplot(res, aes(x=par, y=sRhat)) + 
    geom_point() + 
    ggtitle('Classic split-Rhat') + 
    geom_hline(yintercept = 1.005, linetype = 'dashed') + 
    geom_hline(yintercept = 1) + 
    ylim(c(.99 ,1.26))

  p2 <- ggplot(res, aes(x=par, y=zsRhat)) + 
    geom_point() + 
    ggtitle('Rank normalized split-Rhat') + 
    geom_hline(yintercept = 1.005, linetype = 'dashed') + 
    geom_hline(yintercept = 1) + 
    ylim(c(.99,1.26))
  
  p3 <- ggplot(res, aes(x=par, y=zfsRhat)) + 
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
  
  p1 <- ggplot(res, aes(x=par, y=reff)) + 
    geom_point() + 
    ggtitle('Classic R_eff (reff)') + 
    geom_hline(yintercept = c(0,1)) +
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits) 
  
  p2 <- ggplot(res, aes(x=par, y=zsreff)) + 
    geom_point() + 
    ggtitle('New Bulk-R_eff') + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits)
  
  p3 <- ggplot(res, aes(x=par, y=zfsreff)) + 
    geom_point() + 
    ggtitle('New Tail-R_eff') + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits)
  
  p4 <- ggplot(res, aes(x=par, y=medsreff)) + 
    geom_point() + 
    ggtitle('Median-R_eff') + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits)
  
  p5 <- ggplot(res, aes(x=par, y=madsreff)) + 
    geom_point() + 
    ggtitle('MAD-R_eff') + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_y_continuous(breaks = ybreaks, limits = ylimits)
  
  blank <- grid::grid.rect(gp = grid::gpar(col="white"), draw = FALSE)
  grid.arrange(p1, p2, p3, blank, p4, p5, nrow = 2)
}

plot_local_reff <- function(samp = NULL, par = NULL, nalpha = NULL) {
  delta <- 1 / nalpha
  alphas <- seq(0, 1 - delta, by = delta)
  zsreffs <- rep(NA, length(alphas))
  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    tmp <- samp[, , par]
    I <- tmp >= quantile(tmp, alpha) & tmp < quantile(tmp, alpha + delta)
    rs <- monitor_simple(I)
    zsreffs[i] <- rs$zsreff[1]
  }
  df <- data.frame(
    quantile = seq(0, 1, by = delta), 
    zsreff = c(zsreffs, zsreffs[nalpha])
  )
  ggplot(data=df, aes(x=quantile, y=zsreff)) + 
    geom_step() + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_x_continuous(breaks = seq(0,1, by=0.1)) + 
    scale_y_continuous(
      breaks = seq(0, 1, by=0.25), 
      limits = c(0, 1.1)
    ) +
    labs(x='Quantile', y='R_eff')
}

plot_quantile_reff <- function(samp = NULL, par = NULL, nalpha=NULL) {
  delta <- 1 / nalpha
  alphas <- seq(delta, 1 - delta, by = delta)
  zsreffs <- rep(NA, length(alphas))
  for (i in seq_along(alphas)) {
    alpha <- alphas[i]
    tmp <- samp[, , par]
    I <- tmp < quantile(tmp, alpha)
    rs <- monitor_simple(I)
    zsreffs[i] <- rs$zsreff[1]
  }
  df <- data.frame(
    quantile = seq(delta, 1 - delta, by = delta), 
    zsreff = zsreffs
  )
  ymax <- max(1, round(max(zsreffs, na.rm = TRUE) + 0.15, 1))
  ggplot(df, aes(x=quantile, y=zsreff)) + 
    geom_point() + 
    geom_hline(yintercept = c(0,1)) + 
    geom_hline(yintercept = 0.1, linetype = 'dashed') + 
    scale_x_continuous(breaks = seq(0,1, by=0.1)) + 
    scale_y_continuous(
      breaks = seq(0, ymax, by = 0.25), 
      limits = c(0, ymax)
    ) +
    labs(x='Quantile', y='R_eff')
}
