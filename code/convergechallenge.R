library(tidyverse)
theme_set(bayesplot::theme_default(base_family = "sans"))

blues <- bayesplot::color_scheme_get(scheme = "blue", i = c(4, 2))
blues <- unname(unlist(blues))

set.seed(12345)
dat <- data.frame(
  Iteration = rep(1:1000, 4),
  Chain = factor(rep(c(1:2, 1:2), each = 1000)),
  Cond = factor(rep(1:2, each = 2000))
) %>%
  mutate(
    Simulation = c(
      -2 + arima.sim(list(ar = 0.7), n = 1000, sd = 0.5),
      1 + arima.sim(list(ar = 0.7), n = 1000, sd = 0.5),
      -2 + 0.003 * 1:1000 + arima.sim(list(ar = 0.7), n = 1000, sd = 0.5),
      1 + -0.003 * 1:1000 + arima.sim(list(ar = 0.7), n = 1000, sd = 0.5)
    )
  )

ggplot(dat, aes(Iteration, Simulation, color = Chain)) +
  facet_wrap("Cond") +
  geom_line(cex = 0.5) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + 
  scale_color_manual(values = blues,
                     guide = FALSE)

ggsave(file = "code/convergechallenge.pdf", 
       height = 4, width = 8, units = "in")
