library(bagelR)
library(dplyr)
library(ggplot2)

OUT <- if (file.exists("simulation/speed.R")) "simulation" else "."

# ---- settings ----------------------------------------------------------------
T_END       <- 1000L
N_REP       <- 200L
M_GRID      <- c(50, 100, 200) 
M_FIXED     <- 100                  
H_GRID      <- c(0, 1, 4, 13)
PERIOD      <- 24
PROB_CHANGE <- 1e-5
P0          <- 1 - PROB_CHANGE
DECAY       <- 1
SIGMA       <- 1  
NU          <- 3 
IOTA        <- 3
DELTA_LEVEL  <- 1  
DELTA_SEASON <- 1

set.seed(1)
X <- rnorm(T_END)  
seasonal_row <- function(t, period, harmonics) {
  out <- numeric(2 * harmonics)
  for (k in seq_len(harmonics)) {
    out[2 * k - 1] <- sin(2 * pi * k * t / period)
    out[2 * k]     <- cos(2 * pi * k * t / period)
  }
  out
}
raw_joint_prior_season <- function(harmonics, mu0_level, delta0_level,
                                   delta_level, delta_season) {
  d1  <- 2 * harmonics + 1
  mu0 <- c(rep(0, 2 * harmonics), mu0_level, delta0_level - mu0_level)
  Sig <- matrix(0, d1 + 1, d1 + 1)
  if (harmonics > 0) diag(Sig)[seq_len(2 * harmonics)] <- delta_season
  i <- d1; j <- d1 + 1                                   
  Sig[i, i] <- delta_level
  Sig[i, j] <- Sig[j, i] <- -delta_level
  Sig[j, j] <- 2 * delta_level
  list(mu = mu0, Sig = Sig)
}
transform_matrix_season <- function(harmonics) {
  d1 <- 2 * harmonics + 1
  cbind(diag(d1), c(rep(0, 2 * harmonics), 1))          
}
make_feature_vector_season <- function(period, harmonics) {
  force(period); force(harmonics)
  function(t, tau) {
    sr <- seasonal_row(t, period, harmonics)
    if (tau == 0) return(matrix(c(sr, 1), ncol = 1))
    M <- matrix(c(sr, 1, 0), ncol = 1)
    if (t > tau) M[length(M), 1] <- 1
    M
  }
}
make_prior_season <- function(harmonics, mu0_level, delta0_level,
                              delta_level, delta_season, sigma2 = 1) {
  force(harmonics); force(mu0_level); force(delta0_level)
  force(delta_level); force(delta_season); force(sigma2)
  d1 <- 2 * harmonics + 1
  function(t) {
    if (t == 1)
      return(list(mu    = matrix(c(rep(0, 2 * harmonics), mu0_level), ncol = 1),
                  sigma = diag(c(rep(delta_season, 2 * harmonics), delta_level),
                               d1) / sigma2))
    rj <- raw_joint_prior_season(harmonics, mu0_level, delta0_level,
                                 delta_level, delta_season)
    list(mu = matrix(rj$mu, ncol = 1), sigma = rj$Sig / sigma2)
  }
}
make_transform_season <- function(harmonics) {
  force(harmonics)
  A <- transform_matrix_season(harmonics)
  function(tau) A
}

time_one_run <- function(variance, H, M) {
  fv <- make_feature_vector_season(PERIOD, H)
  pr <- make_prior_season(H, 0, 0, DELTA_LEVEL, DELTA_SEASON)
  tr <- make_transform_season(H)
  b <- if (variance == "known") {
    if (is.finite(M)) sequential_bagel_kv_approximate(fv, pr, tr, P0, DECAY, SIGMA, as.integer(M))
    else              sequential_bagel_kv_exact(fv, pr, tr, P0, DECAY, SIGMA, T_END + 10L)
  } else {
    if (is.finite(M)) sequential_bagel_uv_approximate(fv, pr, tr, P0, DECAY, NU, IOTA, as.integer(M))
    else              sequential_bagel_uv_exact(fv, pr, tr, P0, DECAY, NU, IOTA, T_END + 10L)
  }
  secs <- numeric(T_END)
  t0 <- Sys.time()   
  for (t in seq_len(T_END)) {
    update(b, X[t])
    secs[t] <- as.numeric(Sys.time() - t0, units = "secs")
  }
  secs
}

curves <- rbind(
  expand.grid(panel = "left",  variance = c("known", "unknown"),
              M = c(M_GRID, Inf), H = 0, stringsAsFactors = FALSE),
  expand.grid(panel = "right", variance = c("known", "unknown"),
              M = M_FIXED, H = H_GRID, stringsAsFactors = FALSE),
  expand.grid(panel = "right", variance = c("known", "unknown"),
              M = Inf, H = 13, stringsAsFactors = FALSE)
)

invisible(time_one_run("known", 0, 50))  
res <- list()
for (r in seq_len(N_REP)) {
  for (i in seq_len(nrow(curves))) {
    junk <- numeric(sample.int(30000L, 1)); rm(junk)   # jitter the GC phase
    secs <- time_one_run(curves$variance[i], curves$H[i], curves$M[i])
    res[[length(res) + 1]] <- data.frame(rep = r, panel = curves$panel[i],
                                         variance = curves$variance[i],
                                         M = curves$M[i], H = curves$H[i],
                                         d1 = 2 * curves$H[i] + 1,
                                         t = seq_len(T_END),
                                         us_per_step = diff(c(0, secs)) * 1e6)
  }
  if (r %% 10 == 0) cat(sprintf("rep %d/%d\n", r, N_REP))
}
timings <- bind_rows(res)
saveRDS(timings, file.path(OUT, "speed_timings.rds"))

summ <- timings %>%
  group_by(panel, variance, M, H, d1, t) %>%
  summarise(med = median(us_per_step),
            mn  = min(us_per_step),
            lo  = quantile(us_per_step, 0.1),
            hi  = quantile(us_per_step, 0.9), .groups = "drop")
write.csv(summ, file.path(OUT, "speed_summary.csv"), row.names = FALSE)

plat <- timings %>%
  filter(panel == "right", t >= 300) %>%
  group_by(variance, d1) %>% summarise(cost = median(us_per_step), .groups = "drop")
cat("\nplateau cost (us/step) at M =", M_FIXED, ", d2 = 1:\n")
print(as.data.frame(plat %>% mutate(cost = round(cost))), row.names = FALSE)

panel_labs <- c(left = "(a)", right = "(b)")
summ$panel_lab <- factor(panel_labs[summ$panel], panel_labs)

summ$lab <- ifelse(summ$panel == "left",
                   ifelse(is.finite(summ$M), paste0("M == ", summ$M), '"exact"'),
                   ifelse(is.finite(summ$M),
                          paste0("list(d[1] == ", summ$d1, ", d[2] == 1)"),
                          paste0('list("exact", d[1] == ', summ$d1, ', d[2] == 1)')))
lab_df <- filter(summ, variance == "known", t == max(t))

vlines <- rbind(
  data.frame(panel = "left",  x = M_GRID),
  data.frame(panel = "right", x = M_FIXED))
vlines$panel_lab <- factor(panel_labs[vlines$panel], panel_labs)

cols <- c(known = "#466CA6", unknown = "#D1AE45")  

timings=readRDS("simulation/speed_timings.rds")

library(ggrepel)
make_linear_figure <- function(TMAX, xbreaks, fname, keep = c("left", "right")) {
  s <- timings %>% filter(t <= TMAX, panel %in% keep) %>%
    group_by(panel, variance, M, H, d1, t) %>%
    summarise(med = median(us_per_step), .groups = "drop")
  s$panel_lab <- factor(panel_labs[s$panel], panel_labs)
  s$lab <- ifelse(s$panel == "left",
                  ifelse(is.finite(s$M), paste0("M == ", s$M), '"exact"'),
                  ifelse(is.finite(s$M),
                         paste0("list(d[1] == ", s$d1, ", d[2] == 1)"),
                         paste0('list("exact", d[1] == ', s$d1, ', d[2] == 1)')))
  lab_df <- filter(s, variance == "known", t == max(t))
  vl <- rbind(data.frame(panel = "left", x = M_GRID),
              data.frame(panel = "right", x = M_FIXED))
  vl <- vl[vl$panel %in% keep, ]
  vl$panel_lab <- factor(panel_labs[vl$panel], panel_labs)

  ggplot(s, aes(t, med, colour = variance, linetype = variance,
                group = interaction(variance, M, H))) +
    geom_vline(data = vl, aes(xintercept = x), linetype = "dotted", colour = "grey55") +
    geom_line(linewidth = 1.2) +
    geom_text_repel(data = lab_df,  aes(
    x = TMAX * 1.1,y = med,label = lab),hjust = -4,nudge_x = TMAX * 0.015,direction = "y",size = 3,colour = "grey20",parse = TRUE,show.legend = FALSE,
  segment.colour = NA,box.padding = 0.3,point.padding = 0,force = 1,max.overlaps = Inf,seed = 1) +
    facet_wrap(~panel_lab, strip.position = "bottom",
               labeller = as_labeller(function(x) rep("t", length(x)))) +
    scale_x_continuous(breaks = xbreaks, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous() +
    coord_cartesian(clip = "off") +
    scale_colour_manual(values = cols) +
    scale_linetype_manual(values = c(known = "solid", unknown = "dashed"))  +
    labs(x = NULL, y = expression("per-step runtime (" * mu * "s)")) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(size = 14),
          panel.spacing = unit(6, "lines"),
          plot.margin = margin(6, 100, 6, 6))
  ggsave(file.path(OUT, fname), width = if (length(keep) == 1) 8 else 10,
         height = 4.2, units = "in", device = cairo_pdf)
  cat("Written:", file.path(OUT, fname), "\n")
}

make_linear_figure(T_END, c(1, 250, 500, 750, 1000), "runtime_per_step_linlin.pdf")
make_linear_figure(400L,  c(1, 100, 200, 300, 400),  "runtime_per_step_linlin_T400.pdf")
make_linear_figure(T_END, c(1, 250, 500, 750, 1000), "runtime_per_step_linlin_left.pdf",
                   keep = "left")
