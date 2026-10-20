library(parallel)
library(bagelR)
library(dplyr)

RESULTS_DIR <- "simulation/results_seasonality"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

N <- 2000
TAU_TRUE <- 1000
PERIOD <- 24
HARMONICS <- 1
TRUE_PRE_CHANGE <- c(sin = 1, cos = 1)
TRUE_POST_CHANGE <- TRUE_PRE_CHANGE + rep(0.25, 2)
PRIOR_MEAN <- c(sin = 1, cos = 1)
SIGMA2 <- 1
N_REP <- 500
P0 <- 0.9
DECAY <- 1
N_CORES <- 3
REUSE_THRESH <- FALSE

PRIOR_SETTINGS <- data.frame(
  prior = c("main", "diffuse"),
  sin_variance = c(0.25^2, 10 * 0.25^2),
  cos_variance = c(0.25^2, 10 * 0.25^2),
  stringsAsFactors = FALSE
)
PARTICLE_COUNTS <- c(50, 100, 200)

load_saved_thresh <- function(out_file, pram) {
  if (!file.exists(out_file)) return(NA_real_)
  e <- new.env()
  load(out_file, envir = e)
  if (is.null(e$thresh) || is.null(e$pram)) return(NA_real_)
  if (!isTRUE(all.equal(as.list(e$pram), as.list(pram),
                        check.attributes = FALSE))) return(NA_real_)
  as.numeric(e$thresh)
}

outputs_table <- function(res) {
  df <- bind_rows(res)
  index_out <- which(df$stop <= TAU_TRUE | df$stop > N | is.infinite(df$stop))
  keep <- if (length(index_out)) -index_out else seq_len(nrow(df))
  data.frame(
    coverage_rate = mean(df$coverage[keep] == 1),
    coverage_rate_out = mean(df$coverage[keep] == 1),
    drop_out = length(index_out) / N_REP,
    failure_rate = sum(is.infinite(df$stop) | df$stop > N),
    false_alarm = sum(df$stop <= TAU_TRUE),
    mean_delay = mean(df$stop[keep] - TAU_TRUE),
    sd_delay = sd(df$stop[keep] - TAU_TRUE),
    mean_map = mean(df$map[keep]),
    sd_map = sd(df$map[keep]),
    mean_length = mean(df$length[keep]),
    sd_length = sd(df$length[keep])
  )
}

seasonal_row <- function(t) {
  feature <- numeric(2 * HARMONICS)
  for (harmonic in seq_len(HARMONICS)) {
    feature[2 * harmonic - 1] <- sin(2 * pi * harmonic * t / PERIOD)
    feature[2 * harmonic] <- cos(2 * pi * harmonic * t / PERIOD)
  }
  feature
}

seasonal_feature <- function(t, tau) {
  feature <- seasonal_row(t)
  if (tau == 0) return(matrix(feature, ncol = 1))
  change_feature <- if (t > tau) feature else numeric(length(feature))
  matrix(c(feature, change_feature), ncol = 1)
}

seasonal_transform <- function(tau) {
  dimension <- 2 * HARMONICS
  cbind(diag(dimension), diag(dimension))
}

make_seasonal_prior <- function(pre_mean, post_mean, prior_scale) {
  dimension <- 2 * HARMONICS
  pre_mean <- rep(pre_mean, length.out = dimension)
  post_mean <- rep(post_mean, length.out = dimension)
  scale_matrix <- diag(rep(prior_scale, length.out = dimension), dimension)
  no_change_prior <- list(mu = matrix(pre_mean, ncol = 1),
                          sigma = scale_matrix)
  change_prior <- list(
    mu = matrix(c(pre_mean, post_mean - pre_mean), ncol = 1),
    sigma = rbind(
      cbind(scale_matrix, -scale_matrix),
      cbind(-scale_matrix, 2 * scale_matrix)
    )
  )
  function(t) if (t == 1) no_change_prior else change_prior
}

generate_profile_data <- function(seed, changed) {
  set.seed(seed)
  time <- seq_len(N)
  design <- t(vapply(time, seasonal_row, numeric(2 * HARMONICS)))
  pre_profile <- as.numeric(design %*% TRUE_PRE_CHANGE)
  post_coefficients <- if (changed) TRUE_POST_CHANGE else TRUE_PRE_CHANGE
  post_profile <- as.numeric(design %*% post_coefficients)
  mean <- ifelse(time <= TAU_TRUE, pre_profile, post_profile)
  mean + rnorm(N, 0, sqrt(SIGMA2))
}

find_cover <- function(post, true_change) {
  if (!length(post) || sum(post) <= 0) {
    return(list(cover = 0, len = NA_integer_))
  }
  ord <- order(post, decreasing = TRUE)
  cut_idx <- which(cumsum(post[ord]) >= 0.95)[1]
  selected <- ord[seq_len(cut_idx)]
  list(cover = as.numeric((true_change + 1) %in% selected),
       len = length(selected))
}

posterior_from_fit <- function(fit, pruning) {
  if (pruning == "LP-P") {
    post <- numeric(length(weights_tau_zero(fit)))
    post[as.integer(taus(fit)) + 1] <- as.numeric(active_weights(fit))
  } else {
    post <- as.numeric(weights(fit))
  }
  post[!is.finite(post)] <- 0
  post[1] <- 0
  mass <- sum(post)
  if (mass > 0) post <- post / mass
  post
}

score_fit <- function(fit, pruning, threshold) {
  # Include threshold crossings at the final observation.
  scores <- 1 - as.numeric(weights_tau_zero(fit))
  crossings <- which(scores > threshold)
  stop_time <- if (length(crossings)) crossings[1] else Inf
  post <- posterior_from_fit(fit, pruning)
  map <- which.max(post[-1])
  coverage <- find_cover(post, TAU_TRUE)
  list(stop = stop_time, map = map,
       coverage = coverage$cover, length = coverage$len)
}

fit_model <- function(x, prior, max_particles, threshold, pruning) {
  bagel_kv_approximate(
    x, seasonal_feature, prior, seasonal_transform,
    P0, DECAY, sqrt(SIGMA2), as.integer(max_particles), threshold,
    pruning = pruning
  )
}

return_arl_bagel <- function(seed, prior, max_particles, pruning) {
  x <- generate_profile_data(seed, FALSE)
  fit <- fit_model(x, prior, max_particles, Inf, pruning)
  max(1 - as.numeric(weights_tau_zero(fit)), na.rm = TRUE)
}

run_simulation_bagel <- function(seed, prior, max_particles, threshold, pruning) {
  x <- generate_profile_data(seed, TRUE)
  score_fit(fit_model(x, prior, max_particles, threshold, pruning), pruning, threshold)
}

simulation_grid <- expand.grid(
  method = c("bagel", "LP-M", "LP-P"),
  K = PARTICLE_COUNTS,
  prior = PRIOR_SETTINGS$prior,
  stringsAsFactors = FALSE
)
simulation_grid <- merge(simulation_grid, PRIOR_SETTINGS, by = "prior", sort = FALSE)
simulation_grid <- rbind(
  simulation_grid,
  transform(PRIOR_SETTINGS, method = "bagel", K = 10000)[
    , c("prior", "method", "K", "sin_variance", "cos_variance")]
)
simulation_grid <- simulation_grid[
  , c("method", "K", "prior", "sin_variance", "cos_variance")]
row.names(simulation_grid) <- NULL
start_cluster <- function(n_cores) {
  cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")

  clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(bagelR))
    NULL
  })

  clusterExport(cl, c(
    "N", "TAU_TRUE", "PERIOD", "HARMONICS",
    "TRUE_PRE_CHANGE", "TRUE_POST_CHANGE",
    "SIGMA2", "P0", "DECAY",
    "seasonal_row", "seasonal_feature", "seasonal_transform",
    "generate_profile_data", "find_cover",
    "posterior_from_fit", "score_fit", "fit_model",
    "return_arl_bagel", "run_simulation_bagel"
  ), envir = .GlobalEnv)
  cl
}
run_example_seasonality <- function() {
  cl <- NULL

  on.exit({
    if (!is.null(cl)) {
      try(stopCluster(cl), silent = TRUE)
    }
  }, add = TRUE)
#  for (sim_no in seq_len(nrow(simulation_grid))) {
  for (sim_no in 1:20) {
    pram <- simulation_grid[sim_no, ]
    out_file <- file.path(RESULTS_DIR, sprintf("new_np%d.Rdata", sim_no))
    prior <- make_seasonal_prior(
      PRIOR_MEAN, PRIOR_MEAN,
      c(pram$sin_variance, pram$cos_variance)
    )
    pruning <- pram$method
    max_particles <- pram$K
    no_prune <- pruning == "bagel" && max_particles == 10000
    wanted_cores <- if (no_prune) 2 else N_CORES
    cl <- start_cluster(wanted_cores)

    thresh <- if (REUSE_THRESH) load_saved_thresh(out_file, pram) else NA_real_
    if (!is.finite(thresh)) {
      arl <- parLapplyLB(
        cl, seq_len(N_REP), return_arl_bagel,
        prior = prior, max_particles = max_particles, pruning = pruning
      )
      thresh <- as.numeric(quantile(unlist(arl), 0.95, na.rm = TRUE))
    }

    res <- parLapplyLB(
      cl, seq_len(N_REP), run_simulation_bagel,
      prior = prior, max_particles = max_particles,
      threshold = thresh, pruning = pruning
    )
    save(res, thresh, pram, sim_no, file = out_file)
    stopCluster(cl)
    gc()
  }

  present <- seq_len(nrow(simulation_grid))
  files <- file.path(RESULTS_DIR, sprintf("new_np%d.Rdata", present))
  present <- present[file.exists(files)]
  files <- files[file.exists(files)]
  summaries <- lapply(files, function(file) {
    load(file)
    outputs_table(res)
  })
  final <- bind_cols(simulation_grid[present, ], bind_rows(summaries))
  write.csv(final, file.path(RESULTS_DIR, "example_seasonality.csv"),
            row.names = FALSE)
  final
}

#run_example_seasonality()


source("simulation/example2_postprocessing.R")

postprocess_stopping_times(
  make_data = function(seed) {
    generate_profile_data(seed, changed = TRUE)
  },
  make_prior = function(p) {
    make_seasonal_prior(
      PRIOR_MEAN,
      PRIOR_MEAN,
      c(p$sin_variance, p$cos_variance)
    )
  },
  make_reference = function(prior, p) {
    sequential_bagel_kv_approximate(
      seasonal_feature, prior, seasonal_transform,
      P0, DECAY, sqrt(SIGMA2), N + 10,
      pruning = "bagel"
    )
  },
  make_approximation = function(x, stop_time, prior, p) {
    fit_model(
      x[seq_len(stop_time)], prior, p$K, Inf,
      as.character(p$method)
    )
  },
  output_name = "seasonality",
  n_cores = 2
)
