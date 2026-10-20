library(parallel)
library(bagelR)
library(dplyr)

RESULTS_DIR <- "simulation/results_example1_varunknown"
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

N <- 2000
TAU_TRUE <- 1000
DELTA <- 0.25
SIGMA2 <- 1
N_REP <- 500
P0 <- 0.9
DECAY <- 1
N_CORES <- 3
REUSE_THRESH <- FALSE

PRIOR_SETTINGS <- data.frame(
  prior = c("main", "diffuse"),
  prior_variance = c(0.25^2, 10 * 0.25^2),
  nu = 6,
  iota = 5,
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

feature_vector_mean <- function(t, tau) {
  if (tau == 0) return(matrix(1, 1, 1))
  feature <- matrix(c(1, 0), 2, 1)
  if (t > tau) feature[2, 1] <- 1
  feature
}

make_prior_mean <- function(pre_mean, post_mean, variance) {
  force(pre_mean)
  force(post_mean)
  force(variance)
  function(t) {
    if (t == 1) {
      return(list(mu = matrix(pre_mean, 1, 1),
                  sigma = matrix(variance, 1, 1)))
    }
    list(
      mu = matrix(c(pre_mean, post_mean - pre_mean), 2, 1),
      sigma = variance * matrix(c(1, -1, -1, 2), 2, 2, byrow = TRUE)
    )
  }
}

transform_mean <- function(tau) matrix(c(1, 1), 1, 2)

find_cover <- function(post, true_change) {
  if (!length(post) || sum(post) <= 0) {
    return(list(cover = 0, len = NA_integer_))
  }
  ord <- order(post, decreasing = TRUE)
  cut_idx <- which(cumsum(post[ord]) >= 0.95)[1]
  selected <- ord[seq_len(cut_idx)]
  list(
    cover = as.numeric((true_change + 1) %in% selected),
    len = length(selected)
  )
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

fit_model <- function(x, prior, nu, iota, max_particles, threshold, pruning) {
  bagel_uv_approximate(
    x, feature_vector_mean, prior, transform_mean,
    P0, DECAY, nu, iota, as.integer(max_particles), threshold,
    pruning = pruning
  )
}

return_arl_bagel <- function(seed, prior, nu, iota, max_particles, pruning) {
  set.seed(seed)
  x <- rnorm(N, 0, sqrt(SIGMA2))
  fit <- fit_model(x, prior, nu, iota, max_particles, Inf, pruning)
  max(1 - as.numeric(weights_tau_zero(fit)), na.rm = TRUE)
}

run_simulation_bagel <- function(seed, prior, nu, iota, max_particles,
                                 threshold, pruning) {
  set.seed(seed)
  x <- rnorm(N, 0, sqrt(SIGMA2))
  x[(TAU_TRUE + 1):N] <- x[(TAU_TRUE + 1):N] + DELTA
  score_fit(
    fit_model(x, prior, nu, iota, max_particles, threshold, pruning),
    pruning, threshold
  )
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
    , c("prior", "method", "K", "prior_variance", "nu", "iota")]
)
simulation_grid <- simulation_grid[
  , c("method", "K", "prior", "prior_variance", "nu", "iota")]
row.names(simulation_grid) <- NULL

start_cluster <- function(n_cores) {
  cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")

  clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(bagelR))
    NULL
  })

  clusterExport(cl, c(
    "N", "TAU_TRUE", "DELTA", "SIGMA2", "P0", "DECAY",
    "feature_vector_mean", "transform_mean", "find_cover",
    "posterior_from_fit", "score_fit", "fit_model",
    "return_arl_bagel", "run_simulation_bagel"
  ), envir = .GlobalEnv)

  cl
}
run_example1_varunknown <- function() {
  cl <- start_cluster(N_CORES)
  on.exit(stopCluster(cl), add = TRUE)

  for (sim_no in seq_len(nrow(simulation_grid))) {
    pram <- simulation_grid[sim_no, ]
    out_file <- file.path(RESULTS_DIR, sprintf("new_np%d.Rdata", sim_no))
    prior <- make_prior_mean(0, 0, pram$prior_variance)
    pruning <- pram$method
    max_particles <- pram$K
    no_prune <- pruning == "bagel" && max_particles == 10000
    wanted_cores <- if (no_prune) 1 else N_CORES

    if (length(cl) != wanted_cores) {
      stopCluster(cl)
      cl <- start_cluster(wanted_cores)
    }
    thresh <- if (REUSE_THRESH) load_saved_thresh(out_file, pram) else NA_real_
    if (!is.finite(thresh)) {
      arl <- parLapplyLB(
        cl, seq_len(N_REP), return_arl_bagel,
        prior = prior, nu = pram$nu, iota = pram$iota,
        max_particles = max_particles, pruning = pruning
      )
      thresh <- as.numeric(quantile(unlist(arl), 0.95, na.rm = TRUE))
    }

    res <- parLapplyLB(
      cl, seq_len(N_REP), run_simulation_bagel,
      prior = prior, nu = pram$nu, iota = pram$iota,
      max_particles = max_particles, threshold = thresh, pruning = pruning
    )
    save(res, thresh, pram, sim_no, file = out_file)
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
  write.csv(final, file.path(RESULTS_DIR, "example1_varunknown_new.csv"),
            row.names = FALSE)
  final
}

#run_example1_varunknown()





source("simulation/example2_postprocessing.R")

postprocess_stopping_times(
  make_data = function(seed) {
    set.seed(seed)
    x <- rnorm(N, 0, sqrt(SIGMA2))
    x[(TAU_TRUE + 1):N] <- x[(TAU_TRUE + 1):N] + DELTA
    x
  },
  make_prior = function(p) {
    make_prior_mean(0, 0, p$prior_variance)
  },
  make_reference = function(prior, p) {
    sequential_bagel_uv_approximate(
      feature_vector_mean, prior, transform_mean,
      P0, DECAY, p$nu, p$iota, N + 10,
      pruning = "bagel"
    )
  },
  make_approximation = function(x, stop_time, prior, p) {
    fit_model(
      x[seq_len(stop_time)], prior, p$nu, p$iota,
      p$K, Inf, as.character(p$method)
    )
  },
  output_name = "mean_unknown",
  n_cores = 2
)
