library(dplyr)
library(parallel)
library(bagelR)

RESULTS_DIR <- "simulation/results_calibration"

MODEL_NAME <- "slope_varknown_continuous"
TT <- 1000
NUM_REPLICATES <- 500
MIN_TAU <- 200
MAX_TAU <- 800
NO_CHANGE_PROBABILITY <- 0.9
DECAY <- 1
NOISE_SD <- 1
TIME_SCALE <- TT
PARTICLE_COUNTS <- c(50,100,200)
MAIN_BASELINE_VAR <- c(
  intercept = 0.25^2,
  slope = 0.25^2
)
MAIN_CHANGE_VARIANCE <- 0.25^2
PRIOR_SETTINGS <- list(
  main = list(
    baseline_var = MAIN_BASELINE_VAR,
    change_var = MAIN_CHANGE_VARIANCE
  ),
  diffuse = list(
    baseline_var = 10 * MAIN_BASELINE_VAR,
    change_var = 10 * MAIN_CHANGE_VARIANCE
  )
)
NUM_CORES <- 2
SKIP_EXISTING <- FALSE

particle_cap <- function(particle_count) {
  if (is.na(particle_count)) TT + 10 else as.integer(particle_count)
}

slope_feature <- function(t, tau) {
  normalized_time <- t / TIME_SCALE
  if (tau == 0) return(matrix(c(1, normalized_time), 2, 1))
  z <- matrix(0, 4, 1)
  z[1:2, 1] <- c(1, normalized_time)
  if (t > tau) z[3:4, 1] <- c(1, normalized_time)
  z
}

slope_transform <- function(tau) {
  matrix(c(1, 0, 1, 0, 0, 1, 0, 1), 2, 4, byrow = TRUE)
}

make_slope_prior <- function(d1, d2, d4) function(t) {
  if (t == 1)
    return(list(mu = matrix(0, 2, 1), sigma = diag(c(d1, d2))))
  tau <- (t - 1) / TIME_SCALE
  S <- matrix(0, 4, 4)
  S[1, 1] <- d1
  S[2, 2] <- d2
  S[3, 3] <- tau^2 * d4
  S[3, 4] <- S[4, 3] <- -tau * d4
  S[4, 4] <- d4
  list(mu = matrix(0, 4, 1), sigma = S)
}

normalise_posterior <- function(weights) {
  weights <- as.numeric(weights)
  total_mass <- sum(weights)
  change_mass <- sum(weights[-1])
  list(
    full = weights / total_mass,
    changepoint = weights[-1] / change_mass
  )
}

posterior_distance <- function(approximate, exact) {
  list(
    tv_full = 0.5 * sum(abs(approximate$full - exact$full)),
    tv_change = 0.5 * sum(abs(approximate$changepoint - exact$changepoint)),
    w1_change = sum(abs(cumsum(approximate$changepoint) -
                          cumsum(exact$changepoint)))
  )
}

fit_bagel <- function(data, prior, particle_count) {
  fit <- bagel_kv_approximate(
    data, slope_feature, prior, slope_transform,
    NO_CHANGE_PROBABILITY, DECAY, NOISE_SD,
    particle_cap(particle_count), Inf
  )
  normalise_posterior(weights(fit))
}

fit_reference <- function(data,prior, method,
                          particle_count) {
  fit <- bagel_kv_approximate(
    data, slope_feature, prior, slope_transform,
    NO_CHANGE_PROBABILITY, DECAY, NOISE_SD,
    particle_cap(particle_count), Inf,pruning = method
  )
  if (method == "LP-P") {
    recovered_weights <- numeric(TT)
    recovered_weights[taus(fit)+1] <- active_weights(fit)
    recovered_weights <- recovered_weights / sum(recovered_weights)
  } else {
    recovered_weights <- weights(fit)
  }
  normalise_posterior(recovered_weights)
}


result_row <- function(replicate_id, true_changepoint, method, particle_count,
                       posterior, exact_posterior) {
  distance <- posterior_distance(posterior, exact_posterior)
  data.frame(
    rep = replicate_id, tau_true = true_changepoint,
    method = method, M = particle_count,
    tv_full = distance$tv_full, tv_change = distance$tv_change,
    w1_change = distance$w1_change,
    stringsAsFactors = FALSE
  )
}

run_replicate <- function(replicate_id, baseline_var, change_var) {
  set.seed(replicate_id)
  true_changepoint <- sample(MIN_TAU:MAX_TAU, 1)
  intercept <- rnorm(1, 0, sqrt(baseline_var[1]))
  slope <- rnorm(1, 0, sqrt(baseline_var[2]))
  slope_change <- rnorm(1, 0, sqrt(change_var))
  time_index <- seq_len(TT)
  time <- time_index / TIME_SCALE
  changepoint_time <- true_changepoint / TIME_SCALE
  mean <- intercept + slope * time + ifelse(
    time_index > true_changepoint,
    slope_change * (time - changepoint_time), 0
  )
  observations <- mean + rnorm(TT, 0, NOISE_SD)

  prior <- make_slope_prior(baseline_var[1], baseline_var[2], change_var)
  exact <- fit_bagel(observations, prior, NA_integer_)
  bagel <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_bagel(observations, prior, particle_count)
  })
  probability_only <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_reference(
      observations, prior, "LP-M", particle_count
    )
  })
  top_m <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_reference(
      observations, prior, "LP-P", particle_count
    )
  })

  do.call(rbind, c(
    list(result_row(
      replicate_id, true_changepoint, "exact", NA_integer_, exact, exact
    )),
    Map(function(count, posterior) {
      result_row(replicate_id, true_changepoint, "bagel", count, posterior, exact)
    }, PARTICLE_COUNTS, bagel),
    Map(function(count, posterior) {
      result_row(
        replicate_id, true_changepoint, "prob_only", count, posterior, exact
      )
    }, PARTICLE_COUNTS, probability_only),
    Map(function(count, posterior) {
      result_row(replicate_id, true_changepoint, "topM", count, posterior, exact)
    }, PARTICLE_COUNTS, top_m)
  ))
}

summarise_calibration <- function(df, setting) {
  key <- paste(df$method, ifelse(is.na(df$M), "exact", df$M), sep = ":")
  groups <- split(df, key)
  rows <- lapply(groups, function(g) data.frame(
    model = MODEL_NAME, prior = setting, method = g$method[1], M = g$M[1],
    n_used = nrow(g), tv_full = mean(g$tv_full),
    tv_change = mean(g$tv_change), w1_change = mean(g$w1_change)
  ))
  ans <- do.call(rbind, rows)
  ans[, c("model", "prior", "method", "M", "n_used",
          "tv_full", "tv_change", "w1_change")]
}

make_worker_cluster <- function() {
  root <- normalizePath(".")
  library_paths <- .libPaths()
  cluster <- parallel::makeCluster(NUM_CORES, type = "PSOCK", outfile = "")
  parallel::clusterCall(cluster, function(root, library_paths) {
    .libPaths(library_paths)
    library(bagelR)
    setwd(root)
    invisible(NULL)
  }, root = root, library_paths = library_paths)
  parallel::clusterExport(cluster, varlist = c(
    "run_replicate", "fit_bagel", "fit_reference", "result_row",
    "slope_feature", "make_slope_prior", "slope_transform",
    "normalise_posterior", "posterior_distance", "particle_cap",
    "TT", "MIN_TAU", "MAX_TAU", "NO_CHANGE_PROBABILITY", "DECAY",
    "NOISE_SD", "TIME_SCALE", "PARTICLE_COUNTS"),
    envir = .GlobalEnv)
  cluster
}

run_replicates <- function(make_cluster, worker, ...) {
  cluster <- make_cluster()
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  parallel::parLapplyLB(cluster, seq_len(NUM_REPLICATES), worker, ...)
}

result_file <- function(prior_name) {
  file.path(
    RESULTS_DIR,
    sprintf("calib_slopeconti_sec31_prior10x_m200_ablations_%s.rds", prior_name)
  )
}

run_calibration <- function() {
  for (setting in names(PRIOR_SETTINGS)) {
    out_file <- result_file(setting)
    prior_setting <- PRIOR_SETTINGS[[setting]]
    if (SKIP_EXISTING && file.exists(out_file)) {
      next
    }
    cat(sprintf(
      "Prior '%s', baseline variances=(%g, %g), change variance=%g\n",
      setting,
      prior_setting$baseline_var[1], prior_setting$baseline_var[2],
      prior_setting$change_var
    ))
    reps <- run_replicates(
      make_worker_cluster,
      run_replicate,
      baseline_var = prior_setting$baseline_var,
      change_var = prior_setting$change_var
    )
    df <- do.call(rbind, reps)
    saveRDS(df, out_file)
  }
  aggregate_calibration()
}

aggregate_calibration <- function() {
  all_rows <- lapply(names(PRIOR_SETTINGS), function(setting) {
    summarise_calibration(readRDS(result_file(setting)), setting)
  })
  ans <- do.call(rbind, all_rows)
  write.csv(ans, file.path(RESULTS_DIR, "calibration_slope_continuous.csv"), row.names = FALSE)
  invisible(ans)
}

if (!exists("CALIB_NO_RUN") || !isTRUE(CALIB_NO_RUN)) run_calibration()
