library(parallel)
library(dplyr)
library(bagelR)
RESULTS_DIR <- "simulation/results_calibration"

MODEL_NAME <- "seasonal_varknown"
TT <- 1000
NUM_REPLICATES <- 500
MIN_TAU <- 200
MAX_TAU <- 800
NO_CHANGE_PROBABILITY <- 0.9
DECAY <- 1
NOISE_SD <- 1
NUM_CORES <- 2
SKIP_EXISTING <- TRUE

PERIOD <- 24
HARMONICS <- 1
PRIOR_MEAN <- c(sin = 0, cos = 0)
MAIN_PRIOR_SCALE <- c(sin = 0.25^2, cos = 0.25^2)
PRIOR_SETTINGS <- list(
  main = MAIN_PRIOR_SCALE,
  diffuse = 10 * MAIN_PRIOR_SCALE
)
PARTICLE_COUNTS <- c(200, 100, 50)

particle_cap <- function(particle_count) {
  if (is.na(particle_count)) TT + 10 else as.integer(particle_count)
}

seasonal_row <- function(time, period, harmonics) {
  feature <- numeric(2 * harmonics)
  for (harmonic in seq_len(harmonics)) {
    feature[2 * harmonic - 1] <- sin(2 * pi * harmonic * time / period)
    feature[2 * harmonic] <- cos(2 * pi * harmonic * time / period)
  }
  feature
}

generate_profile_data <- function(n, changepoint, period, harmonics,
                                  pre_change_coefficients,
                                  post_change_coefficients, variance = 1) {
  time <- seq_len(n)
  design <- t(vapply(
    time,
    seasonal_row,
    numeric(2 * harmonics),
    period = period,
    harmonics = harmonics
  ))
  pre_change_profile <- as.numeric(design %*% pre_change_coefficients)
  post_change_profile <- as.numeric(design %*% post_change_coefficients)
  mean <- ifelse(
    time <= changepoint,
    pre_change_profile,
    post_change_profile
  )
  mean + rnorm(n, 0, sqrt(variance))
}

seasonal_feature <- function(time, changepoint) {
  feature <- seasonal_row(time, PERIOD, HARMONICS)
  if (changepoint == 0) {
    return(matrix(feature, ncol = 1))
  }

  change_feature <- if (time > changepoint) {
    feature
  } else {
    rep(0, length(feature))
  }
  matrix(c(feature, change_feature), ncol = 1)
}

make_seasonal_prior <- function(pre_change_mean, post_change_mean, prior_scale) {
  dimension <- 2 * HARMONICS
  pre_change_mean <- rep(pre_change_mean, length.out = dimension)
  post_change_mean <- rep(post_change_mean, length.out = dimension)
  scale_matrix <- diag(rep(prior_scale, length.out = dimension), dimension)

  no_change_prior <- list(
    mu = matrix(pre_change_mean, ncol = 1),
    sigma = scale_matrix
  )
  change_prior <- list(
    mu = matrix(c(
      pre_change_mean,
      post_change_mean - pre_change_mean
    ), ncol = 1),
    sigma = rbind(
      cbind(scale_matrix, -scale_matrix),
      cbind(-scale_matrix, 2 * scale_matrix)
    )
  )

  function(time) {
    if (time == 1) no_change_prior else change_prior
  }
}

seasonal_transform <- function(changepoint) {
  dimension <- 2 * HARMONICS
  cbind(diag(dimension), diag(dimension))
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
    data,
    seasonal_feature,
    prior,
    seasonal_transform,
    NO_CHANGE_PROBABILITY,
    DECAY,
    NOISE_SD,
    particle_cap(particle_count),
    Inf
  )
  normalise_posterior(weights(fit))
}

fit_reference <- function(data, prior, particle_count,method) {
  fit <- bagel_kv_approximate(
    data,
    seasonal_feature,
    prior,
    seasonal_transform,
    NO_CHANGE_PROBABILITY,
    DECAY,
    NOISE_SD,
    particle_cap(particle_count),
    Inf,pruning = method
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

run_replicate <- function(replicate_id, prior_scale) {
  set.seed(replicate_id)
  true_changepoint <- sample(MIN_TAU:MAX_TAU, 1)
  pre_change_coefficients <- rnorm(
    2 * HARMONICS,
    PRIOR_MEAN,
    sqrt(NOISE_SD^2 * prior_scale)
  )
  post_change_coefficients <- rnorm(
    2 * HARMONICS,
    PRIOR_MEAN,
    sqrt(NOISE_SD^2 * prior_scale)
  )
  observations <- generate_profile_data(
    TT,
    true_changepoint,
    PERIOD,
    HARMONICS,
    pre_change_coefficients,
    post_change_coefficients,
    NOISE_SD^2
  )

  prior <- make_seasonal_prior(
    PRIOR_MEAN,
    PRIOR_MEAN,
    prior_scale
  )
  exact <- fit_bagel(observations, prior, NA_integer_)
  bagel <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_bagel(observations, prior, particle_count)
  })
  probability_only <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_reference(observations, prior, particle_count,method = "LP-M"
    )
  })
  top_m <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_reference(observations, prior, particle_count,method= "LP-P")
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

run_replicates <- function(make_cluster, worker, ...) {
  cluster <- make_cluster()
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  parallel::parLapplyLB(cluster, seq_len(NUM_REPLICATES), worker, ...)
}

make_worker_cluster <- function() {
  root <- normalizePath(".")
  library_paths <- .libPaths()
  cluster <- parallel::makeCluster(NUM_CORES, type = "PSOCK", outfile = "")
  parallel::clusterCall(cluster, function(root, library_paths) {
    .libPaths(library_paths)
    suppressMessages(library(bagelR))
    setwd(root)
    invisible(NULL)
  }, root = root, library_paths = library_paths)

  parallel::clusterExport(
    cluster,
    varlist = c(
      "run_replicate",
      "fit_bagel",
      "fit_reference",
      "result_row",
      "seasonal_row",
      "generate_profile_data",
      "seasonal_feature",
      "make_seasonal_prior",
      "seasonal_transform",
      "normalise_posterior",
      "posterior_distance",
      "particle_cap",
      "TT",
      "MIN_TAU",
      "MAX_TAU",
      "NO_CHANGE_PROBABILITY",
      "DECAY",
      "NOISE_SD",
      "PERIOD",
      "HARMONICS",
      "PRIOR_MEAN",
      "PARTICLE_COUNTS"
    ),
    envir = .GlobalEnv
  )

  cluster
}

result_file <- function(prior_name) file.path(
  RESULTS_DIR,
  sprintf("calib_seasonal_iid_profile_prior10x_m200_%s.Rdata", prior_name))

run_calibration <- function() {
  cat(sprintf("\n===== CALIBRATION %s: T=%d R=%d support={%d..%d} cores=%d =====\n",
              MODEL_NAME, TT, NUM_REPLICATES, MIN_TAU, MAX_TAU, NUM_CORES))
  for (setting in names(PRIOR_SETTINGS)) {
    prior_scale <- PRIOR_SETTINGS[[setting]]
    out_file <- result_file(setting)
    if (SKIP_EXISTING && file.exists(out_file)) {
      cat(sprintf("Setting '%s' (skip: exists)\n", setting)); next
    }
    cat(sprintf(
      "Setting '%s', seasonal variances=(%g, %g)\n",
      setting,
      prior_scale[1],
      prior_scale[2]
    ))
    res <- run_replicates(
      make_worker_cluster,
      run_replicate,
      prior_scale = prior_scale
    )
    save(res, setting, prior_scale, file = out_file)
    cat(sprintf("  saved %s\n", basename(out_file)))
  }
  aggregate_calibration()
}

aggregate_calibration <- function() {
  rows <- list()
  for (setting in names(PRIOR_SETTINGS)) {
    f <- result_file(setting)
    load(f)
    df <- bind_rows(res)
    agg <- df %>%
      group_by(method, M) %>%
      summarise(n_used  = n(),
                tv_full = mean(tv_full),
                tv_change = mean(tv_change),
                w1_change = mean(w1_change),
                .groups = "drop") %>%
      mutate(model = MODEL_NAME, prior = setting) %>%
      select(model, prior, method, M, n_used, tv_full, tv_change, w1_change) %>%
      arrange(method != "exact", method, desc(M))
    rows[[setting]] <- agg
  }
  final <- bind_rows(rows)
  out_csv <- file.path(RESULTS_DIR, sprintf("calibration_%s.csv", MODEL_NAME))
  write.csv(final, out_csv, row.names = FALSE)
  cat("\nWritten:", out_csv, "\n")
  print(as.data.frame(final), row.names = FALSE)
  invisible(final)
}

if (!exists("CALIB_NO_RUN") || !isTRUE(CALIB_NO_RUN)) {
  run_calibration()
}
