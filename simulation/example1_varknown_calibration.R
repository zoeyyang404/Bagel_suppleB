library(parallel)
library(dplyr)
library(bagelR)

RESULTS_DIR <- "simulation/results_calibration"

MODEL_NAME <- "mean_varknown"
TT <- 1000
NUM_REPLICATES <- 500
MIN_TAU <- 200
MAX_TAU <- 800
NO_CHANGE_PROBABILITY <- 0.9
DECAY <- 1
NOISE_SD <- 1
NUM_CORES <- 2
SKIP_EXISTING <- TRUE

MAIN_PRIOR_VARIANCE <- 0.25^2
PRIOR_VARIANCES <- list(
  main = MAIN_PRIOR_VARIANCE,
  diffuse = 10 * MAIN_PRIOR_VARIANCE
)
PARTICLE_COUNTS <- c(200, 100, 50)

particle_cap <- function(particle_count) {
  if (is.na(particle_count)) TT + 10 else as.integer(particle_count)
}

mean_feature <- function(time, changepoint) {
  if (changepoint == 0) {
    return(matrix(1, 1, 1))
  }

  feature <- matrix(0, 2, 1)
  feature[1, 1] <- 1
  if (time > changepoint) {
    feature[2, 1] <- 1
  }
  feature
}

make_mean_prior <- function(pre_change_mean, post_change_mean, variance) {
  force(pre_change_mean)
  force(post_change_mean)
  force(variance)
  function(time) {
    if (time == 1) {
      return(list(
        mu = matrix(pre_change_mean, 1, 1),
        sigma = matrix(variance, 1, 1)
      ))
    }

    list(
      mu = matrix(c(pre_change_mean, post_change_mean - pre_change_mean), 2, 1),
      sigma = variance * matrix(c(1, -1, -1, 2), 2, 2, byrow = TRUE)
    )
  }
}

mean_transform <- function(changepoint) {
  matrix(c(1, 1), 1, 2)
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
  approximate_cdf <- cumsum(approximate$changepoint)
  exact_cdf <- cumsum(exact$changepoint)

  list(
    tv_full = 0.5 * sum(abs(approximate$full - exact$full)),
    tv_change = 0.5 * sum(abs(approximate$changepoint - exact$changepoint)),
    w1_change = sum(abs(approximate_cdf - exact_cdf))
  )
}

fit_bagel <- function(data, prior, particle_count) {
  fit <- bagel_kv_exact(
    data,
    mean_feature,
    prior,
    mean_transform,
    NO_CHANGE_PROBABILITY,
    DECAY,
    NOISE_SD,
    particle_cap(particle_count),
    Inf
  )
  normalise_posterior(weights(fit))
}

fit_reference <- function(data, prior, method, particle_count, true_changepoint) {
  fit <- bagel_kv_exact(
    data,
    mean_feature,
    prior,
    mean_transform,
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
    rep = replicate_id,
    tau_true = true_changepoint,
    method = method,
    M = particle_count,
    tv_full = distance$tv_full,
    tv_change = distance$tv_change,
    w1_change = distance$w1_change,
    stringsAsFactors = FALSE
  )
}

run_replicate <- function(replicate_id, prior_variance) {
  set.seed(replicate_id)

  true_changepoint <- sample(MIN_TAU:MAX_TAU, 1L)
  pre_change_mean <- rnorm(1, 0, sqrt(prior_variance))
  post_change_mean <- rnorm(1, 0, sqrt(prior_variance))
  observations <- rnorm(TT, 0, NOISE_SD) + c(
    rep(pre_change_mean, true_changepoint),
    rep(post_change_mean, TT - true_changepoint)
  )

  prior <- make_mean_prior(0, 0, prior_variance)
  exact <- fit_bagel(observations, prior, NA_integer_)
  bagel <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_bagel(observations, prior, particle_count)
  })

  probability_only <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_reference(
      observations,
      prior,
      "LP-M",
      particle_count,
      true_changepoint
    )
  })
  top_m <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_reference(
      observations,
      prior,
      "LP-P",
      particle_count,
      true_changepoint
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

run_replicates <- function(make_cluster, worker, ...) {
  cluster <- make_cluster()
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  parallel::parLapplyLB(
    cluster,
    seq_len(NUM_REPLICATES),
    worker,
    ...
  )
}

result_file <- function(prior_name) {
  file.path(
    RESULTS_DIR,
    sprintf(
      "calib_kv_sec31_prior10x_m200_ablations_%s.Rdata",
      prior_name
    )
  )
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

  parallel::clusterExport(
    cluster,
    varlist = c(
      "run_replicate",
      "fit_bagel",
      "fit_reference",
      "result_row",
      "mean_feature",
      "make_mean_prior",
      "mean_transform",
      "normalise_posterior",
      "posterior_distance",
      "particle_cap",
      "TT",
      "MIN_TAU",
      "MAX_TAU",
      "NO_CHANGE_PROBABILITY",
      "DECAY",
      "NOISE_SD",
      "PARTICLE_COUNTS"
    ),
    envir = .GlobalEnv
  )

  cluster
}

aggregate_calibration <- function() {
  summaries <- list()

  for (prior_name in names(PRIOR_VARIANCES)) {
    file <- result_file(prior_name)
    load(file)
    summaries[[prior_name]] <- bind_rows(res) %>%
      group_by(method, M) %>%
      summarise(
        n_used = n(),
        tv_full = mean(tv_full),
        tv_change = mean(tv_change),
        w1_change = mean(w1_change),
        .groups = "drop"
      ) %>%
      mutate(
        model = MODEL_NAME,
        prior = prior_name,
      ) %>%
      select(
        model,
        prior,
        method,
        M,
        n_used,
        tv_full,
        tv_change,
        w1_change,
      ) %>%
      arrange(method != "exact", method, desc(M))
  }

  final <- bind_rows(summaries)
  output_file <- file.path(
    RESULTS_DIR,
    sprintf("calibration_%s.csv", MODEL_NAME)
  )
  write.csv(final, output_file, row.names = FALSE)
  cat("\nWritten:", output_file, "\n")
  print(as.data.frame(final), row.names = FALSE)
  invisible(final)
}

run_calibration <- function() {
  cat(sprintf(
    "\n===== CALIBRATION %s: T=%d R=%d support={%d..%d} cores=%d =====\n",
    MODEL_NAME,
    TT,
    NUM_REPLICATES,
    MIN_TAU,
    MAX_TAU,
    NUM_CORES
  ))

  for (prior_name in names(PRIOR_VARIANCES)) {
    prior_variance <- PRIOR_VARIANCES[[prior_name]]
    output_file <- result_file(prior_name)

    if (SKIP_EXISTING && file.exists(output_file)) {
      cat(sprintf("Prior '%s' already exists; skipping.\n", prior_name))
      next
    }

    cat(sprintf("Prior '%s', variance=%g\n", prior_name, prior_variance))
    res <- run_replicates(
      make_worker_cluster,
      run_replicate,
      prior_variance = prior_variance
    )
    setting <- prior_name
    tau0 <- prior_variance
    save(res, setting, tau0, file = output_file)
    cat(sprintf("Saved %s\n", basename(output_file)))
  }

  aggregate_calibration()
}

if (!exists("CALIB_NO_RUN") || !isTRUE(CALIB_NO_RUN)) {
  run_calibration()
}
