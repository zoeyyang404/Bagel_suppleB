library(parallel)
library(dplyr)
library(bagelR)

RESULTS_DIR <- "simulation/results_calibration"

MODEL_NAME <- "mean_varunknown"
TT <- 1000
NUM_REPLICATES <- 500
MIN_TAU <- 200
MAX_TAU <- 800
NO_CHANGE_PROBABILITY <- 0.9
DECAY <- 1
NUM_CORES <- 2
SKIP_EXISTING <- TRUE

MAIN_CVAR <- 0.25^2
PRIOR_SETTINGS <- list(
  main    = list(cvar = MAIN_CVAR, nu = 6, iota = 5),
  diffuse = list(cvar = 10 * MAIN_CVAR, nu = 6, iota = 5)
)
PARTICLE_COUNTS <- c(200, 100, 50)

particle_cap <- function(particle_count) {
  if (is.na(particle_count)) TT + 10 else as.integer(particle_count)
}

mean_feature <- function(time, changepoint) {
  if (changepoint == 0) return(matrix(1, 1, 1))
  M <- matrix(0, 2, 1); M[1, 1] <- 1
  if (time > changepoint) M[2, 1] <- 1
  M
}
make_mean_prior <- function(pre_change_mean, post_change_mean, variance) {
  force(pre_change_mean); force(post_change_mean); force(variance)
  function(time) {
    if (time == 1)
      return(list(mu = matrix(pre_change_mean, 1, 1), sigma = matrix(variance, 1, 1)))
    list(mu    = matrix(c(pre_change_mean, post_change_mean - pre_change_mean), 2, 1),
         sigma = variance * matrix(c(1, -1, -1, 2), 2, 2, byrow = TRUE))
  }
}
mean_transform <- function(changepoint) matrix(c(1, 1), 1, 2, byrow = TRUE)



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

fit_bagel <- function(data, prior, particle_count, nu, iota) {
  fit <- bagel_uv_approximate(
    data,
    mean_feature,
    prior,
    mean_transform,
    NO_CHANGE_PROBABILITY,
    DECAY,
    nu,
    iota,
    particle_cap(particle_count),
    Inf
  )
  normalise_posterior(weights(fit))
}

fit_reference <- function(data, prior, nu, iota, method,
                          particle_count, true_changepoint) {
  fit <- bagel_uv_approximate(
    data,
    mean_feature,
    prior,
    mean_transform,
    NO_CHANGE_PROBABILITY,
    DECAY,
    nu,
    iota,
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

run_replicate <- function(replicate_id, prior_variance, nu, iota) {
  set.seed(replicate_id)
  true_changepoint <- sample(MIN_TAU:MAX_TAU, 1)
  noise_variance <- 1 / rgamma(1, shape = nu, rate = iota)
  pre_change_mean <- rnorm(1, 0, sqrt(prior_variance * noise_variance))
  post_change_mean <- rnorm(1, 0, sqrt(prior_variance * noise_variance))
  observations <- rnorm(TT, 0, sqrt(noise_variance)) + c(
    rep(pre_change_mean, true_changepoint),
    rep(post_change_mean, TT - true_changepoint)
  )

  prior <- make_mean_prior(0, 0, prior_variance)
  exact <- fit_bagel(observations, prior, NA_integer_, nu, iota)
  bagel <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_bagel(observations, prior, particle_count, nu, iota)
  })
  probability_only <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_reference(
      observations, prior, nu, iota, "LP-M",
      particle_count, true_changepoint
    )
  })
  top_m <- lapply(PARTICLE_COUNTS, function(particle_count) {
    fit_reference(
      observations, prior, nu, iota, "LP-P",
      particle_count, true_changepoint
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

  parallel::parLapplyLB(cluster, seq_len(NUM_REPLICATES), worker, ...)
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
    "mean_feature", "make_mean_prior", "mean_transform",
    "normalise_posterior", "posterior_distance", "particle_cap",
    "TT", "MIN_TAU", "MAX_TAU", "NO_CHANGE_PROBABILITY", "DECAY",
    "PARTICLE_COUNTS"),
    envir = .GlobalEnv)
  cluster
}

result_file <- function(prior_name) file.path(
  RESULTS_DIR,
  sprintf("calib_uv_sec31_ig6_5_prior10x_m200_ablations_%s.Rdata", prior_name)
)

run_calibration <- function() {
  cat(sprintf("\n===== CALIBRATION %s: T=%d R=%d support={%d..%d} cores=%d =====\n",
              MODEL_NAME, TT, NUM_REPLICATES, MIN_TAU, MAX_TAU, NUM_CORES))
  for (setting in names(PRIOR_SETTINGS)) {
    ps       <- PRIOR_SETTINGS[[setting]]
    out_file <- result_file(setting)
    if (SKIP_EXISTING && file.exists(out_file)) {
      cat(sprintf("Setting '%s' (skip: exists)\n", setting)); next
    }
    cat(sprintf("Setting '%s'  c=%g nu=%g iota=%g  ...\n",
                setting, ps$cvar, ps$nu, ps$iota))
    res <- run_replicates(
      make_worker_cluster,
      run_replicate,
      prior_variance = ps$cvar,
      nu = ps$nu,
      iota = ps$iota
    )
    ps_saved <- ps
    save(res, setting, ps_saved, file = out_file)
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
