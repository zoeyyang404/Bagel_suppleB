library(paletteer)
# read the log version
bagel_dir <- "bagel_logversion/bagelR"
src <- normalizePath(file.path(bagel_dir, "src"), winslash = "/")

units <- c(
  "prior_type",
  "noise_type",
  "normal_cdf",
  "normal_density",
  "density_norm",
  "density_lst",
  "location_scale_t_density"
)

cpp <- c(
  "// [[Rcpp::plugins(cpp20)]]",
  "// [[Rcpp::depends(RcppEigen)]]",
  "// [[Rcpp::depends(BH)]]",
  sprintf('#include "%s/%s.cpp"', src, units),
  readLines(file.path(src, "bagelR.cpp"))
)

Sys.setenv(
  PKG_CXXFLAGS = sprintf(
    '-I"%s" -DNDEBUG -DEIGEN_NO_DEBUG -DEIGEN_NO_STATIC_ASSERT',
    src
  )
)

Rcpp::sourceCpp(code = paste(cpp, collapse = "\n"))

sys.source(
  file.path(bagel_dir, "R", "bagelR.R"),
  envir = globalenv()
)

# we start from period, and we have sort of restart process
# therefore the feature vector should be start + t - 1 because the bagel labels from 1.
make_feature_vector_mean <- function(start) {
  function(t, tau) {
    i <- start + t - 1
    base <- c(1, scaled_time[i], SEAS[i, ])
    if (tau == 0) {
      return(matrix(base, ncol = 1))
    }
    ti <- start + tau - 1
    if (t <= tau) {
      return(matrix(c(base, rep(0, length(base))), ncol = 1))
    }
    matrix(c(base, 1, scaled_time[i] - scaled_time[ti], SEAS[i, ]), ncol = 1)
  }
}

make_transform_mean <- function(start) {
  function(tau) {
    d <- 2 + ncol(SEAS)
    if (tau == 0) return(diag(d))
    ti <- start + tau - 1
    A <- diag(d)
    A[1, 2] <- -scaled_time[ti]
    cbind(diag(d), A)
  }
}
# we assume the prior and change size are independent
make_prior_mean <- function(beta_mean, gamma_mean, beta_cov, gamma_cov) {

  d_beta  <- length(beta_mean)
  d_gamma <- length(gamma_mean)
  d_full  <- d_beta + d_gamma

  prior_null <- list(
    mu = matrix(beta_mean, ncol = 1),
    sigma = beta_cov
  )

  prior_cov_full <- matrix(0, d_full, d_full)
  prior_cov_full[1:d_beta, 1:d_beta] <- beta_cov
  prior_cov_full[(d_beta + 1):d_full, (d_beta + 1):d_full] <- gamma_cov

  prior_change <- list(
    mu = matrix(c(beta_mean, gamma_mean), ncol = 1),
    sigma = prior_cov_full
  )

  function(t) {
    if (t == 1) prior_null else prior_change
  }
}

fit_bagel_mean <- function(y, start, p0, decay,
                           beta_mean, gamma_mean, beta_cov, gamma_cov,
                           nu0, iota0, threshold, K) {

  bagel_uv_approximate(
    Y = y,
    feature_vector = make_feature_vector_mean(start),
    prior = make_prior_mean(beta_mean, gamma_mean, beta_cov, gamma_cov),
    transform = make_transform_mean(start),
    p0 = p0 ,
    p = decay,
    nu = nu0,
    iota = iota0,
    max_particles = as.integer(K),
    threshold = threshold
  )
}

map_posterior_uv <- function(y, start, tau, detected,
                             beta_mean, gamma_mean, beta_cov, gamma_cov,
                             nu0, iota0) {

  feature_vector <- make_feature_vector_mean(start)
  d_beta <- length(beta_mean)
  d_full <- d_beta + length(gamma_mean)

  y_used <- y[1:detected]
  X <- do.call(rbind, lapply(1:detected, function(t) t(feature_vector(t, tau))))

  prior_mean <- matrix(c(beta_mean, gamma_mean), ncol = 1)
  prior_cov <- matrix(0, d_full, d_full)
  prior_cov[1:d_beta, 1:d_beta] <- beta_cov
  prior_cov[(d_beta + 1):d_full, (d_beta + 1):d_full] <- gamma_cov

  prior_precision <- chol2inv(chol(prior_cov))
  posterior_cov <- chol2inv(chol(prior_precision + crossprod(X)))
  posterior_mean <- posterior_cov %*%
    (prior_precision %*% prior_mean + crossprod(X, y_used))

  residual <- y_used - as.numeric(X %*% posterior_mean)
  mean_difference <- posterior_mean - prior_mean

  list(
    mu = posterior_mean,
    Sigma = posterior_cov,
    nu = nu0 + detected / 2,
    iota = iota0 + 0.5 * (sum(residual^2) +
      drop(crossprod(mean_difference, prior_precision %*% mean_difference)))
  )
}

null_posterior_uv <- function(y, start, beta_mean, beta_cov, nu0, iota0) {

  feature_vector <- make_feature_vector_mean(start)
  n_obs <- length(y)
  X <- do.call(rbind, lapply(1:n_obs, function(t) t(feature_vector(t, 0L))))

  prior_mean <- matrix(beta_mean, ncol = 1)
  prior_precision <- chol2inv(chol(beta_cov))
  posterior_cov <- chol2inv(chol(prior_precision + crossprod(X)))
  posterior_mean <- posterior_cov %*%
    (prior_precision %*% prior_mean + crossprod(X, y))

  residual <- y - as.numeric(X %*% posterior_mean)
  mean_difference <- posterior_mean - prior_mean

  list(
    mu = posterior_mean,
    Sigma = posterior_cov,
    nu = nu0 + n_obs / 2,
    iota = iota0 + 0.5 * (sum(residual^2) +
      drop(crossprod(mean_difference, prior_precision %*% mean_difference)))
  )
}

bootstrap_residual_series <- function(residual_clean, block_size, target_length) {
  n <- length(residual_clean)
  starts <- sample.int(n, ceiling(target_length / block_size), replace = TRUE)
  out <- unlist(
    lapply(
      starts,
      function(start) {
        residual_clean[((start - 1L + seq_len(block_size) - 1L) %% n) + 1L]
      }
    ),
    use.names = FALSE
  )
  out[1:target_length]
}

calibrate_variance_prior <- function(fit_season, REP_VARIANCE, block_size, seed) {

  fitted_clean <- as.numeric(fitted(fit_season))
  residual_clean <- as.numeric(residuals(fit_season))
  residual_clean <- residual_clean - mean(residual_clean)
  X_clean <- model.matrix(fit_season)
  n_clean <- length(residual_clean)
  initial_sigma2 <- summary(fit_season)$sigma^2

  set.seed(seed)
  sigma2_boot <- numeric(REP_VARIANCE)
  cf_boot <- matrix(0, REP_VARIANCE, ncol(X_clean))

  for (b in 1:REP_VARIANCE) {
    res_b <- bootstrap_residual_series(residual_clean, block_size, n_clean)
    fit_b <- lm.fit(X_clean, fitted_clean + res_b)
    sigma2_boot[b] <- sum(fit_b$residuals^2) / (n_clean - fit_b$rank)
    cf_boot[b, ] <- fit_b$coefficients
  }

  nu0 <- min(2 + initial_sigma2^2 / var(sigma2_boot), 2 + n_clean / 2)
  iota0 <- (nu0 - 1) * initial_sigma2
  cf_boot_cov <- cov(cf_boot)
  dimnames(cf_boot_cov) <- NULL

  list(
    nu0 = nu0,
    iota0 = iota0,
    initial_sigma2 = initial_sigma2,
    cf_boot_cov = cf_boot_cov,
    sigma2_boot = sigma2_boot,
    residual_clean = residual_clean,
    fitted_clean = fitted_clean,
    X_clean = X_clean
  )
}

calibrate_log_threshold <- function(vc, X_full, prior, p0, monitor_idx,
                                    decay, K, REP_THRESHOLD, block_size, seed) {

  n_clean <- length(vc$residual_clean)
  n_monitor <- length(monitor_idx)
  X_monitor <- X_full[monitor_idx, ]

  feature_vector_null <- make_feature_vector_mean(monitor_idx[1])
  transform_null <- make_transform_mean(monitor_idx[1])

  set.seed(seed)
  max_score <- numeric(REP_THRESHOLD)

  for (b in 1:REP_THRESHOLD) {

    res_train <- bootstrap_residual_series(vc$residual_clean, block_size, n_clean)
    cf_b <- lm.fit(vc$X_clean, vc$fitted_clean + res_train)$coefficients
    mean_b <- as.numeric(X_monitor %*% cf_b)

    res_monitor <- bootstrap_residual_series(vc$residual_clean, block_size, n_monitor)

    fit_b <- bagel_uv_approximate(
      Y = mean_b + res_monitor,
      feature_vector = feature_vector_null,
      prior = prior,
      transform = transform_null,
      p0 = p0,
      p = decay,
      nu = vc$nu0,
      iota = vc$iota0,
      max_particles = as.integer(K),
      threshold = Inf
    )

    max_score[b] <- max(-as.numeric(log_weights_tau_zero(fit_b)))
  }

  max_score
}


plot_detection <- function(detection_result, y, bin_time, period,
                           ylab_text = "Mean generated tokens",
                           xlab_text = "Time", hpd_level = 0.90) {

  par(mar = c(5, 5, 1.5, 5), cex.axis = 1.25, mgp = c(3, 0.9, 0))

  pal <- as.character(paletteer::paletteer_d("lisa::HilmaafKlint"))
  col_train <- pal[1]
  col_stat <- pal[2]
  col_fit <- pal[3]
  col_change <- pal[4]
  col_dot <- pal[5]

  posterior_nonempty <- Filter(
    function(x) is.data.frame(x) && nrow(x) > 0L,
    detection_result$tau_posterior
  )

  score_list <- lapply(detection_result$fit_history, function(fh) {
    list(index = fh$window_start:fh$window_stop,
         score = fh$score,
         alarm = identical(fh$model, "change"))
  })
  log_threshold <- detection_result$log_threshold
  score_max <- max(unlist(lapply(score_list, function(s) s$score)),
                   log_threshold)

  data_range <- range(y, finite = TRUE)
  data_span <- diff(data_range)
  band_height <- 0.16 * data_span
  gap <- 0.04 * data_span

  stat_top <- data_range[1] - gap
  stat_bottom <- stat_top - band_height
  plot_bottom <- stat_bottom - 0.02 * data_span
  score_offset <- stat_bottom + 0.05 * band_height
  score_scale <- 0.90 * band_height / score_max

  plot(
    bin_time, y, type = "n",
    ylim = c(plot_bottom, data_range[2]),
    axes = FALSE, xlab = "", ylab = ""
  )

  rect(bin_time[1], stat_bottom, bin_time[length(y)], stat_top,
       col = adjustcolor(col_stat, alpha.f = 0.10), border = NA)

  lines(bin_time, y, col = "black", lwd = 1.2)

  day_breaks   <- seq(trunc(min(bin_time), "days"), max(bin_time), by = "day")
  minor_breaks <- seq(trunc(min(bin_time), "days"), max(bin_time), by = "6 hours")
  axis.POSIXct(1, at = minor_breaks, labels = FALSE, tcl = -0.25)
  axis.POSIXct(1, at = day_breaks, format = "%b %d", tcl = -0.5)

  y_ticks <- pretty(data_range)
  y_ticks <- y_ticks[y_ticks >= data_range[1] & y_ticks <= data_range[2]]
  axis(2, at = y_ticks, labels = y_ticks, las = 1)

  score_ticks <- pretty(c(0, score_max), n = 3)
  score_ticks <- score_ticks[score_ticks >= 0 & score_ticks <= score_max]
  axis(4, at = score_offset + score_ticks * score_scale,
       labels = score_ticks, las = 1)

  box()
  mtext(xlab_text, side = 1, line = 3, cex = 1.45)
  mtext(ylab_text, side = 2, line = 3.2, cex = 1.45)

  det_times <- detection_result$stop[is.finite(detection_result$stop) &
                                        detection_result$stop < Inf]
  for (dt in det_times) {
    segments(bin_time[dt], stat_bottom, bin_time[dt], stat_top,
             col = "grey65", lty = 3, lwd = 1)
  }

  segments(
    x0 = bin_time[period + 1],
    y0 = score_offset + log_threshold * score_scale,
    x1 = bin_time[length(y)],
    y1 = score_offset + log_threshold * score_scale,
    lty = 2, lwd = 1.4, col = col_dot
  )
  for (s in score_list) {
    score_y <- score_offset + s$score * score_scale
    lines(bin_time[s$index], score_y, col = col_stat, lwd = 1.4)
    if (s$alarm) {
      points(bin_time[s$index[length(s$index)]],
             score_y[length(score_y)],
             pch = 16, cex = 0.9, col = col_stat)
    }
  }

  half_bin <- (as.numeric(bin_time[2]) - as.numeric(bin_time[1])) / 2
  for (posterior_j in posterior_nonempty) {
    probability <- posterior_j$probability_unconditional
    probability <- probability / sum(probability)
    ordered_bins <- order(probability, decreasing = TRUE)
    included_bins <- ordered_bins[seq_len(
      which(cumsum(probability[ordered_bins]) >= hpd_level)[1]
    )]
    bins <- sort(posterior_j$tau_global[included_bins])
    runs <- split(bins, cumsum(c(1L, diff(bins) != 1L)))
    for (run in runs) {
      rect(
        xleft = bin_time[min(run)] - half_bin,
        xright = bin_time[max(run)] + half_bin,
        ybottom = stat_top, ytop = par("usr")[4],
        col = adjustcolor(col_change, alpha.f = 0.40), border = NA
      )
    }
  }

  change_x <- bin_time[detection_result$change[is.finite(detection_result$change)]]
  segments(
    x0 = change_x, y0 = stat_top,
    x1 = change_x, y1 = par("usr")[4],
    lty = 2, lwd = 2, col = col_change
  )

  rect(
    xleft = bin_time[1], xright = bin_time[period],
    ybottom = par("usr")[3], ytop = par("usr")[4],
    col = adjustcolor(col_train, alpha.f = 0.45), border = NA
  )

  X <- cbind(1, scaled_time, SEAS)
  regime_start <- period + 1
  for (fit in detection_result$fit_history) {
    regime_stop <- if (identical(fit$model, "change")) fit$change_global else length(y)
    index <- regime_start:regime_stop
    beta <- fit$posterior_mu[seq_len(ncol(X))]
    lines(bin_time[index], as.numeric(X[index, ] %*% beta),
          lwd = 3, col = col_fit)
    regime_start <- regime_stop + 1
  }

}


run_detection_loop_dynamic <- function(series, log_threshold, p0, decay,
                                       beta_mean0, gamma_mean0, beta_cov0, gamma_cov0,
                                       nu0, iota0, K, start0 = 1) {
  stop_vec <- numeric(0)
  change <- numeric(0)
  start <- as.integer(start0)
  window <- 0
  n <- length(series)

  beta_mean <- beta_mean0
  beta_cov <- beta_cov0
  current_nu <- nu0
  current_iota <- iota0
  current_gamma_cov <- gamma_cov0

  fit_history <- list()
  tau_posterior_history <- list()

  while (start < n) {
    window <- window + 1
    y_window <- series[start:n]

    fit_full <- fit_bagel_mean(y_window, start, p0, decay,
                               beta_mean, gamma_mean0, beta_cov, current_gamma_cov,
                               current_nu, current_iota, threshold = Inf, K)

    log_score <- -as.numeric(log_weights_tau_zero(fit_full))
    crossing <- which(log_score > log_threshold)

    if (length(crossing) == 0) {
      posterior_null <- null_posterior_uv(y_window, start, beta_mean, beta_cov,
                                          current_nu, current_iota)
      fit_history[[window]] <- list(
        model = "no_change_final",
        window_start = start,
        window_stop = n,
        change_global = NA,
        score = log_score,
        posterior_mu = as.numeric(posterior_null$mu),
        posterior_Sigma = posterior_null$Sigma,
        posterior_nu = posterior_null$nu,
        posterior_iota = posterior_null$iota
      )
      stop_vec <- c(stop_vec, Inf)
      change <- c(change, NA)
      break
    }

    detected_local <- crossing[1]
    detected_global <- start + detected_local - 1

    fit_stop <- fit_bagel_mean(y_window[1:detected_local], start, p0,
                               decay, beta_mean, gamma_mean0, beta_cov,
                               current_gamma_cov, current_nu, current_iota,
                               threshold = Inf, K)

    w <- as.numeric(weights(fit_stop))
    w <- w / sum(w)
    tau_local <- seq_along(w)[-1] - 1
    w_pos <- w[-1]
    p_change <- sum(w_pos)

    tau_posterior_history[[window]] <- data.frame(
      window = window,
      window_start = start,
      detection_time = detected_global,
      tau_local = tau_local,
      tau_global = start + tau_local - 1,
      probability_unconditional = w_pos,
      probability_conditional = w_pos / p_change,
      null_probability = w[1],
      change_probability = p_change
    )

    map_tau <- tau_local[which.max(w_pos)]
    change_global <- start + map_tau - 1

    posterior <- map_posterior_uv(y_window, start, map_tau, detected_local,
                                  beta_mean, gamma_mean0, beta_cov, current_gamma_cov,
                                  current_nu, current_iota)

    transformation <- make_transform_mean(start)(map_tau)
    post_beta_mean <- as.numeric(transformation %*% posterior$mu)
    post_beta_cov <- transformation %*% posterior$Sigma %*% t(transformation)
    post_beta_cov <- (post_beta_cov + t(post_beta_cov)) / 2

    fit_history[[window]] <- list(
      model = "change",
      window_start = start,
      window_stop = detected_global,
      map_tau = map_tau,
      change_global = change_global,
      score = log_score[1:detected_local],
      posterior_mu = as.numeric(posterior$mu),
      posterior_Sigma = posterior$Sigma,
      posterior_nu = posterior$nu,
      posterior_iota = posterior$iota
    )

    cat(sprintf("window %d: change %d detected %d (delay %.1f h)\n",
                window, change_global, detected_global,
                (detected_global - change_global) * 5 / 60))

    beta_mean <- post_beta_mean
    beta_cov <- post_beta_cov
    current_nu <- posterior$nu
    current_iota <- posterior$iota
    seasonal_scale <- sum(beta_mean[3:4]^2) * (current_nu - 1) / (4 * current_iota)
    current_gamma_cov <- diag(c(1, 1, seasonal_scale, seasonal_scale))

    change <- c(change, change_global)
    stop_vec <- c(stop_vec, detected_global)
    start <- detected_global + 1
  }

  list(
    change = change,
    stop = stop_vec,
    tau_posterior = tau_posterior_history,
    fit_history = fit_history,
    log_threshold = log_threshold
  )
}
