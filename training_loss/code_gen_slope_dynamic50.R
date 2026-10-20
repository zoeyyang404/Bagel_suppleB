source("training_loss/analysis_slope.R")

use_saved_parameters <- FALSE
parameter_file <- "training_loss/code_gen_slope_paper_parameters.txt"

output_dir <- "training_loss/code_gen_paper_results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

code_bins <- read.csv("training_loss/InferenceTrace_code_5min.csv")
#model is y = beta0 + beta1 * time + beta3 * sin(2*pi*time/period) + 
# beta4 * cos(2*pi*time/period)+gamma0+gamma1*(time-tau)+gamma3*sin(2*pi*time/period)+
# gamma4 * cos(2*pi*time/period) + epsilon
y <- code_bins$m_gen
bin_time <- as.POSIXct(code_bins$bin, tz = "UTC")
n <- length(y)
period <- 288
raw_time <- seq_len(n)
scaled_time <- (raw_time - 1) / period #the estimated coefficient will be per one day. normalize for numerical stability
angles <- 2 * pi * raw_time / period
SEAS <- cbind(sin(angles), cos(angles))
X_full <- cbind(1, scaled_time, SEAS)
train_idx <- seq_len(period)
monitor_idx <- (period + 1):n
p0 <- 1 - 3 / n
decay <- 1
K <- 100
calibration_reps <- 1000
alpha_horizon <- 0.10
block_size <- 50 #chosen by acf

if (!use_saved_parameters) {
  initial_prior <- list(
    beta_mean = rep(0, 4),
    beta_cov = diag(100, 4),
    nu = 2,
    iota = 1
  )
  # estimate the prior parameters
  training_posterior <- null_posterior_uv(
    y[train_idx], start = 1,
    beta_mean = initial_prior$beta_mean,
    beta_cov = initial_prior$beta_cov,
    nu0 = initial_prior$nu,
    iota0 = initial_prior$iota
  )
  # estimate the threshold
  beta_mean <- as.numeric(training_posterior$mu)
  beta_cov <- training_posterior$Sigma
  gamma_mean <- rep(0, 4)
  training_fit <- lm.fit(X_full[train_idx, ], y[train_idx])
  training_residuals <- training_fit$residuals - mean(training_fit$residuals)
  training_variance_mean <- training_posterior$iota / (training_posterior$nu - 1)
  seasonal_amplitude <- sqrt(sum(beta_mean[3:4]^2))
  gamma_cov <- diag(c(1, 1, rep(seasonal_amplitude^2 / (4 * training_variance_mean), 2)))
  residual_scale <- sqrt(mean(training_residuals^2))
  nu <- training_posterior$nu
  iota <- training_posterior$iota
  prior <- make_prior_mean(beta_mean, gamma_mean, beta_cov, gamma_cov)
  feature_vector <- make_feature_vector_mean(monitor_idx[1])
  transform <- make_transform_mean(monitor_idx[1])
  X_monitor <- X_full[monitor_idx, ]
  beta_cholesky <- chol(beta_cov)

  simulate_null_maxima <- function(repetitions, seed, label) {
    set.seed(seed)
    maxima <- numeric(repetitions)
    started <- Sys.time()
    for (b in seq_len(repetitions)) {
      sigma2_b <- 1 / rgamma(1, shape = nu, rate = iota)
      beta_b <- beta_mean + sqrt(sigma2_b) * as.numeric(t(beta_cholesky) %*% rnorm(4))
      y_b <- as.numeric(X_monitor %*% beta_b) +
        sqrt(sigma2_b) / residual_scale *
        bootstrap_residual_series(training_residuals, block_size, length(monitor_idx))
      fit <- bagel_uv_approximate(
        Y = y_b, feature_vector = feature_vector, prior = prior,
        transform = transform, p0 = p0, p = decay,
        nu = nu, iota = iota, max_particles = as.integer(K), threshold = Inf
      )
      maxima[b] <- max(-as.numeric(log_weights_tau_zero(fit)))
      if (b %% 10 == 0) {
        cat(sprintf("%s %d/%d: %.1f minutes\n", label, b, repetitions,
                    as.numeric(difftime(Sys.time(), started, units = "mins"))))
        flush.console()
      }
    }
    maxima
  }

  calibration_maxima <- simulate_null_maxima(calibration_reps, 11, "Calibration")
  saveRDS(calibration_maxima, file.path(output_dir, "calibration_maxima.rds"))
  log_threshold <- as.numeric(quantile(calibration_maxima, 1 - alpha_horizon))

  calibration_summary <- data.frame(
    log_threshold = log_threshold,
    probability_threshold = -expm1(-log_threshold),
    horizon = length(monitor_idx),
    target_alarm_probability = alpha_horizon,
    calibration_reps = calibration_reps
  )
  print(calibration_summary)
  write.csv(calibration_summary, file.path(output_dir, "calibration_summary.csv"),
            row.names = FALSE)

  parameters <- list(
    log_threshold = log_threshold,
    period = period,
    p0 = p0,
    decay = decay,
    beta_mean0 = beta_mean,
    gamma_mean0 = gamma_mean,
    beta_cov0 = beta_cov,
    gamma_cov0 = gamma_cov,
    nu0 = nu,
    iota0 = iota,
    K = K,
    block_size = block_size,
    calibration_reps = calibration_reps,
    calibration_seed = 11,
    alpha_horizon = alpha_horizon
  )
  dput(parameters, file = parameter_file)
}

parameters <- dget(parameter_file)
period <- parameters$period
scaled_time <- (raw_time - 1) / period
angles <- 2 * pi * raw_time / period
SEAS <- cbind(sin(angles), cos(angles))

detection_result <- run_detection_loop_dynamic(
  series = y, log_threshold = parameters$log_threshold,
  p0 = parameters$p0, decay = parameters$decay,
  beta_mean0 = parameters$beta_mean0, gamma_mean0 = parameters$gamma_mean0,
  beta_cov0 = parameters$beta_cov0, gamma_cov0 = parameters$gamma_cov0,
  nu0 = parameters$nu0, iota0 = parameters$iota0,
  K = parameters$K, start0 = period + 1
)

alarm_windows <- which(is.finite(detection_result$stop))
detections <- data.frame(
  change_index = detection_result$change[alarm_windows],
  detection_index = detection_result$stop[alarm_windows],
  change_time = bin_time[detection_result$change[alarm_windows]],
  detection_time = bin_time[detection_result$stop[alarm_windows]],
  estimated_detection_lag_hours = (detection_result$stop[alarm_windows] -
                   detection_result$change[alarm_windows]) * 5 / 60
)
print(detections)
write.csv(detections, file.path(output_dir, "detections.csv"), row.names = FALSE)

saveRDS(
  list(
    parameters = parameters,
    detection_result = detection_result,
    settings = list(period = parameters$period, p0 = parameters$p0,
                    decay = parameters$decay, K = parameters$K,
                    calibration_seed = parameters$calibration_seed,
                    block_size = parameters$block_size,
                    gamma_update = "posterior seasonal amplitude at each restart",
                    likelihood = "independent Gaussian",
                    calibration_noise = "scaled circular block-resampled residuals")
  ),
  file.path(output_dir, "results.rds")
)

pdf(file.path(output_dir, "figure_code_gen_final.pdf"), width = 11, height = 6.5)
plot_detection(detection_result, y, bin_time, period,
               ylab_text = "Mean generated tokens")
dev.off()

cat("Finished:", output_dir, "\n")
