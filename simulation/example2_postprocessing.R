stopping_hpd_set <- function(post, level = 0.95) {
  post[!is.finite(post) | post < 0] <- 0
  mass <- sum(post)
  if (!is.finite(mass) || mass <= 0) return(integer(0))

  post <- post / mass
  ordering <- order(post, decreasing = TRUE)
  cutoff <- which(cumsum(post[ordering]) >= level)[1]
  ordering[seq_len(cutoff)]
}

stopping_reference_summary <- function(fit, true_change, level = 0.95) {
  post <- as.numeric(weights(fit))
  post[!is.finite(post) | post < 0] <- 0
  post[1] <- 0
  mass <- sum(post)

  if (!is.finite(mass) || mass <= 0) {
    return(list(
      posterior = numeric(0), hpd_set = integer(0),
      map = NA_real_, coverage = NA_real_, hpd_size = NA_real_
    ))
  }

  post <- post / mass
  selected <- stopping_hpd_set(post, level)
  list(
    posterior = post,
    hpd_set = selected,
    map = which.max(post[-1]),
    coverage = as.numeric((true_change + 1) %in% selected),
    hpd_size = length(selected)
  )
}

stopping_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

stopping_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1) sd(x) else NA_real_
}

stopping_global_dependencies <- function(functions,
                                         envir = .GlobalEnv) {
  pending <- functions
  inspected <- character(0)
  exports <- character(0)

  while (length(pending)) {
    current <- pending[[1]]
    pending <- pending[-1]
    if (!is.function(current)) next

    globals <- codetools::findGlobals(current, merge = FALSE)
    names_used <- unique(c(globals$functions, globals$variables))
    names_used <- setdiff(names_used, inspected)
    inspected <- c(inspected, names_used)

    local_names <- names_used[vapply(
      names_used,
      exists,
      logical(1),
      envir = envir,
      inherits = FALSE
    )]
    exports <- unique(c(exports, local_names))

    for (name in local_names) {
      value <- get(name, envir = envir, inherits = FALSE)
      if (is.function(value)) pending[[length(pending) + 1]] <- value
    }
  }

  exports
}

stopping_process_replicate <- function(replicate, state) {
  cells <- state$cells
  results <- lapply(cells, function(z) z$res[[replicate]])
  stops <- vapply(results, function(z) as.numeric(z$stop), numeric(1))
  valid <- is.finite(stops) & stops > state$tau_true & stops <= state$n
  snapshots <- list()
  x <- NULL

  if (any(valid)) {
    wanted_times <- sort(unique(as.integer(stops[valid])))
    x <- state$make_data(replicate)
    reference_fit <- state$make_reference(state$prior, state$reference_pram)

    for (t in seq_len(max(wanted_times))) {
#      print(t)
      update(reference_fit, x[t])
      if (t %in% wanted_times) {
        snapshots[[as.character(t)]] <- stopping_reference_summary(
          reference_fit, state$tau_true, state$level
        )
      }
    }
    rm(reference_fit)
    gc(verbose = FALSE)
  }

  rows <- vector("list", length(cells))
  for (j in seq_along(cells)) {
    cell <- cells[[j]]
    result <- results[[j]]
    method <- as.character(cell$pram$method)

    exact_map <- exact_coverage <- exact_hpd_size <- NA_real_
    recomputed_map <- recomputed_coverage <- recomputed_hpd_size <- NA_real_
    hpd_size_ratio <- reference_mass_in_approx_hpd <- NA_real_
    hpd_jaccard <- absolute_map_difference <- posterior_regret <- NA_real_

    if (valid[j]) {
      stop_time <- as.integer(stops[j])
      exact <- snapshots[[as.character(stop_time)]]
      no_pruning_row <-
        method == "bagel" && as.numeric(cell$pram$K) >= state$n
      if (no_pruning_row) {
        approximate_post <- exact$posterior
      } else {
        approximate_fit <- state$make_approximation(
          x, stop_time, state$prior, cell$pram
        )
        approximate_post <- state$posterior_from_fit(approximate_fit, method)
        rm(approximate_fit)
      }
      approximate_set <- stopping_hpd_set(approximate_post, state$level)

      recomputed_map <- which.max(approximate_post[-1])
      recomputed_coverage <-
        as.numeric((state$tau_true + 1) %in% approximate_set)
      recomputed_hpd_size <- length(approximate_set)
      exact_map <- exact$map
      exact_coverage <- exact$coverage
      exact_hpd_size <- exact$hpd_size
      hpd_size_ratio <- recomputed_hpd_size / exact_hpd_size

      # 1. Reference mass captured by the approximate HPD set.
      reference_mass_in_approx_hpd <-
        sum(exact$posterior[approximate_set])

      # 4. Jaccard overlap of approximate and reference HPD sets.
      union_set <- union(approximate_set, exact$hpd_set)
      hpd_jaccard <- if (length(union_set)) {
        length(intersect(approximate_set, exact$hpd_set)) /
          length(union_set)
      } else {
        NA_real_
      }

      # 6. Absolute displacement of the approximate MAP.
      absolute_map_difference <- abs(recomputed_map - exact_map)

      # 8. Reference-posterior regret of the approximate MAP.
      probability_at_approximate_map <-
        exact$posterior[recomputed_map + 1]
      posterior_regret <- log(
        max(exact$posterior) /
          max(probability_at_approximate_map, .Machine$double.xmin)
      )
      status <- "valid"
    } else {
      status <- if (is.finite(stops[j]) && stops[j] <= state$tau_true) {
        "false_alarm"
      } else {
        "failure"
      }
    }

    saved_map <- as.numeric(result$map)
    saved_coverage <- as.numeric(result$coverage)
    saved_hpd_size <- as.numeric(result$length)
    rows[[j]] <- data.frame(
      replicate = replicate,
      sim_no = cell$sim_no,
      method = method,
      K = as.numeric(cell$pram$K),
      prior = state$prior_name,
      stop = stops[j],
      status = status,
      delay = if (status == "valid") {
        stops[j] - state$tau_true
      } else {
        NA_real_
      },
      saved_approx_map = saved_map,
      recomputed_approx_map = recomputed_map,
      exact_map = exact_map,
      saved_approx_coverage = saved_coverage,
      recomputed_approx_coverage = recomputed_coverage,
      exact_coverage = exact_coverage,
      saved_approx_hpd_size = saved_hpd_size,
      recomputed_approx_hpd_size = recomputed_hpd_size,
      exact_hpd_size = exact_hpd_size,
      hpd_size_ratio = hpd_size_ratio,
      reference_mass_in_approx_hpd = reference_mass_in_approx_hpd,
      hpd_jaccard = hpd_jaccard,
      absolute_map_difference = absolute_map_difference,
      posterior_regret = posterior_regret,
      map_matches_saved = if (status == "valid") {
        recomputed_map == saved_map
      } else {
        NA
      },
      coverage_matches_saved = if (status == "valid") {
        recomputed_coverage == saved_coverage
      } else {
        NA
      },
      hpd_size_matches_saved = if (status == "valid") {
        recomputed_hpd_size == saved_hpd_size
      } else {
        NA
      }
    )
  }

  result <- do.call(rbind, rows)
  gc(verbose = FALSE)
  result
}

stopping_process_worker <- function(replicate) {
  stopping_process_replicate(replicate, .stopping_worker_state)
}

stopping_signature_value <- function(value) {
  if (is.function(value)) {
    return(list(
      type = "function",
      formals = paste(deparse(formals(value)), collapse = "\n"),
      body = paste(deparse(body(value)), collapse = "\n")
    ))
  }
  if (is.environment(value) || typeof(value) == "externalptr") {
    return(list(type = typeof(value)))
  }
  value
}

stopping_dependency_signature <- function(names, envir = .GlobalEnv) {
  names <- sort(unique(names))
  values <- lapply(names, function(name) {
    stopping_signature_value(get(name, envir = envir, inherits = FALSE))
  })
  names(values) <- names
  values
}

stopping_load_checkpoint <- function(path, metadata) {
  if (!file.exists(path)) return(NULL)
  checkpoint <- tryCatch(
    readRDS(path),
    error = function(error) {
      warning(
        "Could not read checkpoint ", path, ": ", conditionMessage(error),
        ". Recomputing this batch."
      )
      NULL
    }
  )
  if (is.null(checkpoint)) return(NULL)
  if (!is.list(checkpoint) ||
      !identical(checkpoint$metadata, metadata) ||
      !is.list(checkpoint$output) ||
      length(checkpoint$output) != length(metadata$replicate_ids)) {
    message("  Ignoring incompatible checkpoint: ", basename(path), ".")
    return(NULL)
  }
  checkpoint$output
}

stopping_save_checkpoint <- function(path, metadata, output) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(
    list(metadata = metadata, output = output),
    temporary_path,
    version = 2
  )
  if (!file.copy(temporary_path, path, overwrite = TRUE)) {
    stop("Could not save checkpoint: ", path)
  }
  invisible(path)
}

postprocess_stopping_times <- function(make_data, make_prior,
                                       make_reference, make_approximation,
                                       output_name, level = 0.95,
                                       n_cores = 2,
                                       batch_size = 25,
                                       resume = TRUE,
                                       checkpoint_dir = NULL) {
  files <- list.files(
    RESULTS_DIR,
    pattern = "^new_np[0-9]+\\.Rdata$",
    full.names = TRUE
  )
  if (!length(files)) stop("No new_np*.Rdata files found in ", RESULTS_DIR)

  sim_number <- as.integer(
    sub("^new_np([0-9]+)\\.Rdata$", "\\1", basename(files))
  )
  ordering <- order(sim_number)
  files <- files[ordering]
  sim_number <- sim_number[ordering]
  input_md5 <- unname(tools::md5sum(files))

  cells <- lapply(seq_along(files), function(j) {
    e <- new.env()
    load(files[j], envir = e)
    if (is.null(e$res) || is.null(e$pram)) {
      stop("Missing res or pram in ", files[j])
    }
    list(
      res = e$res, pram = e$pram, sim_no = sim_number[j],
      input_md5 = input_md5[j]
    )
  })

  prior_names <- unique(vapply(
    cells, function(z) as.character(z$pram$prior), character(1)
  ))

  n_cores <- suppressWarnings(as.integer(n_cores))
  if (!length(n_cores) || !is.finite(n_cores) || n_cores < 1) n_cores <- 1
  batch_size <- suppressWarnings(as.integer(batch_size))
  if (!length(batch_size) || !is.finite(batch_size) || batch_size < 1) {
    batch_size <- 25
  }

  if (is.null(checkpoint_dir)) {
    checkpoint_dir <- file.path(
      RESULTS_DIR, paste0(output_name, "_stopping_checkpoints")
    )
  }
  if (isTRUE(resume)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }

  dependency_names <- stopping_global_dependencies(
    list(
      make_data, make_prior, make_reference, make_approximation,
      get("posterior_from_fit", envir = .GlobalEnv, inherits = FALSE)
    )
  )
  dependency_signature <- stopping_dependency_signature(dependency_names)
  callback_signature <- lapply(
    list(
      make_data = make_data,
      make_prior = make_prior,
      make_reference = make_reference,
      make_approximation = make_approximation,
      posterior_from_fit = get(
        "posterior_from_fit", envir = .GlobalEnv, inherits = FALSE
      )
    ),
    stopping_signature_value
  )
  postprocessing_signature <- lapply(
    list(
      stopping_hpd_set = stopping_hpd_set,
      stopping_reference_summary = stopping_reference_summary,
      stopping_process_replicate = stopping_process_replicate
    ),
    stopping_signature_value
  )

  helper_names <- character(0)
  if (n_cores > 1) {
    helper_names <- c(
      "stopping_hpd_set", "stopping_reference_summary",
      "stopping_process_replicate", "stopping_process_worker"
    )
  }

  output <- list()

  for (prior_name in prior_names) {
    cell_ids <- which(vapply(
      cells,
      function(z) as.character(z$pram$prior) == prior_name,
      logical(1)
    ))
    n_rep <- min(vapply(
      cells[cell_ids], function(z) length(z$res), integer(1)
    ))
    reference_pram <- cells[[cell_ids[1]]]$pram
    prior <- make_prior(reference_pram)
    state <- list(
      cells = cells[cell_ids],
      prior_name = prior_name,
      prior = prior,
      reference_pram = reference_pram,
      make_data = make_data,
      make_reference = make_reference,
      make_approximation = make_approximation,
      posterior_from_fit = get(
        "posterior_from_fit", envir = .GlobalEnv, inherits = FALSE
      ),
      tau_true = TAU_TRUE,
      n = N,
      level = level
    )

    workers_used <- min(n_cores, n_rep)
    message(
      "Post-processing ", n_rep, " replicates for prior '", prior_name,
      "' using ", workers_used, " worker", if (workers_used == 1) "" else "s",
      if (workers_used == 1) "." else paste0(
        " in batches of at most ", batch_size, "."
      )
    )

    prior_dependencies <- stopping_global_dependencies(list(prior))
    all_dependencies <- unique(c(dependency_names, prior_dependencies))
    replicate_batches <- split(
      seq_len(n_rep),
      ceiling(seq_len(n_rep) / batch_size)
    )
    prior_output <- list()
    safe_prior_name <- gsub("[^A-Za-z0-9_.-]", "_", prior_name)
    checkpoint_base <- list(
      version = 1,
      output_name = output_name,
      prior_name = prior_name,
      level = level,
      tau_true = TAU_TRUE,
      n = N,
      bagelR_version = as.character(utils::packageVersion("bagelR")),
      callbacks = callback_signature,
      postprocessing = postprocessing_signature,
      dependencies = dependency_signature,
      cells = lapply(state$cells, function(cell) {
        list(
          sim_no = cell$sim_no,
          pram = cell$pram,
          input_md5 = cell$input_md5
        )
      })
    )

    for (batch in seq_along(replicate_batches)) {
      replicate_ids <- replicate_batches[[batch]]
      checkpoint_metadata <- c(
        checkpoint_base,
        list(replicate_ids = as.integer(replicate_ids))
      )
      checkpoint_file <- file.path(
        checkpoint_dir,
        sprintf(
          "%s_%s_rep%04d-%04d.rds",
          output_name, safe_prior_name,
          min(replicate_ids), max(replicate_ids)
        )
      )
      batch_output <- if (isTRUE(resume)) {
        stopping_load_checkpoint(checkpoint_file, checkpoint_metadata)
      } else {
        NULL
      }

      if (!is.null(batch_output)) {
        message(
          "  Batch ", batch, "/", length(replicate_batches),
          ": loaded replicates ", min(replicate_ids), "-",
          max(replicate_ids), " from checkpoint."
        )
      } else {
        batch_workers <- min(n_cores, length(replicate_ids))
        message(
          "  Batch ", batch, "/", length(replicate_batches),
          ": computing replicates ", min(replicate_ids), "-",
          max(replicate_ids), "."
        )

        if (n_cores == 1) {
          batch_output <- lapply(
            replicate_ids, stopping_process_replicate, state = state
          )
        } else {
        cl <- parallel::makeCluster(
          batch_workers, type = "PSOCK", outfile = ""
        )
        batch_output <- tryCatch({
          parallel::clusterEvalQ(cl, {
            suppressPackageStartupMessages(library(bagelR))
            NULL
          })
          parallel::clusterExport(
            cl, helper_names,
            envir = environment(postprocess_stopping_times)
          )
          if (length(all_dependencies)) {
            parallel::clusterExport(
              cl, all_dependencies, envir = .GlobalEnv
            )
          }

          state_environment <- list2env(
            list(.stopping_worker_state = state), parent = emptyenv()
          )
          parallel::clusterExport(
            cl, ".stopping_worker_state", envir = state_environment
          )
          parallel::parLapplyLB(
            cl, replicate_ids, stopping_process_worker
          )
        }, finally = {
          try(parallel::stopCluster(cl), silent = TRUE)
        })
        }

        if (isTRUE(resume)) {
          stopping_save_checkpoint(
            checkpoint_file, checkpoint_metadata, batch_output
          )
          message("    Saved checkpoint: ", basename(checkpoint_file), ".")
        }
      }
      prior_output <- c(prior_output, batch_output)
      gc(verbose = FALSE)
    }
    output <- c(output, prior_output)
  }

  output <- do.call(rbind, output)
  rownames(output) <- NULL
  write.csv(
    output,
    file.path(RESULTS_DIR, paste0(output_name, "_stopping_comparison.csv")),
    row.names = FALSE
  )

  valid_output <- output[output$status == "valid", , drop = FALSE]
  metric_names <- c(
    "delay", "saved_approx_coverage", "recomputed_approx_coverage",
    "exact_coverage", "recomputed_approx_hpd_size", "exact_hpd_size",
    "hpd_size_ratio", "reference_mass_in_approx_hpd", "hpd_jaccard",
    "absolute_map_difference", "posterior_regret"
  )
  group_key <- interaction(
    valid_output$method, valid_output$K, valid_output$prior,
    drop = TRUE, lex.order = TRUE
  )
  summary <- do.call(rbind, lapply(split(valid_output, group_key), function(d) {
    means <- vapply(d[metric_names], stopping_mean, numeric(1))
    sds <- vapply(d[metric_names], stopping_sd, numeric(1))
    data.frame(
      method = d$method[1], K = d$K[1], prior = d$prior[1],
      n_valid = nrow(d),
      as.list(setNames(means, paste0("mean_", names(means)))),
      as.list(setNames(sds, paste0("sd_", names(sds)))),
      row.names = NULL,
      check.names = FALSE
    )
  }))
  rownames(summary) <- NULL
  write.csv(
    summary,
    file.path(RESULTS_DIR, paste0(output_name, "_stopping_summary.csv")),
    row.names = FALSE
  )

  invisible(list(replicates = output, summary = summary))
}
