K_VALUES <- c(50, 100, 200)
EXACT_K <- 10000

scenario_files <- c(
  "\\shortstack[l]{Mean \\\\ change \\\\ known \\\\variance}" =
    "simulation/results_example1_varknown/mean_known_stopping_comparison.csv",
  "\\shortstack[l]{Mean\\\\ change \\\\ unknown \\\\variance}" =
    "simulation/results_example1_varunknown/mean_unknown_stopping_comparison.csv",
  "\\shortstack[l]{Continuous\\\\ slope \\\\ change\\\\ known \\\\ variance}" =
    "simulation/results_example2_conti/continuous_stopping_comparison.csv",
  "\\shortstack[l]{Discontinuous \\\\ slope\\\\ change \\\\ known \\\\ variance}" =
    "simulation/results_example2_dis/discontinuous_stopping_comparison.csv",
  "\\shortstack[l]{Seasonal \\\\ change \\\\ known \\\\ variance}" =
    "simulation/results_seasonality/seasonality_stopping_comparison.csv"
)
method_labels <- c(bagel = "Bagel", "LP-M" = "LP-M", "LP-P" = "LP-P")
prior_labels <- c(main = "Main prior", diffuse = "Diffuse prior")
settings <- expand.grid(
  K = K_VALUES,
  method = names(method_labels),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
settings <- settings[
  order(match(settings$K, K_VALUES), match(settings$method, names(method_labels))),
]
settings$label <- unname(method_labels[settings$method])
settings <- rbind(
  settings,
  data.frame(K = EXACT_K, method = "bagel", label = "Exact")
)

bold_best <- function(text, values, largest = FALSE) {
  best <- if (largest) max(values) else min(values)
  text[values == best] <- paste0("\\textbf{", text[values == best], "}")
  text
}

summary_rows <- list()
for (scenario in names(scenario_files)) {
  results <- read.csv(unname(scenario_files[scenario]))
  results <- results[
    results$K %in% c(K_VALUES, EXACT_K) &
      results$method %in% names(method_labels),
  ]

  for (prior in names(prior_labels)) {
    for (i in seq_len(nrow(settings))) {
      cell <- results[
        results$prior == prior &
          results$K == settings$K[i] &
          results$method == settings$method[i],
      ]
      valid <- cell[cell$status == "valid", ]

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        scenario = scenario,
        prior = prior,
        K = settings$K[i],
        method = settings$label[i],
        invalid = 1 - nrow(valid) / nrow(cell),
        delay = mean(valid$delay),
        map_distance = median(valid$absolute_map_difference),
        jaccard = mean(valid$hpd_jaccard),
        hpd_mass_error = 100 * mean(
          abs(valid$reference_mass_in_approx_hpd - 0.95)
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}
summary_table <- do.call(rbind, summary_rows)

latex_rows <- character(0)
for (scenario in names(scenario_files)) {
  main <- summary_table[
    summary_table$prior == "main" & summary_table$scenario == scenario,
  ]
  diffuse <- summary_table[
    summary_table$prior == "diffuse" & summary_table$scenario == scenario,
  ]

  main_exact <- main[main$K == EXACT_K & main$method == "Exact", ]
  diffuse_exact <- diffuse[diffuse$K == EXACT_K & diffuse$method == "Exact", ]
  latex_rows <- c(
    latex_rows,
    paste0(
      "\\multirow{10}{*}{", scenario, "} & -- & Exact & ",
      sprintf("%.2f", main_exact$invalid), " & ",
      sprintf("%.0f", main_exact$delay), " & ",
      sprintf("%.0f", main_exact$map_distance), " & ",
      sprintf("%.2f", main_exact$jaccard), " & ",
      sprintf("%.2f", main_exact$hpd_mass_error), " & ",
      sprintf("%.2f", diffuse_exact$invalid), " & ",
      sprintf("%.0f", diffuse_exact$delay), " & ",
      sprintf("%.0f", diffuse_exact$map_distance), " & ",
      sprintf("%.2f", diffuse_exact$jaccard), " & ",
      sprintf("%.2f", diffuse_exact$hpd_mass_error), " \\\\"
    ),
    "\\cline{2-13}"
  )

  for (K in K_VALUES) {
    main_K <- main[main$K == K, ]
    main_K <- main_K[match(unname(method_labels), main_K$method), ]
    diffuse_K <- diffuse[diffuse$K == K, ]
    diffuse_K <- diffuse_K[match(unname(method_labels), diffuse_K$method), ]

    main_map_distance <- bold_best(
      sprintf("%.0f", main_K$map_distance),
      main_K$map_distance
    )
    main_jaccard <- bold_best(
      sprintf("%.2f", main_K$jaccard), main_K$jaccard, largest = TRUE
    )
    main_mass_error <- bold_best(
      sprintf("%.2f", main_K$hpd_mass_error), main_K$hpd_mass_error
    )

    diffuse_map_distance <- bold_best(
      sprintf("%.0f", diffuse_K$map_distance),
      diffuse_K$map_distance
    )
    diffuse_jaccard <- bold_best(
      sprintf("%.2f", diffuse_K$jaccard), diffuse_K$jaccard, largest = TRUE
    )
    diffuse_mass_error <- bold_best(
      sprintf("%.2f", diffuse_K$hpd_mass_error), diffuse_K$hpd_mass_error
    )

    for (i in seq_len(nrow(main_K))) {
      K_cell <- if (i == 1) paste0("\\multirow{3}{*}{", K, "}") else ""
      latex_rows <- c(
        latex_rows,
        paste0(
          paste(
            c(
              "",
              K_cell,
              main_K$method[i],
              sprintf("%.2f", main_K$invalid[i]),
              sprintf("%.0f", main_K$delay[i]),
              main_map_distance[i],
              main_jaccard[i],
              main_mass_error[i],
              sprintf("%.2f", diffuse_K$invalid[i]),
              sprintf("%.0f", diffuse_K$delay[i]),
              diffuse_map_distance[i],
              diffuse_jaccard[i],
              diffuse_mass_error[i]
            ),
            collapse = " & "
          ),
          " \\\\"
        )
      )
    }
    if (K != K_VALUES[length(K_VALUES)]) {
      latex_rows <- c(latex_rows, "\\cline{2-13}")
    }
  }
  latex_rows <- c(latex_rows, "\\hline")
}

latex_table <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\renewcommand{\\arraystretch}{0.85}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{l|r|l|r|r|r|r|r|r|r|r|r|r}",
  "\\hline",
  paste0(
    " & & & ",
    "\\multicolumn{5}{c|}{Main prior} & ",
    "\\multicolumn{5}{c}{Diffuse prior} \\\\"
  ),
  "\\cline{4-13}",
  paste0(
    "Example & $K$ & Method & Inv. & Delay & MAP dist. & Jaccard & ",
    "HPD err & Inv. & Delay & MAP dist.  & Jaccard & ",
    "HPD err \\\\"
  ),
  "\\hline",
  latex_rows,
  "\\end{tabular}",
  "}",
  paste0(
    "\\caption{Online detection results under the main and diffuse priors ",
    "over 500 replications. ",
    "Inv. is the proportion of false alarms and detection failures. ",
    "MAP dist. is the median absolute distance between the approximate and ",
    "exact posterior modes, both evaluated at the method's stopping time. ",
    "HPD err is the mean absolute deviation of the exact posterior ",
    "mass in the approximate HPD set from 0.95, expressed in percentage ",
    "points.}"
  ),
  "\\label{tab:example2}",
  "\\end{table}"
)

output_file <- "simulation/example2_summary_table.tex"
writeLines(latex_table, output_file)
message("Written: ", output_file)
