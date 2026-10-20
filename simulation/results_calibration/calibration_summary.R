results_dir <- file.path("simulation", "results_calibration")

mean_known <- read.csv(file.path(results_dir, "calibration_mean_varknown.csv"))
mean_unknown <- read.csv(file.path(results_dir, "calibration_mean_varunknown.csv"))
slope_continuous <- read.csv(file.path(results_dir, "calibration_slope_continuous.csv"))
slope_discontinuous <- read.csv(file.path(results_dir, "calibration_slope_discontinuous.csv"))
seasonal <- read.csv(file.path(results_dir, "calibration_seasonal_varknown.csv"))

# Stack the five cases, retaining a concise label and a fixed display order.
case_data <- list(
  "Mean change known variance" = mean_known,
  "Mean change unknown variance" = mean_unknown,
  "Continuous slope change known variance" = slope_continuous,
  "Discontinuous slope change known variance" = slope_discontinuous,
  "Seasonal change known variance" = seasonal
)

calibration_long <- do.call(
  rbind,
  Map(
    function(dat, label) {
      dat$model_label <- label
      dat
    },
    case_data,
    names(case_data)
  )
)
row.names(calibration_long) <- NULL

# The exact rows are the reference and have no discrepancy values. Retain the
# three approximation methods and the two requested discrepancy measures.
calibration_long <- subset(
  calibration_long,
  method %in% c("bagel", "prob_only", "topM"),
  select = c(model_label, method, prior, M, tv_full, w1_change)
)

# Put the three methods next to one another, leaving prior as a row variable.
method_table <- function(method_name, prefix) {
  out <- subset(
    calibration_long,
    method == method_name,
    select = c(model_label, prior, M, tv_full, w1_change)
  )
  names(out)[4:5] <- paste0(prefix, "_", names(out)[4:5])
  out
}

calibration_table <- Reduce(
  function(x, y) merge(
    x, y, by = c("model_label", "prior", "M"), sort = FALSE
  ),
  list(
    method_table("bagel", "bagel"),
    method_table("prob_only", "lp_merge"),
    method_table("topM", "lp_prune")
  )
)

calibration_table$model_order <- match(
  calibration_table$model_label,
  names(case_data)
)
calibration_table$prior_order <- match(calibration_table$prior, c("main", "diffuse"))
calibration_table <- calibration_table[
  order(
    calibration_table$model_order,
    calibration_table$M,
    calibration_table$prior_order
  ),
  setdiff(names(calibration_table), c("model_order", "prior_order"))
]
row.names(calibration_table) <- NULL

# Put the two prior settings side by side so that each model-M combination
# occupies one row in the LaTeX table.
display_table <- unique(calibration_table[c("model_label", "M")])
row.names(display_table) <- NULL

# Write a compact LaTeX table. The document preamble needs:
#   \usepackage{booktabs,multirow,siunitx,graphicx}
fmt <- function(x, digits, bold = FALSE) {
  value <- sprintf(paste0("%.", digits, "f"), x)
  if (bold) paste0("\\textbf{", value, "}") else value
}
metric_cells <- function(metric, digits, i) {
  columns <- paste0(c("bagel", "lp_merge", "lp_prune"), "_", metric)
  values <- as.numeric(calibration_table[i, columns])
  is_minimum <- values == min(values, na.rm = TRUE)
  mapply(fmt, values, digits, is_minimum, USE.NAMES = FALSE)
}
model_labels <- c(
  "Mean change known variance" =
    "\\shortstack[l]{Mean change \\\\ known variance}",
  "Mean change unknown variance" =
    "\\shortstack[l]{Mean change \\\\ unknown \\\\variance}",
  "Continuous slope change known variance" =
    "\\shortstack[l]{Continuous slope \\\\change  known \\\\variance}",
  "Discontinuous slope change known variance" =
    "\\shortstack[l]{Discontinuous \\\\slope change \\\\known variance}",
  "Seasonal change known variance" =
    "\\shortstack[l]{Seasonal change \\\\ known variance}"
)
latex_rows <- character(nrow(display_table))
for (i in seq_len(nrow(display_table))) {
  first_in_model <- i == 1L ||
    display_table$model_label[i] != display_table$model_label[i - 1L]
  last_in_model <- i == nrow(display_table) ||
    display_table$model_label[i] != display_table$model_label[i + 1L]

  model_name <- display_table$model_label[i]
  particle_count <- display_table$M[i]
  main_i <- which(
    calibration_table$model_label == model_name &
      calibration_table$M == particle_count &
      calibration_table$prior == "main"
  )
  diffuse_i <- which(
    calibration_table$model_label == model_name &
      calibration_table$M == particle_count &
      calibration_table$prior == "diffuse"
  )
  stopifnot(length(main_i) == 1L, length(diffuse_i) == 1L)

  model_cell <- if (first_in_model) {
    paste0(
      "\\multirow{3}{*}{",
      unname(model_labels[model_name]),
      "}"
    )
  } else {
    ""
  }
  tv_main_cells <- metric_cells("tv_full", 3, main_i)
  tv_diffuse_cells <- metric_cells("tv_full", 3, diffuse_i)
  w1_main_cells <- metric_cells("w1_change", 2, main_i)
  w1_diffuse_cells <- metric_cells("w1_change", 2, diffuse_i)

  latex_rows[i] <- paste0(paste(c(
    model_cell,
    particle_count,
    tv_main_cells,
    tv_diffuse_cells,
    w1_main_cells,
    w1_diffuse_cells
  ), collapse = " & "), " \\\\", if (last_in_model) {
    " \\hline"
  } else {
    ""
  })
}

latex_table <- c(
  "\\begin{table}[]",
  "\\centering",
  "\\small",
  "\\scalebox{0.8}{%",
  paste0(
    "\\begin{tabular}{l|l|",
    paste(rep("r@{\\hspace{8pt}}r@{\\hspace{8pt}}r", 4), collapse = "|"),
    "}"
  ),
  "\\hline",
  paste0(
    "& & \\multicolumn{6}{c|}{$\\mathrm{TV}_{\\mathrm{full}}$} & ",
    "\\multicolumn{6}{c}{$W_{1,\\mathrm{change}}$} \\\\\\cline{3-14}"
  ),
  paste0(
    "& & \\multicolumn{3}{c|}{Main} & ",
    "\\multicolumn{3}{c|}{Diffuse} & ",
    "\\multicolumn{3}{c|}{Main} & ",
    "\\multicolumn{3}{c}{Diffuse} \\\\\\cline{3-14}"
  ),
  paste0(
    "Example & $M$ & Bagel & LP-M & LP-P & Bagel & LP-M & LP-P & ",
    "Bagel & LP-M & LP-P & Bagel & LP-M & LP-P \\\\"
  ),
  "\\hline",
  latex_rows,
  "\\end{tabular}",
  "}",
  paste0(
    "\\caption{Simulation on a fixed horizon $T=1000$ with $500$ replicates. W1 stands for the averaged Wasserstein distance and TV is for the averaged total variation distance. The lowest value in each block is in bold.}"
  ),
  "\\label{tab:calibration}",
  "\\end{table}"
)

output_file <- file.path(results_dir, "calibration_summary_table.tex")
writeLines(latex_table, output_file)
message("Written: ", output_file)
