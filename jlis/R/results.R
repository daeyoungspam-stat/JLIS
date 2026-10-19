.method_objects <- list(
  JLIS = c("Final_results_B", "Final_results_gamma"),
  `JLIS-alpha` = c("Final_results_B.alpha", "Final_results_gamma.alpha"),
  ILIS = c("Final_results_B.sep", "Final_results_gamma.sep"),
  IL = c("Final_results_B_onlylow", "Final_results_gamma_onlylow"),
  Lasso = c("Final_results_B_lasso", "Final_results_gamma_lasso"),
  Ridge = c("Final_results_B_ridge", "Final_results_gamma_ridge"))

reference_results <- function(setting = NULL) {
  path <- system.file("extdata", "reference_results.rds", package = "jlis", mustWork = TRUE)
  out <- readRDS(path)
  if (!is.null(setting)) {
    ids <- vapply(setting, function(id) simulation_settings(id)$setting, character(1))
    out <- lapply(out, function(x) x[x$setting %in% ids, , drop = FALSE])
  }
  out
}

.legacy_results <- function(path, config, positions) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  rows <- config$turns[positions]
  errors <- lapply(setdiff(names(.method_objects), "JLIS-alpha"), function(method) {
    vars <- .method_objects[[method]]
    for (v in vars) if (!exists(v, e, inherits = FALSE) || !is.matrix(e[[v]]) ||
        nrow(e[[v]]) < max(rows) || ncol(e[[v]]) != config$G ||
        any(!is.finite(e[[v]][rows, , drop = FALSE])))
      stop("Invalid baseline result: ", v, call. = FALSE)
    do.call(rbind, lapply(seq_len(config$G), function(g) data.frame(
      setting = config$setting, replicate = positions, turn = rows, group = g,
      method = method, B_error = e[[vars[1]]][rows, g], gamma_error = e[[vars[2]]][rows, g])))
  })
  if (!exists("Final_results_rank", e, inherits = FALSE) ||
      length(e$Final_results_rank) < max(rows) || any(!is.finite(e$Final_results_rank[rows])))
    stop("Invalid baseline ranks", call. = FALSE)
  list(errors = do.call(rbind, errors), ranks = data.frame(setting = config$setting,
    replicate = positions, turn = rows, method = "JLIS", rank = e$Final_results_rank[rows]))
}

combine_results <- function(outdir, setting = "A1", replicates = 1:50, baseline = NULL) {
  setting <- simulation_settings(setting)$setting
  if (!is.numeric(replicates) || !length(replicates) || anyNA(replicates) ||
      any(!replicates %in% 1:50) || anyDuplicated(replicates)) stop("Invalid replicate positions", call. = FALSE)
  replicates <- as.integer(replicates)
  files <- file.path(outdir, "parts", sprintf("rep_%03d.RData", replicates))
  if (any(!file.exists(files))) stop("Missing results: ", paste(basename(files[!file.exists(files)]), collapse = ", "), call. = FALSE)
  parts <- lapply(files, .load_part)
  config <- parts[[1]]$config
  if (!identical(config$setting, setting)) stop("Setting does not match result files", call. = FALSE)
  for (i in seq_along(parts)) .check_part(parts[[i]], config, replicates[i])
  r <- lapply(parts, `[[`, "result")
  errors <- do.call(rbind, lapply(r, function(x) data.frame(setting = setting,
    replicate = x$replicate, turn = x$turn, group = seq_len(config$G),
    method = "JLIS-alpha", B_error = x$B_error, gamma_error = x$gamma_error)))
  ranks <- do.call(rbind, lapply(r, function(x) data.frame(setting = setting,
    replicate = x$replicate, turn = x$turn, method = "JLIS-alpha", rank = x$rank)))
  selection <- do.call(rbind, lapply(r, function(x) data.frame(setting = setting,
    replicate = x$replicate, turn = x$turn, as.list(x$selected_values), row.names = NULL)))
  if (!is.null(baseline)) {
    if (!identical(config$profile, "publication")) stop("Baseline comparison requires the publication profile", call. = FALSE)
    if (identical(baseline, "reference")) {
      b <- reference_results(setting)
      b$errors <- b$errors[b$errors$method != "JLIS-alpha" & b$errors$replicate %in% replicates, ]
      b$ranks <- b$ranks[b$ranks$method != "JLIS-alpha" & b$ranks$replicate %in% replicates, ]
    } else {
      if (!is.character(baseline) || length(baseline) != 1L || !file.exists(baseline))
        stop("baseline must be 'reference' or an existing RData file", call. = FALSE)
      b <- .legacy_results(baseline, config, replicates)
    }
    errors <- rbind(b$errors, errors)
    ranks <- rbind(b$ranks, ranks)
  }
  rownames(errors) <- rownames(ranks) <- rownames(selection) <- NULL
  diagnostics <- do.call(rbind, lapply(r, function(x) data.frame(
    replicate = x$replicate, turn = x$turn, as.list(x$final_diagnostics))))
  results <- list(errors = errors, ranks = ranks, selection = selection,
                  config = config, diagnostics = diagnostics)
  e <- list2env(list(results = results), parent = emptyenv())
  for (method in unique(errors$method)) {
    d <- errors[errors$method == method, ]
    for (metric in c("B_error", "gamma_error")) {
      mat <- matrix(NA_real_, length(replicates), config$G,
                    dimnames = list(as.character(replicates), paste0("Group", seq_len(config$G))))
      mat[cbind(match(d$replicate, replicates), d$group)] <- d[[metric]]
      assign(.method_objects[[method]][if (metric == "B_error") 1L else 2L], mat, e)
    }
  }
  for (method in unique(ranks$method)) {
    d <- ranks[ranks$method == method, ]
    rank <- setNames(d$rank[match(replicates, d$replicate)], replicates)
    assign(if (method == "JLIS") "Final_results_rank" else "Final_results_rank.alpha", rank, e)
  }
  e$Final_main_values.alpha <- as.matrix(selection[, c("lambda1", "lambda2", "lambda3", "alpha")])
  e$Final_fit_diagnostics.alpha <- diagnostics
  save(list = ls(e), envir = e, file = file.path(outdir, "combined_results.RData"), version = 3)
  utils::write.csv(errors, file.path(outdir, "coefficient_errors.csv"), row.names = FALSE)
  utils::write.csv(summarize_results(results), file.path(outdir, "summary.csv"), row.names = FALSE)
  results
}

summarize_results <- function(results = reference_results()) {
  data <- results$errors
  groups <- split(data, interaction(data$setting, data$method, data$group, drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(d) data.frame(
    setting = d$setting[1], method = d$method[1], group = d$group[1], n = nrow(d),
    B_mean = mean(d$B_error), B_sd = stats::sd(d$B_error),
    gamma_mean = mean(d$gamma_error), gamma_sd = stats::sd(d$gamma_error))))
  rownames(out) <- NULL
  out
}

plot_results <- function(results = reference_results(), setting = "A1") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 to create plots", call. = FALSE)
  setting <- simulation_settings(setting)$setting
  d <- results$errors[results$errors$setting == setting, , drop = FALSE]
  if (!nrow(d)) stop("Requested setting has zero result rows", call. = FALSE)
  d <- rbind(data.frame(d[c("group", "method")], metric = "B error", error = d$B_error),
             data.frame(d[c("group", "method")], metric = "Gamma error", error = d$gamma_error))
  d$method <- factor(d$method, levels = names(.method_objects))
  d$group <- factor(paste("Group", d$group), levels = paste("Group", 1:4))
  d$metric <- factor(d$metric, levels = c("B error", "Gamma error"))
  colors <- c(JLIS = "#000000", `JLIS-alpha` = "#E69F00", ILIS = "#009E73",
              IL = "#56B4E9", Lasso = "#0072B2", Ridge = "#CC79A7")
  # Local bindings support NSE checks for the plotting expressions.
  method <- error <- group <- metric <- NULL
  ggplot2::ggplot(d, ggplot2::aes(x = method, y = error, fill = method)) +
    ggplot2::geom_boxplot(width = 0.7, outlier.size = 0.45, linewidth = 0.25) +
    ggplot2::facet_grid(metric ~ group, scales = "free_y") +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::labs(title = paste("Setting", setting), x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", linewidth = 0.25),
      plot.title = ggplot2::element_text(size = 12))
}

write_figures <- function(results = reference_results(), directory = "results/figures") {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- character()
  for (setting in unique(results$errors$setting)) {
    plot <- plot_results(results, setting)
    for (ext in c("pdf", "png")) {
      path <- file.path(directory, paste0("JLIS_", setting, ".", ext))
      ggplot2::ggsave(path, plot, width = 12, height = 5.5, units = "in", dpi = 300)
      paths <- c(paths, path)
    }
  }
  utils::write.csv(summarize_results(results), file.path(directory, "summary.csv"), row.names = FALSE)
  invisible(paths)
}
