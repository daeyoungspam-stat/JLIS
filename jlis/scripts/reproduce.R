#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
options <- list(mode = "figures", setting = "A1", profile = "publication", cores = "1", reps = "50")
i <- 1L
while (i <= length(args)) {
  if (!startsWith(args[i], "--") || i == length(args)) stop("Arguments require --name value pairs")
  key <- substring(args[i], 3L)
  if (!key %in% c("mode", "setting", "profile", "cores", "reps", "replicate", "outdir", "baseline"))
    stop("Unknown argument: ", key)
  options[[key]] <- args[i + 1L]
  i <- i + 2L
}
library(jlis)
setting <- toupper(options$setting)
outdir <- if (is.null(options$outdir)) file.path("results", setting) else options$outdir
reps <- as.integer(options$reps)
if (is.na(reps) || reps < 1 || reps > 50) stop("--reps must be between 1 and 50")
if (options$mode == "figures") {
  write_figures(reference_results(), if (is.null(options$outdir)) "results/reference_figures" else outdir)
} else if (options$mode == "run") {
  positions <- if (is.null(options$replicate)) seq_len(reps) else as.integer(options$replicate)
  run_simulation(setting, positions, outdir, cores = as.integer(options$cores), profile = options$profile)
} else if (options$mode == "combine") {
  result <- combine_results(outdir, setting, seq_len(reps), baseline = options$baseline)
  write_figures(result, file.path(outdir, "figures"))
} else if (options$mode == "legacy") {
  if (options$profile != "publication") stop("Legacy mode requires the publication profile")
  run_legacy_simulation(setting, outdir, as.integer(options$cores), reps)
} else stop("--mode must be figures, run, combine, or legacy")
