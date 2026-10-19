run_legacy_simulation <- function(setting = "A1", outdir = file.path("results", setting, "legacy"),
                                  cores = 1L, reps = 50L) {
  setting <- simulation_settings(setting)$setting
  cores <- .positive_integer(cores, "cores")
  reps <- .positive_integer(reps, "reps")
  if (reps > 50L) stop("reps must be at most 50", call. = FALSE)
  needed <- c("foreach", "doParallel", "iterators", "glmnet", "MASS", "pracma")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install simulation dependencies: ", paste(missing, collapse = ", "), call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(outdir, pattern = "^simulation_results.*[.]RData$")))
    stop("Use a fresh legacy output directory", call. = FALSE)
  script <- system.file("legacy", paste0(setting, ".R"), package = "jlis", mustWork = TRUE)
  olddir <- getwd()
  variables <- c("SLURM_CPUS_PER_TASK", "JLIS_LEGACY_REPS")
  oldenv <- Sys.getenv(variables, unset = NA_character_)
  on.exit({
    setwd(olddir)
    for (i in seq_along(variables)) {
      if (is.na(oldenv[i])) Sys.unsetenv(variables[i])
      else do.call(Sys.setenv, setNames(list(oldenv[i]), variables[i]))
    }
  }, add = TRUE)
  Sys.setenv(SLURM_CPUS_PER_TASK = as.character(cores), JLIS_LEGACY_REPS = as.character(reps))
  setwd(outdir)
  e <- new.env(parent = .GlobalEnv)
  on.exit({ if (exists("cl", e, inherits = FALSE)) try(parallel::stopCluster(e$cl), silent = TRUE) }, add = TRUE)
  .with_seed(1L, function() sys.source(script, envir = e))
  files <- list.files(pattern = "^simulation_results.*[.]RData$", full.names = TRUE)
  if (!length(files)) stop("Legacy script completed without a result file", call. = FALSE)
  invisible(normalizePath(files))
}
