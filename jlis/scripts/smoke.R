library(jlis)
root <- "results/smoke"
for (setting in c("A1", "C1")) {
  directory <- file.path(root, setting)
  run_simulation(setting, replicates = 1L, outdir = directory, profile = "smoke", cores = 1L)
  result <- combine_results(directory, setting, replicates = 1L)
  write_figures(result, file.path(directory, "figures"))
}
