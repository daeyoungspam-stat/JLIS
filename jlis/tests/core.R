library(jlis)

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}

set.seed(991)
before <- .Random.seed
d <- simulate_jlis("A1", profile = "smoke")
stopifnot(identical(before, .Random.seed), length(d$X_list) == 4L,
          identical(dim(d$X_list[[1]]), c(40L, 36L)))
stopifnot(identical(d$Y_list, simulate_jlis("A1", profile = "smoke")$Y_list))
stopifnot(simulation_settings("D1")$turns[10] == 11L,
          simulation_settings("D2")$turns[1] == 51L,
          simulation_settings("C1")$control$cv_max_iter == 250L)
fixture <- readRDS(system.file("extdata", "dgp_reference.rds", package = "jlis"))
for (id in names(fixture)) {
  generated <- simulate_jlis(id)
  for (n in names(fixture[[id]])) {
    x <- generated[[n]]
    actual <- if (is.list(x)) lapply(x, function(v) as.numeric(v)[seq_len(min(length(v), 25L))])
    else as.numeric(x)[seq_len(min(length(x), 25L))]
    stopifnot(isTRUE(all.equal(actual, fixture[[id]][[n]], tolerance = 1e-10)))
  }
}
rm(generated)

fit <- fit_jlis(d$Y_list, d$X_list, d$Z_list, d$p, d$q,
                0.1, 0.1, 0.1, alpha = 0.1,
                control = list(max_iter = 40L, dykstra_max_iter = 100L))
stopifnot(all(is.finite(unlist(coef(fit)))), max(abs(coef(fit)$L)) <= 0.1 / 6 + 1e-8)
pred <- predict(fit, d$X_list, d$Z_list)
stopifnot(length(pred) == 4L, all(lengths(pred) == 40L))
expect_error(fit_jlis(d$Y_list, d$X_list, d$Z_list, 3L, 6L, 0.1, 0.1, 0.1))
expect_error(fit_jlis(d$Y_list, d$X_list, d$Z_list, 6L, 6L, -1, 0.1, 0.1))
expect_error(fit_jlis(d$Y_list, d$X_list, d$Z_list, 6L, 6L, 1, 1, 1,
                      control = list(max_iter = 0)))

set.seed(18)
X <- list(matrix(rnorm(24 * 12), 24, 12))
Z <- list(matrix(rnorm(24 * 2), 24, 2))
Y <- list(matrix(rnorm(24), 24, 1))
rectangular <- fit_jlis(Y, X, Z, 3, 4, 0.1, 0.1, 0.1,
                        control = list(max_iter = 20L))
stopifnot(identical(dim(coef(rectangular)$L), c(3L, 4L)))
M <- diag(c(3, 1, 0))
prox <- getFromNamespace(".alpha_linear_b", "jlis")$prox
out <- prox(M, tau = 0.5, alpha = 3, p = 3, q = 3)
stopifnot(isTRUE(all.equal(out$L, diag(c(1, 0.5, 0)), tolerance = 1e-6)))

d <- simulate_jlis("C1", profile = "smoke")
fit <- fit_jlis(d$Y_list, d$X_list, d$Z_list, 6, 6, 0.1, 0.1, 0.1,
                alpha = 1, family = "binomial",
                control = list(max_iter = 15L, pgd_max_inner = 15L))
p <- unlist(predict(fit, d$X_list, d$Z_list))
stopifnot(all(is.finite(p)), all(p >= 0 & p <= 1), max(abs(fit$L)) <= 1 / 6 + 1e-8)

reference <- reference_results()
stopifnot(nrow(reference$errors) == 9600L, nrow(reference$ranks) == 800L,
          nrow(reference$selection) == 400L,
          all(table(reference$errors$setting, reference$errors$method) == 200L),
          all(is.finite(reference$errors$B_error)), all(is.finite(reference$errors$gamma_error)))

for (setting in c("A1", "C1")) {
  directory <- tempfile(paste0("jlis-", setting))
  paths <- run_simulation(setting, 1L, directory, profile = "smoke")
  modified <- file.info(paths)$mtime
  stopifnot(file.exists(paths))
  run_simulation(setting, 1L, directory, profile = "smoke")
  stopifnot(identical(modified, file.info(paths)$mtime))
  combined <- combine_results(directory, setting, 1L)
  stopifnot(nrow(combined$errors) == 4L, nrow(combined$ranks) == 1L,
            file.exists(file.path(directory, "combined_results.RData")))
  expect_error(combine_results(directory, setting, 1:2))
  expect_error(combine_results(directory, setting, 1L, baseline = "reference"))
  expect_error(run_simulation(setting, 1L, directory, profile = "publication"))
  unlink(directory, recursive = TRUE)
}
cat("Core estimation and reproduction checks passed.\n")
