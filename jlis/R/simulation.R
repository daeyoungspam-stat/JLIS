.simulation_engine <- function(config) {
  if (config$family == "binomial") .alpha_logistic else
    if (config$alpha_svd == "standard") .alpha_linear_a else .alpha_linear_b
}

.simulation_fit <- function(data, config, lambda, alpha, final = FALSE,
                            init = NULL, steps = NULL) {
  ctl <- config$control
  args <- list(Y_list = data$Y_list, X_list = data$X_list, Z_list = data$Z_list,
    px = config$p, qx = config$q, lambda1 = lambda[1], lambda2 = lambda[2],
    lambda3 = lambda[3], alpha = alpha, rho = 1,
    max_iter = if (final) ctl$final_max_iter else ctl$cv_max_iter, tol = ctl$tol,
    dykstra_max_iter = ctl$dykstra_max_iter, dykstra_tol = ctl$dykstra_tol, init = init)
  if (config$family == "binomial") args <- c(args, list(
    pgd_max_inner = ctl$pgd_max_inner, pgd_tol = ctl$pgd_tol, pgd_t_init = steps))
  do.call(.simulation_engine(config)$fit, args)
}

.fold_data <- function(data, fold) {
  config <- data$config
  train <- test <- list()
  for (name in c("X_list", "Z_list", "Y_list")) {
    test[[name]] <- train[[name]] <- vector("list", config$G)
    for (g in seq_len(config$G)) {
      idx <- seq.int(config$n_per_fold[g] * (fold - 1L) + 1L, config$n_per_fold[g] * fold)
      train[[name]][[g]] <- data[[name]][[g]][-idx, , drop = FALSE]
      test[[name]][[g]] <- data[[name]][[g]][idx, , drop = FALSE]
    }
  }
  list(train = train, test = test)
}

.validation_loss <- function(fit, data, family) {
  losses <- vapply(seq_along(data$Y_list), function(g) {
    if (family == "gaussian") return(sum((data$Y_list[[g]] -
      data$X_list[[g]] %*% fit$B[[g]] - data$Z_list[[g]] %*% fit$gamma[[g]])^2))
    eta <- as.numeric(data$X_list[[g]] %*% fit$B[[g]] + data$Z_list[[g]] %*% fit$gamma[[g]])
    y <- as.numeric(data$Y_list[[g]])
    prob <- pmin(pmax(.alpha_logistic$prob(eta), 1e-15), 1 - 1e-15)
    -mean(y * log(prob) + (1 - y) * log1p(-prob))
  }, numeric(1))
  if (family == "gaussian") sum(losses) / length(losses) else mean(losses)
}

.cv_row <- function(id, context) {
  config <- context$config
  indices <- as.integer(context$grid[id, ])
  loss <- rep(Inf, length(config$alpha))
  diagnostics <- vector("list", length(config$alpha))
  init <- NULL
  for (a in rev(seq_along(config$alpha))) {
    fit <- tryCatch(.simulation_fit(context$split$train, config,
      config$lambda[indices], config$alpha[a], init = init, steps = context$steps),
      error = identity)
    if (inherits(fit, "error")) {
      diagnostics[[a]] <- list(error = conditionMessage(fit))
      init <- NULL
    } else {
      score <- .validation_loss(fit, context$split$test, config$family)
      if (is.finite(score)) {
        loss[a] <- score
        diagnostics[[a]] <- fit$diagnostics
        init <- fit
      } else {
        diagnostics[[a]] <- list(error = "Non-finite validation loss")
        init <- NULL
      }
    }
  }
  list(id = id, loss = loss, diagnostics = diagnostics)
}

.run_one <- function(config, position, cluster = NULL) {
  data <- simulate_jlis(config$setting, position, config$profile)
  data$config <- config
  grid <- expand.grid(i1 = seq_along(config$lambda), i2 = seq_along(config$lambda),
                      i3 = seq_along(config$lambda))
  cv <- matrix(0, nrow(grid), length(config$alpha))
  diagnostic_rows <- list()
  failures <- list()
  for (fold in seq_len(config$folds)) {
    split <- .fold_data(data, fold)
    steps <- if (config$family == "binomial")
      .alpha_logistic$steps(split$train$X_list, split$train$Z_list) else NULL
    context <- list(config = config, grid = grid, split = split, steps = steps)
    if (is.null(cluster)) rows <- lapply(seq_len(nrow(grid)), .cv_row, context = context)
    else {
      export <- list2env(list(.jlis_cv_context = context), parent = emptyenv())
      parallel::clusterExport(cluster, ".jlis_cv_context", envir = export)
      rows <- parallel::parLapply(cluster, seq_len(nrow(grid)), function(id) {
        utils::getFromNamespace(".cv_row", "jlis")(id, get(".jlis_cv_context", envir = .GlobalEnv))
      })
    }
    for (row in rows) cv[row$id, ] <- cv[row$id, ] + row$loss
    for (a in seq_along(config$alpha)) {
      d <- lapply(rows, function(row) row$diagnostics[[a]])
      ok <- !vapply(d, function(x) !is.null(x$error), logical(1))
      for (bad in which(!ok)) failures[[length(failures) + 1L]] <- data.frame(
        fold = fold, grid_id = rows[[bad]]$id, alpha_id = a, message = d[[bad]]$error)
      good <- d[ok]
      avg <- function(n) if (length(good)) mean(vapply(good, function(x) x[[n]], numeric(1))) else NA_real_
      violation <- if (length(good)) max(vapply(good, function(x)
        x$max_abs_L - x$alpha_bound, numeric(1))) else NA_real_
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
        fold = fold, alpha_id = a, alpha = config$alpha[a], n_fits = length(d),
        failures = sum(!ok), mean_iter = avg("iter"), mean_primal = avg("primal"),
        mean_rel_change = avg("rel_change"),
        mean_dykstra_active_count = avg("dykstra_active_count"),
        max_abs_L_minus_bound = violation)
    }
    message(config$setting, " replicate ", position, ": fold ", fold, "/", config$folds)
  }
  if (!any(is.finite(cv))) stop("Every tuning combination failed; inspect the data and controls", call. = FALSE)
  best <- arrayInd(which.min(cv), dim(cv))[1, ]
  indices <- as.integer(grid[best[1], ])
  selected <- c(config$lambda[indices], config$alpha[best[2]])
  names(selected) <- c("lambda1", "lambda2", "lambda3", "alpha")
  steps <- if (config$family == "binomial") .alpha_logistic$steps(data$X_list, data$Z_list) else NULL
  fit <- .simulation_fit(data, config, selected[1:3], selected[4], final = TRUE, steps = steps)
  fit$diagnostics$converged <- max(fit$diagnostics$primal, fit$diagnostics$rel_change) < config$control$tol
  list(replicate = position, turn = config$turns[position], seed = config$seeds[position],
       selected_indices = setNames(c(indices, best[2]), c("lambda1_idx", "lambda2_idx", "lambda3_idx", "alpha_idx")),
       selected_values = selected,
       B_error = vapply(seq_len(config$G), function(g) sqrt(sum((fit$B[[g]] - data$B_list[[g]])^2)), numeric(1)),
       gamma_error = vapply(seq_len(config$G), function(g) sqrt(sum((fit$gamma[[g]] - data$gamma_list[[g]])^2)), numeric(1)),
       rank = .simulation_engine(config)$rank(fit$L, config$p, config$q),
       final_diagnostics = fit$diagnostics, cv_min = min(cv),
       alpha_cv_diagnostics = do.call(rbind, diagnostic_rows),
       failures = if (length(failures)) do.call(rbind, failures) else data.frame(),
       runtime = list(R = as.character(getRversion()), pracma = as.character(utils::packageVersion("pracma"))))
}

.load_part <- function(file) {
  e <- new.env(parent = emptyenv())
  load(file, envir = e)
  if (!all(c("rep_result", "config") %in% ls(e))) stop("Invalid result file: ", basename(file), call. = FALSE)
  list(result = e$rep_result, config = e$config)
}

.check_part <- function(part, config, position) {
  r <- part$result
  if (!identical(part$config, config) || !identical(as.integer(r$replicate), as.integer(position)) ||
      r$turn != config$turns[position] || r$seed != config$seeds[position])
    stop("Result configuration or repetition mismatch", call. = FALSE)
  if (length(r$B_error) != config$G || length(r$gamma_error) != config$G ||
      any(!is.finite(c(r$B_error, r$gamma_error, r$rank, r$selected_values, r$cv_min))))
    stop("Result contains incomplete coefficient errors", call. = FALSE)
  d <- r$final_diagnostics
  if (!is.finite(d$max_abs_L) || d$max_abs_L > d$alpha_bound + 1e-8)
    stop("Result violates the entrywise constraint", call. = FALSE)
  invisible(TRUE)
}

run_simulation <- function(setting = "A1", replicates = 1:50,
                           outdir = file.path("results", setting), cores = 1L,
                           profile = c("publication", "smoke"), resume = TRUE) {
  config <- simulation_settings(setting, match.arg(profile))
  if (!is.numeric(replicates) || !length(replicates) || anyNA(replicates) ||
      any(!replicates %in% 1:50) || anyDuplicated(replicates))
    stop("replicates must contain distinct integers between 1 and 50", call. = FALSE)
  replicates <- as.integer(replicates)
  cores <- .positive_integer(cores, "cores")
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) stop("resume must be TRUE or FALSE", call. = FALSE)
  part_dir <- file.path(outdir, "parts")
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  cluster <- NULL
  if (cores > 1L) {
    cluster <- parallel::makePSOCKcluster(min(cores, length(config$lambda)^3))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(paths) { .libPaths(paths); loadNamespace("jlis"); NULL }, .libPaths())
  }
  paths <- character(length(replicates))
  for (i in seq_along(replicates)) {
    position <- replicates[i]
    path <- file.path(part_dir, sprintf("rep_%03d.RData", position))
    paths[i] <- path
    if (file.exists(path) && resume) {
      .check_part(.load_part(path), config, position)
      message("Validated existing result: ", basename(path))
      next
    }
    rep_result <- .run_one(config, position, cluster)
    .check_part(list(result = rep_result, config = config), config, position)
    temporary <- tempfile("part-", tmpdir = part_dir, fileext = ".RData")
    save(rep_result, config, file = temporary, version = 3)
    if (!file.rename(temporary, path)) {
      unlink(temporary)
      stop("Could not finalize result: ", basename(path), call. = FALSE)
    }
  }
  invisible(paths)
}
