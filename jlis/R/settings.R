.positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x))
    stop(name, " must be a positive integer", call. = FALSE)
  as.integer(x)
}

.with_seed <- function(seed, fun) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  fun()
}

simulation_settings <- function(setting = "A1", profile = c("publication", "smoke")) {
  setting <- match.arg(toupper(setting), c("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2"))
  profile <- match.arg(profile)
  ids <- c("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2")
  i <- match(setting, ids)
  family <- if (i <= 4L) "gaussian" else "binomial"
  p <- c(20L, 20L, 40L, 40L, 15L, 15L, 30L, 30L)[i]
  r <- c(25L, 25L, 40L, 40L, 25L, 25L, 40L, 40L)[i]
  rank <- c(2L, 3L, 4L, 5L, 2L, 3L, 3L, 4L)[i]
  n_fold <- c(40L, 40L, 200L, 200L, 20L, 20L, 60L, 60L)[i]
  turns <- if (setting == "D1") c(1:9, 11:51) else if (setting == "D2") 51:100 else 1:50
  seed_index <- if (setting == "A1") 51:100 else turns
  master <- c(726L, 726L, 55455L, 123L, 78256L, 55414L, 4321L, 4141L)[i]
  draws <- c(150L, 150L, 50L, 50L, 50L, 50L, 200L, 200L)[i]
  upper <- if (i >= 7L) 1000L else 100000L
  seeds <- .with_seed(master, function() sample(seq_len(upper), draws, replace = TRUE))
  lambda <- 10^seq(-3, if (setting %in% c("B1", "B2")) 0.5 else 3, by = 0.5)
  alpha <- exp(seq(log(if (family == "gaussian") 2 else 1), log(p), length.out = 10))
  control <- list(cv_max_iter = if (family == "gaussian") 500L else 250L,
                  final_max_iter = if (family == "gaussian") 800L else 500L,
                  tol = 1e-4, pgd_max_inner = 200L, pgd_tol = 1e-4,
                  dykstra_max_iter = if (family == "gaussian") 200L else 150L,
                  dykstra_tol = 1e-5)
  if (profile == "smoke") {
    p <- 6L; r <- 6L; rank <- 2L; n_fold <- 8L
    lambda <- 0.1
    alpha <- c(1, 6)
    control$cv_max_iter <- 12L
    control$final_max_iter <- 24L
    control$pgd_max_inner <- 12L
    control$dykstra_max_iter <- 20L
  }
  list(setting = setting, profile = profile, family = family,
       G = 4L, p = p, q = p, r = r, rank = rank,
       folds = 5L, n_per_fold = rep(n_fold, 4L), n = rep(n_fold * 5L, 4L),
       positions = seq_len(50L), turns = as.integer(turns),
       seeds = as.integer(seeds[seed_index]),
       lambda = lambda, alpha = alpha, control = control,
       alpha_svd = if (setting %in% c("A1", "A2")) "standard" else "safe")
}

simulate_jlis <- function(setting = "A1", replicate = 1L,
                          profile = c("publication", "smoke")) {
  config <- simulation_settings(setting, match.arg(profile))
  replicate <- .positive_integer(replicate, "replicate")
  if (replicate > 50L) stop("replicate must be at most 50", call. = FALSE)
  out <- .with_seed(config$seeds[replicate], function() {
    common <- list(G = config$G, n = config$n_per_fold, p = config$p,
                   q = config$q, r = config$r, rank_L = config$rank,
                   sparsegam = 0.7, fold = config$folds)
    if (config$family == "gaussian") {
      is_a <- config$setting %in% c("A1", "A2")
      singular <- switch(config$setting,
        A1 = c(0.5, 0.1, 0.1), A2 = c(0.5, 0.5, 0.1),
        B1 = c(0.5, 0.5, 0.1, 0.1), B2 = c(0.5, 0.5, 0.5, 0.1, 0.1))
      generator <- if (is_a) .dgp_linear_a$generate else .dgp_linear_b$generate
      extra <- if (is_a) list(scale = 1) else list()
      do.call(generator, c(common, extra, list(
        signal = if (is_a) 0.7 else 0.2, gamma = if (is_a) 0.75 else 0.2,
        singular_vals = singular, umingam = 0.5, umaxgam = 1, sigma = 1)))
    } else {
      is_c <- config$setting %in% c("C1", "C2")
      do.call(.dgp_logistic$generate, c(common, list(
        uminL = if (is_c) 1 else 0.75, umaxL = if (is_c) 2 else 1.5,
        uminS = if (is_c) 0.01 else 0.0001, umaxS = if (is_c) 0.1 else 0.01,
        umingam = if (is_c) 0.5 else 0.3, umaxgam = if (is_c) 1 else 0.6,
        Sblockrow = 3, Sblockcol = 3, Sblocknum = 1)))
    }
  })
  out$config <- config
  out$replicate <- replicate
  out$turn <- config$turns[replicate]
  out
}
