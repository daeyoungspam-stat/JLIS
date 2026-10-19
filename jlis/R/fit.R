.validate_design <- function(X_list, Z_list, p, q, Y_list = NULL, family = "gaussian") {
  if (!is.list(X_list) || !length(X_list) || !is.list(Z_list) ||
      length(X_list) != length(Z_list)) stop("X_list and Z_list require matching nonempty lists", call. = FALSE)
  if (!is.null(Y_list) && (!is.list(Y_list) || length(Y_list) != length(X_list)))
    stop("Y_list must have one element per group", call. = FALSE)
  r <- ncol(Z_list[[1]])
  for (g in seq_along(X_list)) {
    X <- X_list[[g]]; Z <- Z_list[[g]]
    if (!is.matrix(X) || !is.numeric(X) || nrow(X) < 1L || ncol(X) != p * q ||
        any(!is.finite(X))) stop("Invalid X_list element in group ", g, call. = FALSE)
    if (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != nrow(X) ||
        ncol(Z) != r || any(!is.finite(Z))) stop("Invalid Z_list element in group ", g, call. = FALSE)
    if (!is.null(Y_list)) {
      y <- Y_list[[g]]
      if (!is.numeric(y) || length(y) != nrow(X) || any(!is.finite(y)) ||
          (!is.null(dim(y)) && (!is.matrix(y) || ncol(y) != 1L)))
        stop("Y_list elements require one numeric response per row", call. = FALSE)
      if (family == "binomial" && any(!y %in% c(0, 1)))
        stop("Binomial responses must be 0 or 1", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.fit_control <- function(control, family) {
  defaults <- list(rho = 1, max_iter = if (family == "gaussian") 800L else 500L,
    tol = 1e-4, dykstra_max_iter = if (family == "gaussian") 200L else 150L,
    dykstra_tol = 1e-5, svd = "safe")
  if (family == "gaussian") defaults <- c(defaults, list(gamma_max_iter = 500L, gamma_tol = 1e-6))
  else defaults <- c(defaults, list(rankthres = 1e-3, pgd_max_inner = 200L, pgd_tol = 1e-4))
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown control: ", paste(unknown, collapse = ", "), call. = FALSE)
  if (length(control) && (is.null(names(control)) || any(!nzchar(names(control)))))
    stop("control entries must be named", call. = FALSE)
  defaults[names(control)] <- control
  defaults$svd <- match.arg(defaults$svd, c("safe", "standard"))
  for (n in setdiff(names(defaults), "svd")) {
    x <- defaults[[n]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 0 || (n != "rankthres" && x == 0)) stop("Invalid control: ", n, call. = FALSE)
    if (grepl("iter|inner", n)) defaults[[n]] <- .positive_integer(x, n)
  }
  defaults
}

fit_jlis <- function(Y_list, X_list, Z_list, px, qx,
                     lambda1, lambda2, lambda3, alpha = Inf,
                     family = c("gaussian", "binomial"),
                     omega1 = rep(1, length(Y_list)), omega2 = rep(1, length(Y_list)),
                     control = list(), init = NULL) {
  family <- match.arg(family)
  px <- .positive_integer(px, "px"); qx <- .positive_integer(qx, "qx")
  .validate_design(X_list, Z_list, px, qx, Y_list, family)
  penalty <- list(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3)
  if (!all(vapply(penalty, function(x) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= 0, logical(1)))) stop("Penalties must be finite nonnegative scalars", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0)
    stop("alpha must be positive or Inf", call. = FALSE)
  for (w in list(omega1, omega2)) if (!is.numeric(w) || length(w) != length(Y_list) ||
    any(!is.finite(w)) || any(w < 0)) stop("Weights require one nonnegative value per group", call. = FALSE)
  control <- .fit_control(control, family)
  engine <- if (family == "binomial") .alpha_logistic else
    if (control$svd == "standard") .alpha_linear_a else .alpha_linear_b
  control$svd <- NULL
  if (!is.null(init)) {
    required <- c("B", "L", "S", "gamma", "U")
    if (!all(required %in% names(init))) stop("init requires B, L, S, gamma, and U", call. = FALSE)
    for (n in c("B", "S", "U", "gamma")) {
      d <- if (n == "gamma") ncol(Z_list[[1]]) else px * qx
      if (length(init[[n]]) != length(Y_list) ||
          any(!vapply(init[[n]], function(x) is.numeric(x) && length(x) == d &&
                      all(is.finite(x)), logical(1)))) stop("Invalid init component: ", n, call. = FALSE)
      init[[n]] <- lapply(init[[n]], matrix, nrow = d, ncol = 1)
    }
    if (length(init$L) != px * qx || any(!is.finite(init$L))) stop("Invalid init L", call. = FALSE)
  }
  fit <- do.call(engine$fit, c(list(Y_list = lapply(Y_list, matrix, ncol = 1),
    X_list = X_list, Z_list = Z_list, px = px, qx = qx,
    alpha = alpha, omega1 = omega1, omega2 = omega2, init = init), penalty, control))
  fit$family <- family; fit$px <- px; fit$qx <- qx; fit$r <- ncol(Z_list[[1]])
  fit$penalties <- unlist(penalty)
  fit$diagnostics$converged <- max(fit$diagnostics$primal, fit$diagnostics$rel_change) < control$tol
  class(fit) <- c("jlis_fit", "list")
  fit
}

predict.jlis_fit <- function(object, X_list, Z_list, type = c("response", "link"), ...) {
  type <- match.arg(type)
  .validate_design(X_list, Z_list, object$px, object$qx)
  if (length(X_list) != length(object$B) || ncol(Z_list[[1]]) != object$r)
    stop("Prediction groups and covariate dimensions must match the fit", call. = FALSE)
  lapply(seq_along(X_list), function(g) {
    eta <- as.numeric(X_list[[g]] %*% object$B[[g]] + Z_list[[g]] %*% object$gamma[[g]])
    if (type == "response" && object$family == "binomial") stats::plogis(eta) else eta
  })
}

coef.jlis_fit <- function(object, ...) {
  reshape <- function(x) matrix(x, object$px, object$qx)
  list(L = reshape(object$L), B = lapply(object$B, reshape),
       S = lapply(object$S, reshape), gamma = lapply(object$gamma, as.numeric))
}

print.jlis_fit <- function(x, ...) {
  cat("JLIS ", x$family, " fit: ", length(x$B), " groups, ", x$px, " x ", x$qx, " matrices\n", sep = "")
  cat("Iterations:", x$diagnostics$iter, "  Converged:", x$diagnostics$converged,
      "  Entrywise bound:", x$diagnostics$alpha_bound, "\n")
  invisible(x)
}
