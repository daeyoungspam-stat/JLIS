.alpha_linear_b <- local({
safe_svd <- function(A) {
    if (!all(is.finite(A))) {
        stop("non-finite matrix entries before SVD")
    }
    tryCatch(svd(A), error = function(e) {
        eig <- eigen(crossprod(A), symmetric = TRUE)
        d <- sqrt(pmax(eig$values, 0))
        v <- eig$vectors
        u <- matrix(0, nrow = nrow(A), ncol = ncol(A))
        keep <- d > .Machine$double.eps * max(dim(A)) * max(max(d), 1)
        if (any(keep)) {
            u[, keep] <- sweep(A %*% v[, keep, drop = FALSE], 2, d[keep], "/")
        }
        list(d = d, u = u, v = v)
    })
}
soft_threshold <- function(x, thresh) {
    sign(x) * pmax(abs(x) - thresh, 0)
}
svt <- function(A, tau) {
    sv <- safe_svd(A)
    d <- pmax(sv$d - tau, 0)
    keep <- which(d > 0)
    if (length(keep) == 0) {
        return(matrix(0, nrow = nrow(A), ncol = ncol(A)))
    }
    sv$u[, keep, drop = FALSE] %*% diag(d[keep], nrow = length(keep), ncol = length(keep)) %*% t(sv$v[, 
        keep, drop = FALSE])
}
project_inf_ball <- function(A, bound) {
    if (is.infinite(bound)) {
        return(A)
    }
    pmin(pmax(A, -bound), bound)
}
prox_nuclear_box <- function(M, tau, alpha, p, q, dykstra_max_iter = 200, dykstra_tol = 1e-05) {
    bound <- if (is.infinite(alpha)) 
        Inf
    else alpha/sqrt(p * q)
    L_svt <- svt(M, tau)
    if (is.infinite(bound) || max(abs(L_svt)) <= bound + 1e-10) {
        return(list(L = L_svt, active = FALSE, iter = 0, bound = bound))
    }
    R <- M
    P <- matrix(0, nrow = p, ncol = q)
    Q <- matrix(0, nrow = p, ncol = q)
    for (iter in seq_len(dykstra_max_iter)) {
        Y <- svt(R + P, tau)
        P <- R + P - Y
        R_next <- project_inf_ball(Y + Q, bound)
        Q <- Y + Q - R_next
        rel_change <- sqrt(sum((R_next - R)^2))/(1 + sqrt(sum(R^2)))
        R <- R_next
        if (rel_change < dykstra_tol && max(abs(R)) <= bound + 1e-08) {
            break
        }
    }
    list(L = R, active = TRUE, iter = iter, bound = bound)
}
coordinate_descent_lasso <- function(y, Z, lambda, beta_init = NULL, tol = 1e-06, max_iter = 1000, eps = 1e-08) {
    n <- length(y)
    p <- ncol(Z)
    if (is.null(beta_init)) {
        beta <- rep(0, p)
    }
    else {
        beta <- as.numeric(beta_init)
    }
    z_norm2 <- colSums(Z^2)/n
    z_norm2[z_norm2 < eps] <- eps
    resid <- as.numeric(y - Z %*% beta)
    for (iter in seq_len(max_iter)) {
        beta_old <- beta
        for (j in seq_len(p)) {
            resid <- resid + Z[, j] * beta[j]
            zj_tr <- sum(Z[, j] * resid)/n
            beta[j] <- soft_threshold(zj_tr, lambda/2)/z_norm2[j]
            resid <- resid - Z[, j] * beta[j]
        }
        rel_change <- sqrt(sum((beta - beta_old)^2))/(sqrt(sum(beta_old^2)) + eps)
        if (rel_change < tol) {
            break
        }
    }
    beta
}
linear_admm_alpha <- function(Y_list, X_list, Z_list, px, qx, lambda1, lambda2, lambda3, alpha = Inf, 
    omega1, omega2, rho = 1, max_iter = 500, tol = 1e-04, gamma_max_iter = 500, gamma_tol = 1e-06, dykstra_max_iter = 200, 
    dykstra_tol = 1e-05, init = NULL, return_diagnostics = TRUE) {
    G <- length(Y_list)
    if (missing(omega1)) {
        omega1 <- rep(1, G)
    }
    if (missing(omega2)) {
        omega2 <- rep(1, G)
    }
    p_times_q <- ncol(X_list[[1]])
    p <- px
    q <- qx
    r <- ncol(Z_list[[1]])
    alpha_bound <- if (is.infinite(alpha)) 
        Inf
    else alpha/sqrt(p * q)
    precomp <- vector("list", G)
    for (g in seq_len(G)) {
        Xg <- X_list[[g]]
        Zg <- Z_list[[g]]
        Yg <- Y_list[[g]]
        n_g <- nrow(Xg)
        A_mat <- crossprod(Xg)/n_g + (rho/2) * diag(p_times_q)
        precomp[[g]] <- list(n = n_g, chol_A = chol(A_mat), XtY = crossprod(Xg, Yg), XtZ = crossprod(Xg, 
            Zg))
    }
    if (is.null(init)) {
        B <- replicate(G, matrix(0, nrow = p_times_q, ncol = 1), simplify = FALSE)
        S <- replicate(G, matrix(0, nrow = p_times_q, ncol = 1), simplify = FALSE)
        gamma <- replicate(G, matrix(0, nrow = r, ncol = 1), simplify = FALSE)
        U <- replicate(G, matrix(0, nrow = p_times_q, ncol = 1), simplify = FALSE)
        L <- matrix(0, nrow = p_times_q, ncol = 1)
    }
    else {
        B <- init$B
        S <- init$S
        gamma <- init$gamma
        U <- init$U
        L <- matrix(init$L, nrow = p_times_q, ncol = 1)
        if (!is.infinite(alpha_bound)) {
            L <- as.vector(project_inf_ball(matrix(L, nrow = p, ncol = q), alpha_bound))
            L <- matrix(L, nrow = p_times_q, ncol = 1)
        }
    }
    dykstra_active_count <- 0
    dykstra_iter_total <- 0
    last_primal <- Inf
    last_rel_change <- Inf
    for (iter in seq_len(max_iter)) {
        L_prev <- L
        B_prev <- B
        gamma_prev <- gamma
        B <- lapply(seq_len(G), function(g) {
            gamg <- gamma[[g]]
            rhs <- precomp[[g]]$XtY/precomp[[g]]$n - precomp[[g]]$XtZ %*% gamg/precomp[[g]]$n + (rho/2) * 
                (L + S[[g]] - U[[g]])
            backsolve(precomp[[g]]$chol_A, forwardsolve(t(precomp[[g]]$chol_A), rhs))
        })
        gamma <- lapply(seq_len(G), function(g) {
            v_vec <- Y_list[[g]] - X_list[[g]] %*% B[[g]]
            coordinate_descent_lasso(y = v_vec, Z = Z_list[[g]], lambda = lambda3 * omega2[g], beta_init = gamma[[g]], 
                tol = gamma_tol, max_iter = gamma_max_iter)
        })
        gamma <- lapply(gamma, function(x) matrix(x, nrow = r, ncol = 1))
        M <- Reduce("+", lapply(seq_len(G), function(g) {
            matrix(B[[g]] - S[[g]] + U[[g]], nrow = p, ncol = q)
        }))/G
        prox_L <- prox_nuclear_box(M = M, tau = lambda2/(rho * G), alpha = alpha, p = p, q = q, dykstra_max_iter = dykstra_max_iter, 
            dykstra_tol = dykstra_tol)
        L <- matrix(as.vector(prox_L$L), nrow = p_times_q, ncol = 1)
        dykstra_active_count <- dykstra_active_count + as.integer(prox_L$active)
        dykstra_iter_total <- dykstra_iter_total + prox_L$iter
        S <- lapply(seq_len(G), function(g) {
            soft_threshold(B[[g]] - L + U[[g]], (lambda1 * omega1[g])/rho)
        })
        U <- lapply(seq_len(G), function(g) {
            U[[g]] + (B[[g]] - L - S[[g]])
        })
        diff_B <- max(sapply(seq_len(G), function(g) {
            norm(B[[g]] - B_prev[[g]], type = "F")/(norm(B_prev[[g]], type = "F") + 1e-08)
        }))
        diff_gamma <- max(sapply(seq_len(G), function(g) {
            norm(gamma[[g]] - gamma_prev[[g]], type = "F")/(norm(gamma_prev[[g]], type = "F") + 1e-08)
        }))
        diff_L <- norm(L - L_prev, type = "F")/(norm(L_prev, type = "F") + 1e-08)
        last_primal <- max(sapply(seq_len(G), function(g) {
            norm(B[[g]] - L - S[[g]], type = "F")/(1 + norm(B[[g]], type = "F"))
        }))
        last_rel_change <- max(diff_B, diff_gamma, diff_L)
        if (max(last_rel_change, last_primal) < tol) {
            break
        }
    }
    diagnostics <- list(iter = iter, rel_change = last_rel_change, primal = last_primal, alpha = alpha, 
        alpha_bound = alpha_bound, max_abs_L = max(abs(L)), dykstra_active_count = dykstra_active_count, 
        dykstra_iter_total = dykstra_iter_total)
    out <- list(B = B, L = L, S = S, gamma = gamma, U = U)
    if (return_diagnostics) {
        out$diagnostics <- diagnostics
    }
    out
}
estimate_rank <- function(L, p, q, tol = 0.001) {
    sum(safe_svd(matrix(L, nrow = p, ncol = q))$d > tol)
}
list(fit = linear_admm_alpha, rank = estimate_rank, prox = prox_nuclear_box)
})
