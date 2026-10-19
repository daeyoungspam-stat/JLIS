.alpha_logistic <- local({
soft_threshold <- function(x, thresh) {
    sign(x) * pmax(abs(x) - thresh, 0)
}
softplus_stable <- function(x) {
    out <- numeric(length(x))
    pos <- x > 0
    out[pos] <- x[pos] + log1p(exp(-x[pos]))
    out[!pos] <- log1p(exp(x[!pos]))
    out
}
sigmoid_stable <- function(x) {
    out <- numeric(length(x))
    pos <- x >= 0
    out[pos] <- 1/(1 + exp(-x[pos]))
    ex <- exp(x[!pos])
    out[!pos] <- ex/(1 + ex)
    out
}
g_value <- function(eta, B, M, y, rho) {
    n <- length(y)
    loss_logit <- mean(softplus_stable(eta) - y * eta)
    penalty <- 0.5 * rho * sum((B - M)^2)
    loss_logit + penalty
}
pgd_logit <- function(X, Z, y, L_vec, S_vec, U_vec, lambda3, omega2g, rho, B_init, gamma_init, t_init = NULL, 
    eta_bt = 2, max_inner = 1000, tol_inner = 1e-04) {
    n <- length(y)
    M <- L_vec + S_vec - U_vec
    B <- B_init
    gam <- gamma_init
    if (is.null(t_init)) {
        normX.sq <- sum(X^2)
        normZ.sq <- sum(Z^2)
        Lip.const <- rho + 0.25 * (normX.sq + normZ.sq)/n
        t <- 1/(2 * Lip.const)
    }
    else t <- t_init
    eta <- as.vector(X %*% B + Z %*% gam)
    p <- sigmoid_stable(eta)
    py <- (p - y)/n
    gB <- crossprod(X, py) + rho * (B - M)
    ggam <- crossprod(Z, py)
    gx <- g_value(eta, B, M, y, rho)
    for (it in 1:max_inner) {
        Bold <- B
        gamold <- gam
        repeat {
            Bc <- B - t * gB
            gamc0 <- gam - t * ggam
            gamc <- soft_threshold(gamc0, t * lambda3 * omega2g)
            dB <- Bc - B
            dgam <- gamc - gam
            eta_new <- eta + as.vector(X %*% dB) + as.vector(Z %*% dgam)
            gy <- g_value(eta_new, Bc, M, y, rho)
            quad <- as.numeric(sum(gB * dB) + sum(ggam * dgam)) + (1/(2 * t)) * (sum(dB^2) + sum(dgam^2))
            if (gy <= gx + quad) {
                B <- Bc
                gam <- gamc
                eta <- eta_new
                gx <- gy
                break
            }
            else {
                t <- t/eta_bt
                if (t < 1e-16) 
                  break
            }
        }
        p <- sigmoid_stable(eta)
        py <- (p - y)/n
        gB_new <- crossprod(X, py) + rho * (B - M)
        ggam_new <- crossprod(Z, py)
        gx_new <- g_value(eta, B, M, y, rho)
        gB <- gB_new
        ggam <- ggam_new
        gx <- gx_new
        dB <- B - Bold
        dgam <- gam - gamold
        step_norm <- sqrt(sum(dB^2) + sum(dgam^2))
        base_norm <- sqrt(sum(Bold^2) + sum(gamold^2)) + 1e-12
        if (step_norm/base_norm < tol_inner) 
            break
    }
    list(B = B, gamma = gam, step = t, inner_iter = it)
}
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
pgd_initial_steps <- function(X_list, Z_list, rho = 1) {
    if (length(X_list) != length(Z_list)) {
        stop("X_list and Z_list must have the same number of groups")
    }
    vapply(seq_along(X_list), function(g) {
        n <- nrow(X_list[[g]])
        lip_const <- rho + 0.25 * (sum(X_list[[g]]^2) + sum(Z_list[[g]]^2))/n
        1/(2 * lip_const)
    }, numeric(1))
}
logistic_admm_alpha <- function(Y_list, X_list, Z_list, px, qx, lambda1, lambda2, lambda3, alpha = Inf, 
    omega1, omega2, rho = 1, max_iter = 300, tol = 1e-04, rankthres = 0.001, pgd_max_inner = 300, pgd_tol = 1e-04, 
    pgd_t_init = NULL, dykstra_max_iter = 200, dykstra_tol = 1e-05, init = NULL, return_diagnostics = TRUE) {
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
    if (is.null(pgd_t_init)) {
        pgd_t_init <- pgd_initial_steps(X_list, Z_list, rho)
    }
    if (length(pgd_t_init) != G || any(!is.finite(pgd_t_init)) || any(pgd_t_init <= 0)) {
        stop("pgd_t_init must contain one finite positive step per group")
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
    last_inner_iter <- 0
    for (iter in seq_len(max_iter)) {
        L_prev <- L
        B_prev <- B
        gamma_prev <- gamma
        BGam_updates <- lapply(seq_len(G), function(g) {
            pgd_logit(X = X_list[[g]], Z = Z_list[[g]], y = Y_list[[g]], L_vec = L, S_vec = S[[g]], U_vec = U[[g]], 
                lambda3 = lambda3, omega2g = omega2[g], rho = rho, B_init = B[[g]], gamma_init = gamma[[g]], 
                t_init = pgd_t_init[g], eta_bt = 2, max_inner = pgd_max_inner, tol_inner = pgd_tol)
        })
        B <- lapply(BGam_updates, `[[`, "B")
        gamma <- lapply(BGam_updates, `[[`, "gamma")
        last_inner_iter <- mean(vapply(BGam_updates, `[[`, numeric(1), "inner_iter"))
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
    if (rankthres > 0) {
        post <- prox_nuclear_box(M = matrix(L, nrow = p, ncol = q), tau = rankthres, alpha = alpha, p = p, 
            q = q, dykstra_max_iter = dykstra_max_iter, dykstra_tol = dykstra_tol)
        L <- matrix(as.vector(post$L), nrow = p_times_q, ncol = 1)
    }
    diagnostics <- list(iter = iter, rel_change = last_rel_change, primal = last_primal, mean_inner_iter = last_inner_iter, 
        alpha = alpha, alpha_bound = alpha_bound, max_abs_L = max(abs(L)), dykstra_active_count = dykstra_active_count, 
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
list(fit = logistic_admm_alpha, rank = estimate_rank, steps = pgd_initial_steps, prob = sigmoid_stable)
})
