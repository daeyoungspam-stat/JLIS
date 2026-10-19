.legacy_linear <- local({
soft_threshold <- function(x, thresh) {
    sign(x) * pmax(abs(x) - thresh, 0)
}
coordinate_descent_lasso <- function(y, Z, lambda, tol = 1e-06, max_iter = 1000, eps = 1e-08) {
    n.y <- length(y)
    p.z <- ncol(Z)
    fits <- rep(0, p.z)
    Z_norm2 <- colSums(Z^2)/n.y
    r <- y
    for (iter in 1:max_iter) {
        fits_old <- fits
        for (j in 1:p.z) {
            r <- r + Z[, j] * fits[j]
            zjTr <- sum(Z[, j] * r)/n.y
            fits[j] <- soft_threshold(zjTr, lambda/2)/Z_norm2[j]
            r <- r - Z[, j] * fits[j]
        }
        rel_change <- sqrt(sum((fits - fits_old)^2))/(sqrt(sum(fits_old^2)) + eps)
        if (rel_change < tol) {
            break
        }
    }
    return(fits)
}
linear_admm <- function(Y_list, X_list, Z_list, px, qx, lambda1, lambda2, lambda3, omega1, omega2, rho = 1, 
    max_iter = 1000, tol = 1e-04, rankthres = 0.001) {
    G <- length(Y_list)
    if (missing(omega1)) {
        omega1 <- rep(1, G)
    }
    if (missing(omega2)) {
        omega2 <- rep(1, G)
    }
    p.times.q <- ncol(X_list[[1]])
    p <- px
    q <- qx
    r <- ncol(Z_list[[1]])
    precomp <- vector("list", G)
    for (g in 1:G) {
        Xg <- X_list[[g]]
        n_g <- nrow(Xg)
        XtX <- crossprod(Xg)
        A_mat <- (1/n_g) * XtX + (rho/2) * diag(p.times.q)
        precomp[[g]] <- list(n = n_g, chol_A = chol(A_mat))
    }
    B <- vector("list", G)
    S <- vector("list", G)
    gamma <- vector("list", G)
    U <- vector("list", G)
    for (g in 1:G) {
        B[[g]] <- matrix(0, nrow = p.times.q, ncol = 1)
        S[[g]] <- matrix(0, nrow = p.times.q, ncol = 1)
        gamma[[g]] <- matrix(0, nrow = r, ncol = 1)
        U[[g]] <- matrix(0, nrow = p.times.q, ncol = 1)
    }
    L <- matrix(0, nrow = p.times.q, ncol = 1)
    for (iter in 1:max_iter) {
        L_prev <- L
        B_prev <- B
        gamma_prev <- gamma
        B <- lapply(1:G, function(g) {
            Xg <- X_list[[g]]
            Yg <- Y_list[[g]]
            Zg <- Z_list[[g]]
            gamg <- gamma[[g]]
            Ug <- U[[g]]
            Sg <- S[[g]]
            n_g <- precomp[[g]]$n
            Y_mod <- Yg - Zg %*% gamg
            rhs <- (1/n_g) * crossprod(Xg, Y_mod) + (rho/2) * (L + Sg - Ug)
            chol_A <- precomp[[g]]$chol_A
            backsolve(chol_A, forwardsolve(t(chol_A), rhs))
        })
        gamma <- lapply(1:G, function(g) {
            Xg <- X_list[[g]]
            Yg <- Y_list[[g]]
            Zg <- Z_list[[g]]
            Bg <- B[[g]]
            v.vec <- Yg - Xg %*% Bg
            coordinate_descent_lasso(y = v.vec, Z = Zg, lambda = lambda3 * omega2[g], tol = 1e-06, max_iter = 1000, 
                eps = 1e-08)
        })
        gamma <- lapply(gamma, function(x) matrix(x, nrow = r, ncol = 1))
        M <- Reduce("+", lapply(1:G, function(g) {
            matrix(B[[g]] - S[[g]] + U[[g]], nrow = p, ncol = q)
        }))/G
        svdM <- svd(M)
        thresh_val <- lambda2/(rho * G)
        D_thresholded <- pmax(svdM$d - thresh_val, 0)
        L_mat <- svdM$u %*% diag(D_thresholded, nrow = length(D_thresholded)) %*% t(svdM$v)
        L <- as.vector(L_mat)
        S <- lapply(1:G, function(g) {
            soft_threshold(B[[g]] - L + U[[g]], (lambda1 * omega1[g])/rho)
        })
        U <- lapply(1:G, function(g) {
            U[[g]] + (B[[g]] - L - S[[g]])
        })
        diff_B <- max(sapply(1:G, function(g) norm(B[[g]] - B_prev[[g]], type = "F")/(norm(B_prev[[g]], 
            type = "F") + 1e-08)))
        diff_gamma <- max(sapply(1:G, function(g) norm(gamma[[g]] - gamma_prev[[g]], type = "F")/(norm(gamma_prev[[g]], 
            type = "F") + 1e-08)))
        if (max(diff_B, diff_gamma) < tol) 
            break
    }
    L <- matrix(L, nrow = p, ncol = q)
    svdL <- svd(L)
    LD_thresholded <- pmax(svdL$d - rankthres, 0)
    L_mat <- svdL$u %*% diag(LD_thresholded, nrow = length(LD_thresholded)) %*% t(svdL$v)
    L <- as.vector(L_mat)
    return(list(B = B, L = L, S = S, gamma = gamma, U = U))
}
block_cd_lowrank_single <- function(X, Y, Z, lam2, lam3, p, q, diagthres = 0.001, max_iter_L = 1000, 
    tol_L = 1e-04, max_outer = 1000, tol_outer = 1e-04) {
    n <- nrow(X)
    r <- ncol(Z)
    L_old <- matrix(0, nrow = p, ncol = q)
    gamma_old <- rep(0, r)
    Lf <- 2 * sum(rowSums(X^2))/n
    alpha_L <- 1/(2 * Lf)
    for (outer_iter in seq_len(max_outer)) {
        L_current <- L_old
        for (inner_L in seq_len(max_iter_L)) {
            vecL <- as.vector(L_current)
            residZ <- as.vector(Z %*% gamma_old)
            predL <- X %*% vecL
            r_vec <- Y - residZ - predL
            grad_vec <- -(2/n) * crossprod(X, r_vec)
            grad_L <- matrix(grad_vec, nrow = p, ncol = q)
            L_tilde <- L_current - alpha_L * grad_L
            svd_res <- svd(L_tilde)
            U <- svd_res$u
            V <- svd_res$v
            s <- svd_res$d
            s_new <- pmax(s - alpha_L * lam2, 0)
            L_next <- U %*% diag(s_new, length(s_new), length(s_new)) %*% t(V)
            diff_norm <- sqrt(sum((L_next - L_current)^2))
            rel_norm <- diff_norm/(1 + sqrt(sum(L_current^2)))
            if (rel_norm < tol_L) {
                L_current <- L_next
                break
            }
            L_current <- L_next
        }
        L_new <- L_current
        vecL <- as.vector(L_new)
        predL <- X %*% vecL
        r_vec_g <- Y - predL
        gamma_new <- coordinate_descent_lasso(y = r_vec_g, Z = Z, lambda = lam3, tol = 1e-06, max_iter = 1000, 
            eps = 1e-08)
        delta_L <- sqrt(sum((L_new - L_old)^2))
        rel_L <- delta_L/(1 + sqrt(sum(L_old^2)))
        delta_gamma <- sqrt(sum((gamma_new - gamma_old)^2))/(1 + sqrt(sum(gamma_old^2)))
        if (rel_L < tol_outer && delta_gamma < tol_outer) {
            L_old <- L_new
            gamma_old <- gamma_new
            break
        }
        L_old <- L_new
        gamma_old <- gamma_new
    }
    svdL_old <- svd(L_old)
    oldLD_thresholded <- pmax(svdL_old$d - diagthres, 0)
    L_mat <- svdL_old$u %*% diag(oldLD_thresholded, nrow = length(oldLD_thresholded)) %*% t(svdL_old$v)
    L_old <- as.vector(L_mat)
    list(L = L_old, gamma = gamma_old)
}
list(fit = linear_admm, il = block_cd_lowrank_single)
})
