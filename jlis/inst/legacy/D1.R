library(foreach)
library(doParallel)
library(glmnet)
library(iterators)
library(MASS)
library(pracma)
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
logistic_admm <- function(Y_list, X_list, Z_list, px, qx, lambda1, lambda2, lambda3, omega1, omega2, 
    rho = 1, max_iter = 1000, tol = 1e-04, rankthres = 0.001) {
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
        BGam_updates <- lapply(1:G, function(g) {
            Xg <- X_list[[g]]
            Yg <- Y_list[[g]]
            Zg <- Z_list[[g]]
            Sg <- S[[g]]
            Ug <- U[[g]]
            out <- pgd_logit(X = Xg, Z = Zg, y = Yg, L_vec = L, S_vec = Sg, U_vec = Ug, lambda3 = lambda3, 
                omega2g = omega2[g], rho = rho, B_init = B[[g]], gamma_init = gamma[[g]], t_init = NULL, 
                eta_bt = 2, max_inner = 1000, tol_inner = 1e-04)
            out
        })
        B <- lapply(BGam_updates, `[[`, "B")
        gamma <- lapply(BGam_updates, `[[`, "gamma")
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
sigmoid <- function(eta) {
    out <- eta
    pos <- eta >= 0
    out[pos] <- 1/(1 + exp(-eta[pos]))
    e <- exp(eta[!pos])
    out[!pos] <- e/(1 + e)
    out
}
logistic_loss_vec <- function(eta, y) {
    mean(softplus_stable(eta) - y * eta)
}
prox_l1 <- function(v, tlam) {
    sign(v) * pmax(abs(v) - tlam, 0)
}
prox_nuclear <- function(M, tlam) {
    sv <- svd(M)
    s <- pmax(sv$d - tlam, 0)
    if (all(s == 0)) 
        return(matrix(0, nrow = nrow(M), ncol = ncol(M)))
    sv$u %*% diag(s, length(s)) %*% t(sv$v)
}
block_cd_logistic_single <- function(X, Y, Z, lam2, lam3, p, q, diagthres = 0.001, max_outer = 1000, 
    tol_outer = 1e-04, max_inner_B = 1000, max_inner_g = 1000, use_backtracking = TRUE) {
    n <- nrow(X)
    r <- if (is.null(Z)) 
        0L
    else ncol(Z)
    B <- matrix(0, p, q)
    g <- if (r > 0) 
        rep(0, r)
    else numeric(0)
    LB <- 0.25 * sum(X^2)/n
    stepB <- 1/(LB * 2)
    Lg <- 0.25 * sum(Z^2)/n
    stepg <- 1/(Lg * 2)
    for (outer in seq_len(max_outer)) {
        B_old <- B
        g_old <- g
        for (it in seq_len(max_inner_B)) {
            eta <- as.vector(X %*% as.vector(B))
            if (r > 0) 
                eta <- eta + as.vector(Z %*% g)
            pvec <- sigmoid(eta)
            resid <- (pvec - Y)/n
            grad_vec <- crossprod(X, resid)
            gradB <- matrix(grad_vec, p, q)
            tB <- stepB
            repeat {
                B_try <- prox_nuclear(B - tB * gradB, tB * lam2)
                if (!use_backtracking) 
                  break
                eta_new <- as.vector(X %*% as.vector(B_try))
                if (r > 0) 
                  eta_new <- eta_new + as.vector(Z %*% g)
                f_old <- logistic_loss_vec(eta, Y)
                f_new <- logistic_loss_vec(eta_new, Y)
                quad <- sum(gradB * (B_try - B)) + (1/(2 * tB)) * sum((B_try - B)^2)
                if (f_new <= f_old + quad + 1e-12) 
                  break
                tB <- tB * 0.5
                if (tB < 1e-12) 
                  break
            }
            rel <- sqrt(sum((B_try - B)^2))/(1e-12 + sqrt(sum(B^2)))
            B <- B_try
            if (rel < 1e-04) 
                break
        }
        if (r > 0) {
            for (it in seq_len(max_inner_g)) {
                eta <- as.vector(X %*% as.vector(B)) + as.vector(Z %*% g)
                pvec <- sigmoid(eta)
                resid <- (pvec - Y)/n
                gradg <- as.vector(crossprod(Z, resid))
                tg <- stepg
                repeat {
                  g_try <- prox_l1(g - tg * gradg, tg * lam3)
                  if (!use_backtracking) 
                    break
                  eta_old <- eta
                  eta_new <- as.vector(X %*% as.vector(B)) + as.vector(Z %*% g_try)
                  f_old <- logistic_loss_vec(eta_old, Y)
                  f_new <- logistic_loss_vec(eta_new, Y)
                  quad <- sum(gradg * (g_try - g)) + (1/(2 * tg)) * sum((g_try - g)^2)
                  if (f_new <= f_old + quad + 1e-12) 
                    break
                  tg <- tg * 0.5
                  if (tg < 1e-12) 
                    break
                }
                rel <- sqrt(sum((g_try - g)^2))/(1e-12 + sqrt(sum(g^2)))
                g <- g_try
                if (rel < 1e-04) 
                  break
            }
        }
        delta_B <- sqrt(sum((B - B_old)^2))
        rel_B <- delta_B/(1e-12 + sqrt(sum(B_old^2)))
        delta_g <- if (r > 0) 
            sqrt(sum((g - g_old)^2))
        else 0
        rel_g <- if (r > 0) 
            delta_g/(1e-12 + sqrt(sum(g_old^2)))
        else 0
        if (max(rel_B, rel_g) < tol_outer) {
            break
        }
    }
    if (diagthres > 0) {
        sv <- svd(B)
        s <- pmax(sv$d - diagthres, 0)
        B <- if (all(s == 0)) 
            matrix(0, p, q)
        else sv$u %*% diag(s, length(s)) %*% t(sv$v)
    }
    B <- as.vector(B)
    list(L = B, gamma = g)
}
make_sparse_matrix <- function(p, q, sparsity, unif_min, unif_max) {
    num_nonzero <- ceiling((1 - sparsity) * p * q)
    mat <- matrix(0, nrow = p, ncol = q)
    idx <- sample.int(p * q, size = num_nonzero, replace = FALSE)
    mat[idx] <- runif(num_nonzero, min = unif_min, max = unif_max)
    return(mat)
}
make_sparse_matrix2 <- function(p, q, sparsity, unif_min, unif_max) {
    num_nonzero <- ceiling((1 - sparsity) * p * q)
    mat <- matrix(0, nrow = p, ncol = q)
    idx <- sample.int(p * q, size = num_nonzero, replace = FALSE)
    mat[idx] <- runif(num_nonzero, min = unif_min, max = unif_max) * (2 * rbinom(n = num_nonzero, size = 1, 
        prob = 0.5) - 1)
    return(mat)
}
make_low_rank_matrix <- function(p, q, rk, unif_min, unif_max) {
    Up_full <- qr.Q(qr(matrix(rnorm(p * p), nrow = p, ncol = p)))
    Vq_full <- qr.Q(qr(matrix(rnorm(q * q), nrow = q, ncol = q)))
    sigmas <- runif(rk, min = unif_min, max = unif_max)
    L <- Up_full[, 1:rk] %*% diag(sigmas) %*% t(Vq_full[, 1:rk])
    return(L)
}
make_block_sparse_matrix <- function(p, q, block_size_row, block_size_col, num_nonzero_blocks, unif_min, 
    unif_max) {
    S <- matrix(0, nrow = p, ncol = q)
    nblocks_row <- floor(p/block_size_row)
    nblocks_col <- floor(q/block_size_col)
    total_blocks <- nblocks_row * nblocks_col
    if (num_nonzero_blocks > total_blocks) {
        stop("Error!!! due to num_nonzero_blocks > total_blocks")
    }
    selected_indices <- sample(1:total_blocks, num_nonzero_blocks, replace = FALSE)
    for (index in selected_indices) {
        row_block_index <- (index - 1)%/%nblocks_col + 1
        col_block_index <- (index - 1)%%nblocks_col + 1
        row_start <- (row_block_index - 1) * block_size_row + 1
        row_end <- row_block_index * block_size_row
        col_start <- (col_block_index - 1) * block_size_col + 1
        col_end <- col_block_index * block_size_col
        S[row_start:row_end, col_start:col_end] <- matrix(runif(block_size_row * block_size_col, min = unif_min, 
            max = unif_max), nrow = block_size_row, ncol = block_size_col)
    }
    return(S)
}
make_block_sparse_matrix2 <- function(p, q, block_size_row, block_size_col, num_nonzero_blocks, unif_min, 
    unif_max) {
    S <- matrix(0, nrow = p, ncol = q)
    nblocks_row <- floor(p/block_size_row)
    nblocks_col <- floor(q/block_size_col)
    total_blocks <- nblocks_row * nblocks_col
    if (num_nonzero_blocks > total_blocks) {
        stop("Error!!! due to num_nonzero_blocks > total_blocks")
    }
    selected_indices <- sample(1:total_blocks, num_nonzero_blocks, replace = FALSE)
    for (index in selected_indices) {
        row_block_index <- (index - 1)%/%nblocks_col + 1
        col_block_index <- (index - 1)%%nblocks_col + 1
        row_start <- (row_block_index - 1) * block_size_row + 1
        row_end <- row_block_index * block_size_row
        col_start <- (col_block_index - 1) * block_size_col + 1
        col_end <- col_block_index * block_size_col
        S[row_start:row_end, col_start:col_end] <- matrix(runif(block_size_row * block_size_col, min = unif_min, 
            max = unif_max) * (2 * rbinom(n = block_size_row * block_size_col, size = 1, prob = 0.5) - 
            1), nrow = block_size_row, ncol = block_size_col)
    }
    return(S)
}
simulate_data_logit <- function(G, n, p, q, r, signal = 0.7, gamma = 0.25, scale, rank_L, sparsegam, 
    singular_vals = c(0.5, 0.1, 0.05), umingam, umaxgam, fold) {
    L_true <- matrix(0, p, q)
    if (p == q) {
        L_basis <- randortho(p)
        for (j in 1:rank_L) {
            L_true = L_true + singular_vals[j] * (L_basis[, j] %*% t(L_basis[, j]))
        }
    }
    else {
        L_basis <- randortho(p)
        L_basis2 <- randortho(q)
        for (j in 1:rank_L) {
            L_true = L_true + singular_vals[j] * (L_basis[, j] %*% t(L_basis2[, j]))
        }
    }
    L_true <- (signal * gamma/max(L_true)) * L_true
    L <- as.vector(L_true)
    S_list <- vector("list", G)
    B_list <- vector("list", G)
    for (g in seq_len(G)) {
        S_mat <- matrix(0, p, q)
        off <- (-1)^(g - 1) * ceiling(g/2)
        for (j in seq_len(p)) {
            col <- j + off
            if (col >= 1 && col <= q) {
                S_mat[j, col] <- (-1)^g * signal
            }
        }
        S_list[[g]] <- as.vector(S_mat)
        B_list[[g]] <- L + S_list[[g]]
    }
    max_eigen = numeric(G)
    for (g in seq_len(G)) {
        max_eigen[g] <- max(abs(eigen(matrix(B_list[[g]], p, q))$values))
    }
    grand_max_eigen <- max(max_eigen)
    for (g in seq_len(G)) {
        B_list[[g]] <- B_list[[g]] * scale/grand_max_eigen
        S_list[[g]] <- S_list[[g]] * scale/grand_max_eigen
    }
    L <- L * scale/grand_max_eigen
    gamma_list <- vector("list", G)
    for (g in 1:G) {
        gamma_list[[g]] <- make_sparse_matrix2(p = r, q = 1, sparsity = sparsegam, unif_min = umingam, 
            unif_max = umaxgam)
    }
    X_list <- vector("list", G)
    Z_list <- vector("list", G)
    Y_list <- vector("list", G)
    for (g in 1:G) {
        Xg <- matrix(rnorm(n[g] * fold * p * q), nrow = n[g] * fold, ncol = p * q)
        Zg <- matrix(rnorm(r * n[g] * fold), nrow = n[g] * fold, ncol = r)
        Bg <- B_list[[g]]
        gammag <- gamma_list[[g]]
        eta_g <- Xg %*% Bg + Zg %*% gammag
        p_g <- 1/(1 + exp(-eta_g))
        Yg <- rbinom(n[g] * fold, size = 1, prob = p_g)
        Yg <- matrix(Yg, ncol = 1)
        X_list[[g]] <- Xg
        Z_list[[g]] <- Zg
        Y_list[[g]] <- Yg
    }
    return(list(G = G, n = ngs, p = p, q = q, r = r, L = L, S_list = S_list, B_list = B_list, gamma_list = gamma_list, 
        X_list = X_list, Z_list = Z_list, Y_list = Y_list))
}
simulate_data_diff_logit <- function(G, n, p, q, r, rank_L, sparsegam, uminL, umaxL, uminS, umaxS, umingam, 
    umaxgam, Sblockrow, Sblockcol, Sblocknum, fold) {
    L.mat <- make_low_rank_matrix(p = p, q = q, rk = rank_L, unif_min = uminL, unif_max = umaxL)
    L <- as.vector(L.mat)
    S_list <- vector("list", G)
    B_list <- vector("list", G)
    for (g in 1:G) {
        S.mat <- make_block_sparse_matrix(p = p, q = q, block_size_row = Sblockrow, block_size_col = Sblockcol, 
            num_nonzero_blocks = Sblocknum, unif_min = uminS, unif_max = umaxS)
        S_list[[g]] <- as.vector(S.mat)
        B_list[[g]] <- L + S_list[[g]]
    }
    gamma_list <- vector("list", G)
    for (g in 1:G) {
        gamma_list[[g]] <- make_sparse_matrix2(p = r, q = 1, sparsity = sparsegam, unif_min = umingam, 
            unif_max = umaxgam)
    }
    X_list <- vector("list", G)
    Z_list <- vector("list", G)
    Y_list <- vector("list", G)
    for (g in 1:G) {
        Xg <- matrix(rnorm(n[g] * fold * p * q), nrow = n[g] * fold, ncol = p * q)
        Zg <- matrix(rnorm(r * n[g] * fold), nrow = n[g] * fold, ncol = r)
        Bg <- B_list[[g]]
        gammag <- gamma_list[[g]]
        eta_g <- Xg %*% Bg + Zg %*% gammag
        p_g <- 1/(1 + exp(-eta_g))
        Yg <- rbinom(n[g] * fold, size = 1, prob = p_g)
        Yg <- matrix(Yg, ncol = 1)
        X_list[[g]] <- Xg
        Z_list[[g]] <- Zg
        Y_list[[g]] <- Yg
    }
    return(list(G = G, n = ngs, p = p, q = q, r = r, L = L, S_list = S_list, B_list = B_list, gamma_list = gamma_list, 
        X_list = X_list, Z_list = Z_list, Y_list = Y_list))
}
lam1.vec = 10^(seq(from = -3, to = 3, by = 0.5))
lam2.vec = 10^(seq(from = -3, to = 3, by = 0.5))
lam3.vec = 10^(seq(from = -3, to = 3, by = 0.5))
total.length = length(lam1.vec) * length(lam2.vec) * length(lam3.vec)
param_grid <- expand.grid(i1 = 1:length(lam1.vec), i2 = 1:length(lam2.vec), i3 = 1:length(lam3.vec))
grid23 <- expand.grid(i2 = seq_along(lam2.vec), i3 = seq_along(lam3.vec))
reps <- 200
num.Groups <- 4
num.Folds <- 5
ngs = rep(60, num.Groups)
rowx <- 30
colx <- 30
result.vec = numeric(total.length)
sum.result.vec = numeric(total.length)
result.vec.lasso = matrix(0, nrow = num.Groups, ncol = length(lam1.vec))
sum.result.vec.lasso = matrix(0, nrow = num.Groups, ncol = length(lam1.vec))
result.vec.ridge = matrix(0, nrow = num.Groups, ncol = length(lam1.vec))
sum.result.vec.ridge = matrix(0, nrow = num.Groups, ncol = length(lam1.vec))
Final_main_lambdas = matrix(NA, nrow = reps, ncol = 3)
Final_results_B = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_gamma = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_rank = rep(NA, reps)
Final_lasso_lambdas = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_B_lasso = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_gamma_lasso = matrix(NA, nrow = reps, ncol = num.Groups)
Final_ridge_lambdas = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_B_ridge = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_gamma_ridge = matrix(NA, nrow = reps, ncol = num.Groups)
Final_onlylow_lambdas = matrix(NA, nrow = reps, ncol = 2 * num.Groups)
Final_results_B_onlylow = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_gamma_onlylow = matrix(NA, nrow = reps, ncol = num.Groups)
error_array <- array(0, dim = c(num.Groups, length(lam2.vec), length(lam3.vec)), dimnames = list(paste0("Group", 
    seq_len(num.Groups)), paste0("i2_", seq_len(length(lam2.vec))), paste0("i3_", seq_len(length(lam3.vec)))))
sum_error_array <- array(0, dim = c(num.Groups, length(lam2.vec), length(lam3.vec)), dimnames = list(paste0("Group", 
    seq_len(num.Groups)), paste0("i2_", seq_len(length(lam2.vec))), paste0("i3_", seq_len(length(lam3.vec)))))
Final_main_lambdas.sep = matrix(NA, nrow = reps, ncol = 3 * num.Groups)
Final_results_B.sep = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_gamma.sep = matrix(NA, nrow = reps, ncol = num.Groups)
Final_results_rank.sep = matrix(NA, nrow = reps, ncol = num.Groups)
error_array.sep <- array(0, dim = c(num.Groups, length(lam1.vec), length(lam2.vec), length(lam3.vec)), 
    dimnames = list(paste0("Group", seq_len(num.Groups)), paste0("i1_", seq_len(length(lam1.vec))), paste0("i2_", 
        seq_len(length(lam2.vec))), paste0("i3_", seq_len(length(lam3.vec)))))
sum_error_array.sep <- array(0, dim = c(num.Groups, length(lam1.vec), length(lam2.vec), length(lam3.vec)), 
    dimnames = list(paste0("Group", seq_len(num.Groups)), paste0("i1_", seq_len(length(lam1.vec))), paste0("i2_", 
        seq_len(length(lam2.vec))), paste0("i3_", seq_len(length(lam3.vec)))))
bce_clipped <- function(p, y, eps = 1e-08) {
    y <- as.numeric(y)
    p <- pmin(pmax(p, eps), 1 - eps)
    -mean(y * log(p) + (1 - y) * log1p(-p))
}
set.seed(4321)
random_numbers <- sample(1:1000, size = 200, replace = TRUE)
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cores) || num_cores < 1) {
    num_cores <- 1
}
cl <- makeCluster(num_cores)
registerDoParallel(cl)
takes = c(1:9, 11:50, 51)
for (turn in head(takes, as.integer(Sys.getenv("JLIS_LEGACY_REPS", "50")))) {
    set.seed(random_numbers[turn])
    gen.data = simulate_data_diff_logit(G = 4, n = ngs, p = rowx, q = colx, r = 40, rank_L = 3, sparsegam = 0.7, 
        uminL = 0.75, umaxL = 1.5, uminS = 1e-04, umaxS = 0.01, umingam = 0.3, umaxgam = 0.6, Sblockrow = 3, 
        Sblockcol = 3, Sblocknum = 1, fold = num.Folds)
    for (flds in 1:num.Folds) {
        test.idx = (ngs[1] * (flds - 1) + 1):(ngs[1] * flds)
        Y_train = lapply(gen.data$Y_list, function(mat) mat[-test.idx, ])
        Y_test = lapply(gen.data$Y_list, function(mat) mat[test.idx, ])
        X_train = lapply(gen.data$X_list, function(mat) mat[-test.idx, ])
        X_test = lapply(gen.data$X_list, function(mat) mat[test.idx, ])
        Z_train = lapply(gen.data$Z_list, function(mat) mat[-test.idx, ])
        Z_test = lapply(gen.data$Z_list, function(mat) mat[test.idx, ])
        results_parallel <- foreach(idx = iter(param_grid, by = "row")) %dopar% {
            logistic_admm(Y_train, X_train, Z_train, lambda1 = lam1.vec[as.numeric(idx$i1)], lambda2 = lam2.vec[as.numeric(idx$i2)], 
                lambda3 = lam3.vec[as.numeric(idx$i3)], rho = 1, max_iter = 1000, tol = 1e-04, px = rowx, 
                qx = colx)
        }
        results_parallel_lasso <- foreach(gps = 1:num.Groups, .packages = "glmnet") %:% foreach(idx = seq_along(lam1.vec)) %dopar% 
            {
                lambda_val <- lam1.vec[idx]/nrow(X_train[[gps]])
                glmnet(x = cbind(X_train[[gps]], Z_train[[gps]]), y = Y_train[[gps]], alpha = 1, lambda = lambda_val, 
                  family = "binomial", intercept = FALSE)
            }
        results_parallel_ridge <- foreach(gps = 1:num.Groups, .packages = "glmnet") %:% foreach(idx = seq_along(lam1.vec)) %dopar% 
            {
                lambda_val <- lam1.vec[idx]/nrow(X_train[[gps]])
                glmnet(x = cbind(X_train[[gps]], Z_train[[gps]]), y = Y_train[[gps]], alpha = 0, lambda = lambda_val, 
                  family = "binomial", intercept = FALSE)
            }
        result.vec <- sapply(results_parallel, function(res) {
            nll_list <- sapply(1:num.Groups, function(i) {
                eta_i <- as.vector(X_test[[i]] %*% res$B[[i]] + Z_test[[i]] %*% res$gamma[[i]])
                p_i <- 1/(1 + exp(-eta_i))
                nll <- bce_clipped(p_i, Y_test[[i]])
            })
            mean(nll_list)
        })
        sum.result.vec <- sum.result.vec + result.vec
        result.vec.lasso <- sapply(seq_along(lam1.vec), function(idx) {
            sapply(1:num.Groups, function(i) {
                beta <- as.vector(coef(results_parallel_lasso[[i]][[idx]]))
                eta_i <- cbind(X_test[[i]], Z_test[[i]]) %*% beta[-1]
                p_i <- 1/(1 + exp(-eta_i))
                nll <- bce_clipped(p_i, Y_test[[i]])
            })
        })
        sum.result.vec.lasso <- sum.result.vec.lasso + result.vec.lasso
        result.vec.ridge <- sapply(seq_along(lam1.vec), function(idx) {
            sapply(1:num.Groups, function(i) {
                beta <- as.vector(coef(results_parallel_ridge[[i]][[idx]]))
                eta_i <- cbind(X_test[[i]], Z_test[[i]]) %*% beta[-1]
                p_i <- 1/(1 + exp(-eta_i))
                nll <- bce_clipped(p_i, Y_test[[i]])
            })
        })
        sum.result.vec.ridge <- sum.result.vec.ridge + result.vec.ridge
        for (gg in seq_len(num.Groups)) {
            results_parallel <- foreach(idx = iter(grid23, by = "row")) %dopar% {
                block_cd_logistic_single(Y = Y_train[[gg]], X = X_train[[gg]], Z = Z_train[[gg]], lam2 = lam2.vec[as.numeric(idx$i2)], 
                  lam3 = lam3.vec[as.numeric(idx$i3)], p = rowx, q = colx)
            }
            for (kk in seq_len(nrow(grid23))) {
                i2 <- grid23$i2[kk]
                i3 <- grid23$i3[kk]
                fit_k <- results_parallel[[kk]]
                eta_i <- as.vector(X_test[[gg]] %*% fit_k$L + Z_test[[gg]] %*% fit_k$gamma)
                p_i <- 1/(1 + exp(-eta_i))
                nll <- bce_clipped(p_i, Y_test[[gg]])
                error_array[gg, i2, i3] <- nll
            }
            sum_error_array <- sum_error_array + error_array
            results_parallel.ILIS <- foreach(idx = iter(param_grid, by = "row")) %dopar% {
                logistic_admm(Y_list = list(Y_train[[gg]]), X_list = list(X_train[[gg]]), Z_list = list(Z_train[[gg]]), 
                  lambda1 = lam1.vec[as.numeric(idx$i1)], lambda2 = lam2.vec[as.numeric(idx$i2)], lambda3 = lam3.vec[as.numeric(idx$i3)], 
                  px = rowx, qx = colx)
            }
            for (kk in seq_len(nrow(param_grid))) {
                i1 <- param_grid$i1[kk]
                i2 <- param_grid$i2[kk]
                i3 <- param_grid$i3[kk]
                fit_i <- results_parallel.ILIS[[kk]]
                eta_i <- as.vector(X_test[[gg]] %*% fit_i$B[[1]] + Z_test[[gg]] %*% fit_i$gamma[[1]])
                p_i <- 1/(1 + exp(-eta_i))
                nll <- bce_clipped(p_i, Y_test[[gg]])
                error_array.sep[gg, i1, i2, i3] <- nll
            }
        }
        sum_error_array.sep <- sum_error_array.sep + error_array.sep
        cat("fold number ", flds, " in the turn number ", turn, " is done now :)", "\n")
    }
    best_index <- which.min(sum.result.vec)
    best_params <- param_grid[best_index, ]
    Final_fit <- logistic_admm(gen.data$Y_list, gen.data$X_list, gen.data$Z_list, lambda1 = lam1.vec[best_params$i1], 
        lambda2 = lam2.vec[best_params$i2], lambda3 = lam3.vec[best_params$i3], rho = 1, max_iter = 1000, 
        tol = 1e-04, px = rowx, qx = colx)
    Final_B = Final_fit$B
    Final_gamma = Final_fit$gamma
    Final_main_lambdas[turn, ] = c(best_params$i1, best_params$i2, best_params$i3)
    Final_results_B[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((Final_B[[i]] - gen.data$B_list[[i]])^2)))
    Final_results_gamma[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((Final_gamma[[i]] - gen.data$gamma_list[[i]])^2)))
    Final_results_rank[turn] = qr(matrix(Final_fit$L, nrow = rowx, ncol = colx))$rank
    cat("selected lambdas indices are ", best_params$i1, " and ", best_params$i2, " and ", best_params$i3, 
        " the rank is ", Final_results_rank[turn])
    min_indices_lasso <- apply(sum.result.vec.lasso, 1, which.min)
    coeff_mat_lasso <- do.call(rbind, lapply(seq_len(num.Groups), function(gps) {
        fit <- glmnet(x = cbind(gen.data$X_list[[gps]], gen.data$Z_list[[gps]]), y = gen.data$Y_list[[gps]], 
            alpha = 1, family = "binomial", lambda = lam1.vec[min_indices_lasso[gps]]/nrow(gen.data$X_list[[gps]]), 
            intercept = FALSE)
        as.vector(coef(fit))[-1]
    }))
    Final_lasso_lambdas[turn, ] = min_indices_lasso
    Final_results_B_lasso[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((coeff_mat_lasso[i, 1:ncol(gen.data$X_list[[i]])] - 
        gen.data$B_list[[i]])^2)))
    Final_results_gamma_lasso[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((coeff_mat_lasso[i, 
        -(1:ncol(gen.data$X_list[[i]]))] - gen.data$gamma_list[[i]])^2)))
    min_indices_ridge <- apply(sum.result.vec.ridge, 1, which.min)
    coeff_mat_ridge <- do.call(rbind, lapply(seq_len(num.Groups), function(gps) {
        fit <- glmnet(x = cbind(gen.data$X_list[[gps]], gen.data$Z_list[[gps]]), y = gen.data$Y_list[[gps]], 
            alpha = 0, family = "binomial", lambda = lam1.vec[min_indices_ridge[gps]]/nrow(gen.data$X_list[[gps]]), 
            intercept = FALSE)
        as.vector(coef(fit))[-1]
    }))
    Final_ridge_lambdas[turn, ] = min_indices_ridge
    Final_results_B_ridge[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((coeff_mat_ridge[i, 1:ncol(gen.data$X_list[[i]])] - 
        gen.data$B_list[[i]])^2)))
    Final_results_gamma_ridge[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((coeff_mat_ridge[i, 
        -(1:ncol(gen.data$X_list[[i]]))] - gen.data$gamma_list[[i]])^2)))
    best_idx_mat <- matrix(NA, nrow = num.Groups, ncol = 2)
    colnames(best_idx_mat) <- c("best_i2", "best_i3")
    for (gg in seq_len(num.Groups)) {
        error_mat_g <- sum_error_array[gg, , ]
        current_min <- min(error_mat_g, na.rm = TRUE)
        min_pos <- which(error_mat_g == current_min, arr.ind = TRUE)
        best_i2 <- min_pos[1, 1]
        best_i3 <- min_pos[1, 2]
        best_idx_mat[gg, ] <- c(best_i2, best_i3)
        Final_onlylow_lambdas[turn, (gg * 2 - 1)] = best_i2
        Final_onlylow_lambdas[turn, (gg * 2)] = best_i3
    }
    best_lambda2 <- lam2.vec[best_idx_mat[, 1]]
    best_lambda3 <- lam3.vec[best_idx_mat[, 2]]
    list_of_results <- lapply(seq_len(num.Groups), function(gg) {
        lam2_opt <- best_lambda2[gg]
        lam3_opt <- best_lambda3[gg]
        fit_g <- block_cd_logistic_single(Y = gen.data$Y_list[[gg]], X = gen.data$X_list[[gg]], Z = gen.data$Z_list[[gg]], 
            lam2 = lam2_opt, lam3 = lam3_opt, p = rowx, q = colx)
        list(L = fit_g$L, gamma = fit_g$gamma)
    })
    onlylow_L_list <- lapply(list_of_results, "[[", "L")
    onlylow_gamma_list <- lapply(list_of_results, "[[", "gamma")
    Final_results_B_onlylow[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((onlylow_L_list[[i]] - 
        gen.data$B_list[[i]])^2)))
    Final_results_gamma_onlylow[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((onlylow_gamma_list[[i]] - 
        gen.data$gamma_list[[i]])^2)))
    best_idx_mat <- matrix(NA, nrow = num.Groups, ncol = 3)
    colnames(best_idx_mat) <- c("best_i1", "best_i2", "best_i3")
    for (gg in seq_len(num.Groups)) {
        error_mat_g <- sum_error_array.sep[gg, , , ]
        current_min <- min(error_mat_g, na.rm = TRUE)
        min_pos <- which(error_mat_g == current_min, arr.ind = TRUE)
        best_i1 <- min_pos[1, 1]
        best_i2 <- min_pos[1, 2]
        best_i3 <- min_pos[1, 3]
        best_idx_mat[gg, ] <- c(best_i1, best_i2, best_i3)
        Final_main_lambdas.sep[turn, (gg * 3 - 2)] = best_i1
        Final_main_lambdas.sep[turn, (gg * 3 - 1)] = best_i2
        Final_main_lambdas.sep[turn, (gg * 3)] = best_i3
    }
    best_lambda1 <- lam1.vec[best_idx_mat[, 1]]
    best_lambda2 <- lam2.vec[best_idx_mat[, 2]]
    best_lambda3 <- lam3.vec[best_idx_mat[, 3]]
    list_of_results <- lapply(seq_len(num.Groups), function(gg) {
        lam1_opt <- best_lambda1[gg]
        lam2_opt <- best_lambda2[gg]
        lam3_opt <- best_lambda3[gg]
        fit_g <- logistic_admm(Y_list = list(gen.data$Y_list[[gg]]), X_list = list(gen.data$X_list[[gg]]), 
            Z_list = list(gen.data$Z_list[[gg]]), lambda1 = lam1_opt, lambda2 = lam2_opt, lambda3 = lam3_opt, 
            px = rowx, qx = colx)
        list(B = fit_g$B[[1]], L = fit_g$L[[1]], gamma = fit_g$gamma[[1]])
    })
    sep_B_list <- lapply(list_of_results, "[[", "B")
    sep_L_list <- lapply(list_of_results, "[[", "L")
    sep_gamma_list <- lapply(list_of_results, "[[", "gamma")
    Final_results_B.sep[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((sep_B_list[[i]] - gen.data$B_list[[i]])^2)))
    Final_results_gamma.sep[turn, ] = sapply(1:num.Groups, function(i) sqrt(sum((sep_gamma_list[[i]] - 
        gen.data$gamma_list[[i]])^2)))
    Final_results_rank.sep[turn, ] = sapply(1:num.Groups, function(i) qr(matrix(sep_L_list[[i]], nrow = rowx, 
        ncol = colx))$rank)
    cat("turn ", turn, " is done now :)", "\n")
    save(Final_main_lambdas, Final_results_B, Final_results_gamma, Final_results_rank, Final_lasso_lambdas, 
        Final_results_B_lasso, Final_results_gamma_lasso, Final_ridge_lambdas, Final_results_B_ridge, 
        Final_results_gamma_ridge, Final_onlylow_lambdas, Final_results_gamma_onlylow, Final_results_B_onlylow, 
        Final_main_lambdas.sep, Final_results_B.sep, Final_results_gamma.sep, Final_results_rank.sep, 
        file = "simulation_results_Logistic_D1_main_high.RData")
}
stopCluster(cl)
