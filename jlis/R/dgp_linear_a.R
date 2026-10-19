.dgp_linear_a <- local({
make_sparse_matrix2 <- function(p, q, sparsity, unif_min, unif_max) {
    num_nonzero <- ceiling((1 - sparsity) * p * q)
    mat <- matrix(0, nrow = p, ncol = q)
    idx <- sample.int(p * q, size = num_nonzero, replace = FALSE)
    signs <- 2 * rbinom(n = num_nonzero, size = 1, prob = 0.5) - 1
    mat[idx] <- runif(num_nonzero, min = unif_min, max = unif_max) * signs
    mat
}
simulate_data1 <- function(G, n, p, q, r, signal = 0.7, gamma = 0.25, scale, rank_L, sparsegam, singular_vals = c(0.5, 
    0.1, 0.05), umingam, umaxgam, sigma = 1, fold) {
    L_true <- matrix(0, p, q)
    if (p == q) {
        L_basis <- pracma::randortho(p)
        for (j in seq_len(rank_L)) {
            L_true <- L_true + singular_vals[j] * (L_basis[, j] %*% t(L_basis[, j]))
        }
    }
    else {
        L_basis <- pracma::randortho(p)
        L_basis2 <- pracma::randortho(q)
        for (j in seq_len(rank_L)) {
            L_true <- L_true + singular_vals[j] * (L_basis[, j] %*% t(L_basis2[, j]))
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
    max_eigen <- numeric(G)
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
    for (g in seq_len(G)) {
        gamma_list[[g]] <- make_sparse_matrix2(p = r, q = 1, sparsity = sparsegam, unif_min = umingam, 
            unif_max = umaxgam)
    }
    X_list <- vector("list", G)
    Z_list <- vector("list", G)
    Y_list <- vector("list", G)
    for (g in seq_len(G)) {
        Xg <- matrix(rnorm(n[g] * fold * p * q), nrow = n[g] * fold, ncol = p * q)
        Zg <- matrix(rnorm(r * n[g] * fold), nrow = n[g] * fold, ncol = r)
        noise <- matrix(rnorm(n[g] * fold, mean = 0, sd = sigma), nrow = n[g] * fold, ncol = 1)
        Yg <- Xg %*% B_list[[g]] + Zg %*% gamma_list[[g]] + noise
        X_list[[g]] <- Xg
        Z_list[[g]] <- Zg
        Y_list[[g]] <- Yg
    }
    list(G = G, n = n, p = p, q = q, r = r, L = L, S_list = S_list, B_list = B_list, gamma_list = gamma_list, 
        X_list = X_list, Z_list = Z_list, Y_list = Y_list)
}
list(generate = simulate_data1)
})
