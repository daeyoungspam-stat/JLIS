.dgp_logistic <- local({
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
    return(list(G = G, n = n, p = p, q = q, r = r, L = L, S_list = S_list, B_list = B_list, gamma_list = gamma_list, 
        X_list = X_list, Z_list = Z_list, Y_list = Y_list))
}
list(generate = simulate_data_diff_logit)
})
