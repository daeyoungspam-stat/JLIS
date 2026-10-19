.legacy_logistic <- local({
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
list(fit = logistic_admm, il = block_cd_logistic_single)
})
