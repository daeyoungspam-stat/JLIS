\name{fit_jlis}
\alias{fit_jlis}
\alias{predict.jlis_fit}
\alias{coef.jlis_fit}
\alias{print.jlis_fit}
\title{Joint Low-Rank and Individual Sparse Regression}
\description{Fits a shared low-rank matrix and group-specific sparse matrix and vector coefficients.}
\usage{
fit_jlis(Y_list, X_list, Z_list, px, qx,
  lambda1, lambda2, lambda3, alpha = Inf,
  family = c("gaussian", "binomial"),
  omega1 = rep(1, length(Y_list)), omega2 = rep(1, length(Y_list)),
  control = list(), init = NULL)
\method{predict}{jlis_fit}(object, X_list, Z_list,
  type = c("response", "link"), ...)
\method{coef}{jlis_fit}(object, ...)
\method{print}{jlis_fit}(x, ...)
}
\arguments{
\item{Y_list}{List of numeric response vectors or one-column matrices.}
\item{X_list}{List of numeric matrices with one vectorized matrix predictor per row.}
\item{Z_list}{List of numeric covariate matrices with matching row counts.}
\item{px, qx}{Row and column dimensions of each matrix predictor.}
\item{lambda1, lambda2, lambda3}{Nonnegative penalties for the sparse matrices, common low-rank matrix, and vector coefficients.}
\item{alpha}{Positive entrywise-constraint parameter. The bound is \code{alpha / sqrt(px * qx)}.}
\item{family}{Gaussian or binomial response model.}
\item{omega1, omega2}{Nonnegative penalty weights, one per group.}
\item{control}{Named solver controls: \code{rho}, \code{max_iter}, \code{tol}, \code{dykstra_max_iter}, \code{dykstra_tol}, \code{svd}; Gaussian controls \code{gamma_max_iter}, \code{gamma_tol}; binomial controls \code{rankthres}, \code{pgd_max_inner}, \code{pgd_tol}.}
\item{init}{Optional initial state containing \code{B}, \code{L}, \code{S}, \code{gamma}, and \code{U}.}
\item{object, x}{A fitted \code{jlis_fit} object.}
\item{type}{Predictions on the response or linear-predictor scale.}
\item{...}{Additional S3 arguments.}
}
\details{Vectorization follows R column order. The Gaussian solver stores dense feature-space linear systems. Convergence diagnostics refer to the iteration stopping criterion before logistic post-thresholding.}
\value{\code{fit_jlis} returns a \code{jlis_fit} object with coefficients and diagnostics. \code{predict} returns one prediction vector per group. \code{coef} returns reshaped coefficient matrices and vector coefficients. \code{print} returns the fit invisibly.}
\examples{
d <- simulate_jlis("A1", profile = "smoke")
fit <- fit_jlis(d$Y_list, d$X_list, d$Z_list, d$p, d$q,
                0.1, 0.1, 0.1, alpha = 2,
                control = list(max_iter = 10L))
dim(coef(fit)$L)
length(predict(fit, d$X_list, d$Z_list))
}
