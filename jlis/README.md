# JLIS

Joint low-rank and individual sparse matrix regression for Gaussian and binary responses.
The package includes JLIS-alpha estimation, the A1-D2 simulation settings, completed
simulation results, and coefficient-error plots for JLIS, JLIS-alpha, ILIS, IL,
lasso, and ridge.

## Installation

From the repository directory:

```sh
Rscript scripts/install_dependencies.R
R CMD INSTALL .
```

After publication to a GitHub repository:

```r
remotes::install_github("REPOSITORY/jlis")
```

`REPOSITORY/jlis` denotes the repository location. Package metadata contains the
reserved contact address `package@example.invalid`; a project mailbox can replace
this field. License terms appear in `LICENSE`.

## Model

For group $g$, the linear predictor is

$$\eta_{ig}=\langle X_{ig},L+S_g\rangle+z_{ig}^{\mathsf T}\gamma_g.$$

The estimation problem is

$$\min_{L,S,\gamma}\;\sum_{g=1}^G\mathcal L_g(L+S_g,\gamma_g)
+\lambda_1\sum_g\omega_{1g}\|S_g\|_1
+\lambda_2\|L\|_*
+\lambda_3\sum_g\omega_{2g}\|\gamma_g\|_1,
\qquad\|L\|_\infty\le\alpha/\sqrt{pq}.$$

Here $\mathcal L_g=n_g^{-1}\sum_i(y_{ig}-\eta_{ig})^2$ for Gaussian responses
and $\mathcal L_g=n_g^{-1}\sum_i\{\log(1+e^{\eta_{ig}})-y_{ig}\eta_{ig}\}$
for binary responses. Setting `alpha = Inf` gives the unconstrained objective.
Tuning minimizes

$$\widehat\theta\in\operatorname*{arg\,min}_{\theta\in\Theta}
\sum_{k=1}^{5}\sum_{g=1}^{G}
\mathcal E\!\left(Y_{g,k},\widehat\eta_{g,k}^{(-k)}(\theta)\right),
\qquad\theta=(\lambda_1,\lambda_2,\lambda_3,\alpha).$$

The validation error $\mathcal E$ is squared-error sum for Gaussian responses
and mean binary cross-entropy for binary responses.

## Input and Estimation

Each list element represents a group. `X_list[[g]]` is an $n_g\times pq$ matrix,
`Z_list[[g]]` is an $n_g\times r$ matrix, and `Y_list[[g]]` contains $n_g$ responses.
Each image contributes `as.vector(image)` as one row of the design matrix.
`px` and `qx` specify its row and column dimensions.

```r
library(jlis)
d <- simulate_jlis("A1", replicate = 1, profile = "smoke")
fit <- fit_jlis(d$Y_list, d$X_list, d$Z_list,
                px = d$p, qx = d$q,
                lambda1 = 0.1, lambda2 = 0.1, lambda3 = 0.1,
                alpha = 2, family = "gaussian")
coef(fit)$L
predict(fit, d$X_list, d$Z_list)
fit$diagnostics
```

`family = "binomial"` selects logistic regression. The Gaussian implementation
stores dense $pq\times pq$ linear-system matrices. Its memory requirement grows
quadratically in the number of pixels.

## Simulation Settings

| Setting | Response | p = q | r | Rank of L | Observations per group |
|---|---|---:|---:|---:|---:|
| A1 | Gaussian | 20 | 25 | 2 | 200 |
| A2 | Gaussian | 20 | 25 | 3 | 200 |
| B1 | Gaussian | 40 | 40 | 4 | 1000 |
| B2 | Gaussian | 40 | 40 | 5 | 1000 |
| C1 | Binary | 15 | 25 | 2 | 100 |
| C2 | Binary | 15 | 25 | 3 | 100 |
| D1 | Binary | 30 | 40 | 3 | 300 |
| D2 | Binary | 30 | 40 | 4 | 300 |

Each setting contains four groups and 50 replicates. `profile = "publication"`
selects the final simulation configuration. `profile = "smoke"` selects a small
execution test with shortened iteration limits.

## Reproduction

Generate the eight figure sets from the bundled completed results:

```sh
Rscript scripts/reproduce.R --mode figures
```

Run the small execution examples:

```sh
Rscript scripts/smoke.R
```

Run JLIS-alpha and combine its output with the archived comparison results:

```sh
Rscript scripts/reproduce.R --mode run --setting A1 --cores 1
Rscript scripts/reproduce.R --mode combine --setting A1 --baseline reference
```

Outputs include per-replicate `.RData` files, `combined_results.RData`, numerical
summaries, and PDF/PNG figures. Existing replicate files are validated before reuse.

The preserved comparison scripts compute JLIS, ILIS, IL, lasso, and ridge:

```sh
Rscript scripts/reproduce.R --mode legacy --setting A1 --cores 1 --outdir results/A1/legacy
Rscript scripts/reproduce.R --mode combine --setting A1 \
  --baseline results/A1/legacy/simulation_results_Linear_A1_main_high.RData
```

The legacy path executes the preserved sequential experiment. Its historical
accumulation of validation scores across repetitions remains part of that path.
`run_simulation()` executes the final JLIS-alpha experiment. Completed numerical
results provide a fixed reference for the figures; fresh fits can vary with
numerical-library versions and finite solver tolerances.

Execution checks and archive interpretation appear in [validation](docs/validation.md).

## Slurm

From the repository directory, install the package in the cluster R library:

```sh
mkdir -p log
export JLIS_R_MODULE=R/4.5.1
sbatch --partition=compute1 slurm/install.sh
```

After the installation job finishes:

```sh
sbatch --partition=compute1 slurm/simulation.sh A1
```

The template requests 40 CPUs, all allocated-node memory, and three days per
array task, with up to ten concurrent tasks. Set the partition, R module, and
resource requests to values supported by the cluster. Each array task writes
one replicate. After completion:

```sh
Rscript scripts/reproduce.R --mode combine --setting A1 --baseline reference
```

`slurm/legacy.sh` runs the preserved comparison experiment as a single sequential
batch with parallel tuning fits. The R functions default to one worker.

## Package Checks

```sh
R CMD build .
R CMD check --no-manual jlis_0.1.0.tar.gz
```
