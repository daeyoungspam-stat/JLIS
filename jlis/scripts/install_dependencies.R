packages <- c("pracma", "ggplot2", "glmnet", "foreach", "doParallel", "iterators", "MASS")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
