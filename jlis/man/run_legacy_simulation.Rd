\name{run_legacy_simulation}
\alias{run_legacy_simulation}
\title{Execute the Preserved Comparison Experiment}
\description{Runs the preserved JLIS, ILIS, IL, lasso, and ridge experiment for a simulation setting.}
\usage{
run_legacy_simulation(setting = "A1",
  outdir = file.path("results", setting, "legacy"), cores = 1L, reps = 50L)
}
\arguments{
\item{setting}{Simulation setting identifier.}
\item{outdir}{A fresh output directory.}
\item{cores}{Number of parallel tuning workers.}
\item{reps}{Number of sequential experiment repetitions, at most 50.}
}
\details{The historical sequential score accumulation remains part of the preserved comparison experiment. Dependencies include glmnet, foreach, doParallel, iterators, MASS, and pracma.}
\value{Result-file paths, invisibly.}
