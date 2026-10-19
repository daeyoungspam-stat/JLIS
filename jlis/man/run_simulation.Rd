\name{run_simulation}
\alias{run_simulation}
\title{Execute the JLIS-Alpha Simulation}
\description{Fits the final constrained estimator and writes one result file per replicate.}
\usage{
run_simulation(setting = "A1", replicates = 1:50,
  outdir = file.path("results", setting), cores = 1L,
  profile = c("publication", "smoke"), resume = TRUE)
}
\arguments{
\item{setting}{Simulation setting identifier.}
\item{replicates}{Distinct replication positions between 1 and 50.}
\item{outdir}{Output directory.}
\item{cores}{Number of tuning workers.}
\item{profile}{Final simulation or small execution-test configuration.}
\item{resume}{Validate and reuse existing result files.}
}
\value{Result-file paths, invisibly.}
\details{Result files contain configuration metadata, selected penalties, coefficient errors, and fit diagnostics.}
