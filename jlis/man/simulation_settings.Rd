\name{simulation_settings}
\alias{simulation_settings}
\alias{simulate_jlis}
\title{Simulation Settings and Generated Data}
\description{Provides the final A1-D2 JLIS-alpha experiment configurations and their generated matrix-regression data.}
\usage{
simulation_settings(setting = "A1", profile = c("publication", "smoke"))
simulate_jlis(setting = "A1", replicate = 1L,
  profile = c("publication", "smoke"))
}
\arguments{
\item{setting}{One of A1, A2, B1, B2, C1, C2, D1, or D2.}
\item{profile}{Final simulation configuration or small execution-test configuration.}
\item{replicate}{Replication position from 1 to 50.}
}
\value{\code{simulation_settings} returns a configuration list. \code{simulate_jlis} returns predictor and response lists, true coefficients, dimensions, and configuration metadata.}
\examples{
simulation_settings("C2")$p
d <- simulate_jlis("C2", profile = "smoke")
dim(d$X_list[[1]])
}
