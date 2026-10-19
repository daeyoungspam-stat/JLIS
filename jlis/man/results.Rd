\name{reference_results}
\alias{reference_results}
\alias{combine_results}
\alias{summarize_results}
\alias{plot_results}
\alias{write_figures}
\title{Simulation Results and Figures}
\description{Loads archived results, combines completed replicate files, and produces coefficient-error summaries and plots.}
\usage{
reference_results(setting = NULL)
combine_results(outdir, setting = "A1", replicates = 1:50, baseline = NULL)
summarize_results(results = reference_results())
plot_results(results = reference_results(), setting = "A1")
write_figures(results = reference_results(), directory = "results/figures")
}
\arguments{
\item{setting}{Setting identifier; \code{NULL} loads all archived settings.}
\item{outdir}{Directory containing the \code{parts} subdirectory.}
\item{replicates}{Expected replication positions.}
\item{baseline}{\code{NULL}, \code{"reference"} for archived comparison results, or a legacy RData file.}
\item{results}{Result list containing \code{errors}, \code{ranks}, and \code{selection}.}
\item{directory}{Figure output directory.}
}
\value{\code{reference_results} and \code{combine_results} return result lists. \code{summarize_results} returns a data frame. \code{plot_results} returns a ggplot object. \code{write_figures} returns figure paths invisibly.}
\details{The archive contains 50 replicates per setting and six methods. Coefficient errors are Frobenius norms for matrices and Euclidean norms for vectors. \code{combine_results} also writes a combined RData file and CSV tables. Stored convergence diagnostics distinguish completed fitting from attainment of the stopping tolerance.}
\examples{
r <- reference_results("A1")
head(summarize_results(r))
}
