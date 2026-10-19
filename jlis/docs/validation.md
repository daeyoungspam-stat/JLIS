# Validation

The validation scope comprises source parity, reduced execution tests, archive
comparison, figure generation, and package checks.

| Check | Result |
|---|---|
| Final-source data components, eight settings | 56 exact comparisons passed |
| Final-source estimator states, three configurations | 15 exact comparisons passed |
| Reduced A1 and C1 experiments | Selected parameters, validation loss, coefficient errors, and ranks matched |
| Archived simulation files | 400 replicate files validated |
| Combined coefficient errors | 9,600 rows matched the completed result archives |
| Gaussian and binomial comparison estimators | Fixed-penalty JLIS, ILIS, IL, lasso, and ridge checks passed |
| Public fit and reproduction API | Dimension checks, entrywise feasibility, predictions, resume, and combination checks passed |
| Figures | Eight PDF and eight PNG files generated |
| Package installation and checks | R 4.6.1 package check completed successfully |
| R 4.5.1 execution | Core tests and comparison-estimator checks passed |

The numerical archive supplies the 50-replicate reference results. Execution
validation uses reduced experiments and fixed-penalty fits. A1/A2 archive
comparisons represent marginal error distributions: the constrained and
comparison scripts generate distinct realizations of part of the data.

`inst/extdata/source_manifest.csv` records source-file checksums.
`inst/extdata/dgp_reference.rds` contains data snapshots for the regression tests.
