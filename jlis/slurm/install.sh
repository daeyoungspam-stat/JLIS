#!/bin/bash
#SBATCH --job-name=jlis_install
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export R_LIBS_USER="${R_LIBS_USER:-$HOME/Rlibs}"
mkdir -p "$R_LIBS_USER"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
if [[ -n "${JLIS_R_MODULE:-}" ]]; then
  module load "$JLIS_R_MODULE"
fi
Rscript --vanilla scripts/install_dependencies.R
R CMD INSTALL --library="$R_LIBS_USER" .
