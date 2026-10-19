#!/bin/bash
#SBATCH --job-name=jlis_legacy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=3-00:00:00
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
setting="${1:?Pass a setting such as A1 or C1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export R_LIBS_USER="${R_LIBS_USER:-$HOME/Rlibs}"
if [[ -n "${JLIS_R_MODULE:-}" ]]; then
  module load "$JLIS_R_MODULE"
fi
Rscript --vanilla scripts/reproduce.R \
  --mode legacy --setting "$setting" --cores "${SLURM_CPUS_PER_TASK:-1}" \
  --outdir "results/$setting/legacy"
