#!/bin/bash
# SLURM submission script for multiple serial jobs on Niagara
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=01:00:00
#SBATCH -p debug
#SBATCH --job-name specs
#SBATCH --output=spec_%j.txt
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8

export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true
export MPLCONFIGDIR=/scratch/s/steinman/macdo708/.config/matplotlib
export NUMBA_CACHE_DIR=/scratch/s/steinman/macdo708/.config/numba

module purge
module load NiaEnv/2019b intelpython3 gnu-parallel
source activate aneurisk_librosa

SPECS="/scratch/s/steinman/macdo708/pereira_rupture_study_analysis/spectrograms"

RESULTS_FOLDER="/scratch/s/steinman/macdo708/pereira_rupture_study/results"

SPECS_DATA=$SPECS/data
SPECS_IMGS=$SPECS/imgs

mkdir -p $SPECS_IMGS
mkdir -p $SPECS_DATA

FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d)

parallel -j10 python ./make_spectrograms/00_compute_spectrograms.py {} $SPECS_DATA $SPECS_IMGS ::: $FOLDERS