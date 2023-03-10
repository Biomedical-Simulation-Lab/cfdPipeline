#!/bin/bash
# SLURM submission script for multiple serial jobs on Niagara
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=01:00:00
#SBATCH --job-name spec_
#SBATCH --output=spec_%j.txt
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true
export MPLCONFIGDIR=/scratch/s/steinman/macdo708/.config/matplotlib
export NUMBA_CACHE_DIR=/scratch/s/steinman/macdo708/.config/numba
 
module load CCEnv StdEnv/2020 gcc/9.3.0 vtk/9.0.1 python/3.7.7
module load openmpi/4.0.3

source ~/.virtualenvs/vtk9/bin/activate

SEQ="seq_006"

MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_analysis"
PROJ_FOLDER="/scratch/s/steinman/macdo708/neuromorph_analysis/$SEQ"
RESULTS_FOLDER="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/results"
MESHING_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/meshing"

# PROJ_FOLDER="/scratch/s/steinman/macdo708/surge_cfd_analysis"
OUT_DATA_FOLDER=$PROJ_FOLDER/spectrogram_data
OUT_IMG_FOLDER=$PROJ_FOLDER/spectrogram_imgs

cd $MAIN_DIR

# RESULTS_FOLDER="/scratch/s/steinman/macdo708/surge_cfd/results"
FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )
parallel -j20 python ./make_spectrograms/00_compute_spectrograms.py {} $OUT_DATA_FOLDER $OUT_IMG_FOLDER $MESHING_DIR ::: $FOLDERS

# python ./make_spectrograms/01_merge_spectrograms.py
# ~/xvfb-run-safe python ./make_figures/00_make_fig.py