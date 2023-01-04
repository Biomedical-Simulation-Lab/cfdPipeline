#!/bin/bash
# SLURM submission script for multiple serial jobs on Niagara
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=01:00:00
#SBATCH -p debug
#SBATCH --job-name viz
#SBATCH --output=viz_%j.txt
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8

export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true
export MPLCONFIGDIR=/scratch/s/steinman/macdo708/.config/matplotlib
export NUMBA_CACHE_DIR=/scratch/s/steinman/macdo708/.config/numba
 
module load CCEnv StdEnv/2020 gcc/9.3.0 vtk/9.0.1 python/3.7.7
module load openmpi/4.0.3

source ~/.virtualenvs/vtk9/bin/activate

PROJ_FOLDER="/scratch/s/steinman/macdo708/pereira_rupture_study_analysis"

cd $PROJ_FOLDER

RESULTS_FOLDER="/scratch/s/steinman/macdo708/pereira_rupture_study/results"

# for data_type in "qcriterion_nd" "u_mag"; do
#     OUT_FOLDER=$PROJ_FOLDER/viz/$data_type
#     mkdir -p $OUT_FOLDER
#     FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )
#     parallel -j10 ~/xvfb-run-safe python ./make_viz/00_viz.py {} $OUT_FOLDER $data_type 8 ::: $FOLDERS
# done

# for data_type in "qcriterion_nd" "u_mag"; do
OUT_FOLDER=$PROJ_FOLDER/viz/u_mag
mkdir -p $OUT_FOLDER
FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )
parallel -j10 ~/xvfb-run-safe python ./make_viz/00_viz.py {} $OUT_FOLDER u_mag 8 ::: $FOLDERS
# done