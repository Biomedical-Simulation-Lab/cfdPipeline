#!/bin/bash
# SLURM submission script for multiple serial jobs on Niagara
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=01:00:00
#SBATCH -p debug
#SBATCH --job-name wss
#SBATCH --output=wss_%j.txt
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

PROJ_FOLDER="/scratch/s/steinman/macdo708/surge_cfd_analysis"
OUT_FOLDER=$PROJ_FOLDER/data/wss_params_per_case

cd $PROJ_FOLDER

RESULTS_FOLDER="/scratch/s/steinman/macdo708/surge_cfd/results"
FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )
parallel -j20 python ./make_wss/02_compute_wss_hemodynamics.py {} $OUT_FOLDER ::: $FOLDERS

python ./make_wss/03_merge_wss_params.py