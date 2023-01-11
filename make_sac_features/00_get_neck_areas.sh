#!/bin/bash
# SLURM submission script for multiple serial jobs on Niagara
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:15:00
#SBATCH -p debug
#SBATCH --job-name neck
#SBATCH --output=neck_%j.txt
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export NUMEXPR_NUM_THREADS=2

export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true
export MPLCONFIGDIR=/scratch/s/steinman/macdo708/.config/matplotlib
export NUMBA_CACHE_DIR=/scratch/s/steinman/macdo708/.config/numba
 
module load CCEnv StdEnv/2020 gcc/9.3.0 vtk/9.0.1 python/3.7.7
module load openmpi/4.0.3

source ~/.virtualenvs/vtk9/bin/activate

SEQ="seq_005"

MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_analysis"
PROJ_FOLDER="/scratch/s/steinman/macdo708/neuromorph_analysis/$SEQ"
RESULTS_FOLDER="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/results"
MESHING_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/meshing"

cd $MAIN_DIR

FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )

############
SURF_OUT_FOLDER=$PROJ_FOLDER/data/neck_surfs
CSV_OUT_FOLDER=$PROJ_FOLDER/data/neck_areas

mkdir -p $SURF_OUT_FOLDER
mkdir -p $CSV_OUT_FOLDER

parallel -j40 python ./make_sac_features/00_get_neck_area.py {} $SURF_OUT_FOLDER $CSV_OUT_FOLDER $MESHING_DIR ::: $FOLDERS

# python ./make_sac_features/01_merge_neck_areas.py