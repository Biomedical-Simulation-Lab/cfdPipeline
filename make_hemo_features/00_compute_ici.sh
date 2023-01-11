#!/bin/bash
# SLURM submission script for multiple serial jobs on Niagara
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=01:00:00
#SBATCH --job-name ici_
#SBATCH --output=ici_%j.txt
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true
export MPLCONFIGDIR=/scratch/s/steinman/macdo708/.config/matplotlib
export NUMBA_CACHE_DIR=/scratch/s/steinman/macdo708/.config/numba
 
# module load CCEnv StdEnv/2020 gcc/9.3.0 vtk/9.0.1 python/3.7.7
# module load openmpi/4.0.3

# source ~/.virtualenvs/vtk9/bin/activate

module load NiaEnv/2019b intelpython3 gnu-parallel
source activate aneurisk_librosa

SEQ="seq_006"

MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_analysis"
PROJ_FOLDER="/scratch/s/steinman/macdo708/neuromorph_analysis/$SEQ"
RESULTS_FOLDER="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/results"
MESHING_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/meshing"

OUT_DATA_FOLDER=$PROJ_FOLDER/data/ici_data
mkdir -p $OUT_DATA_FOLDER

NECK_SURF_DIR="$PROJ_FOLDER/data/neck_surfs"
PARENT_SLICE_DIR="$PROJ_FOLDER/data/parent_slices"
mkdir -p $PARENT_SLICE_DIR

FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )

cd $MAIN_DIR

parallel -j20 python ./make_hemo_features/00_make_ici.py {} $OUT_DATA_FOLDER $NECK_SURF_DIR $MESHING_DIR $PARENT_SLICE_DIR ::: $FOLDERS

python ./make_hemo_features/01_merge_hemo_features.py $OUT_DATA_FOLDER
