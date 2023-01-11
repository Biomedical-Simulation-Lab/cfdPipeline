#!/bin/bash
# SLURM submission script for multiple serial jobs on Niagara
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=01:00:00
#SBATCH -p debug
#SBATCH --job-name band_
#SBATCH --output=band_%j.txt
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true
export MPLCONFIGDIR=/scratch/s/steinman/macdo708/.config/matplotlib
export NUMBA_CACHE_DIR=/scratch/s/steinman/macdo708/.config/numba

module load NiaEnv/2019b intelpython3 gnu-parallel
source activate aneurisk_librosa

SEQ="seq_005"

MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_analysis"
PROJ_FOLDER="/scratch/s/steinman/macdo708/neuromorph_analysis/$SEQ"

HEMO_DATA_FOLDER=$PROJ_FOLDER/data/params_from_tec
SBI_DATA_FOLDER=$PROJ_FOLDER/data/SBI

mkdir - p $SBI_DATA_FOLDER

cd $MAIN_DIR

python ./make_spectrograms/02_calc_bandedness.py $PROJ_FOLDER $HEMO_DATA_FOLDER $SBI_DATA_FOLDER