#!/bin/bash
# SLURM submission script for multiple serial jobs on Niagara
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=11
#SBATCH --time=12:00:00
#SBATCH --job-name viz_q_all_sequences
#SBATCH --output=viz_%j.txt
#SBATCH --mail-type=FAIL

# Use this for debug: -p debug

export OMP_NUM_THREADS=7
export MKL_NUM_THREADS=7
export NUMEXPR_NUM_THREADS=7

export PYVISTA_OFF_SCREEN=true
export PYVISTA_USE_PANEL=true
export MPLCONFIGDIR=/scratch/s/steinman/macdo708/.config/matplotlib
export NUMBA_CACHE_DIR=/scratch/s/steinman/macdo708/.config/numba
 
module load CCEnv StdEnv/2020 gcc/9.3.0 vtk/9.0.1 python/3.7.7
module load openmpi/4.0.3

source ~/.virtualenvs/vtk9/bin/activate

################################################################################
NSTEPS=2
SEQ="seq_001"
VTYPE="u_mag"

MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_analysis"
PROJ_FOLDER="/scratch/s/steinman/macdo708/neuromorph_analysis/$SEQ"
RESULTS_FOLDER="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/results"
MESHING_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/meshing"
CPOS_FILE="/scratch/s/steinman/macdo708/neuromorph_analysis/cpos/$SEQ.npz"

FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )
cd $MAIN_DIR

# First, u_mag
OUT_FOLDER=$PROJ_FOLDER/viz/"$VTYPE"
MERGE_FOLDER=$PROJ_FOLDER/viz/"$VTYPE"_merge
VID_FOLDER=$PROJ_FOLDER/vids/"$VTYPE"

parallel -j11 ~/xvfb-run-safe python ./make_viz/00_viz.py {} $OUT_FOLDER $VTYPE $NSTEPS $MESHING_DIR $CPOS_FILE ::: $FOLDERS

python ./make_viz/01_merge_viz.py $OUT_FOLDER $MERGE_FOLDER $VID_FOLDER 

################################################################################
NSTEPS=2
SEQ="seq_005"
VTYPE="u_mag"

MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_analysis"
PROJ_FOLDER="/scratch/s/steinman/macdo708/neuromorph_analysis/$SEQ"
RESULTS_FOLDER="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/results"
MESHING_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/meshing"
CPOS_FILE="/scratch/s/steinman/macdo708/neuromorph_analysis/cpos/$SEQ.npz"

FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )
cd $MAIN_DIR

# First, u_mag
OUT_FOLDER=$PROJ_FOLDER/viz/"$VTYPE"
MERGE_FOLDER=$PROJ_FOLDER/viz/"$VTYPE"_merge
VID_FOLDER=$PROJ_FOLDER/vids/"$VTYPE"

parallel -j11 ~/xvfb-run-safe python ./make_viz/00_viz.py {} $OUT_FOLDER $VTYPE $NSTEPS $MESHING_DIR $CPOS_FILE ::: $FOLDERS

python ./make_viz/01_merge_viz.py $OUT_FOLDER $MERGE_FOLDER $VID_FOLDER 

################################################################################
NSTEPS=2
SEQ="seq_006"
VTYPE="u_mag"

MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_analysis"
PROJ_FOLDER="/scratch/s/steinman/macdo708/neuromorph_analysis/$SEQ"
RESULTS_FOLDER="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/results"
MESHING_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/meshing"
CPOS_FILE="/scratch/s/steinman/macdo708/neuromorph_analysis/cpos/$SEQ.npz"

FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )
cd $MAIN_DIR

# First, u_mag
OUT_FOLDER=$PROJ_FOLDER/viz/"$VTYPE"
MERGE_FOLDER=$PROJ_FOLDER/viz/"$VTYPE"_merge
VID_FOLDER=$PROJ_FOLDER/vids/"$VTYPE"

parallel -j11 ~/xvfb-run-safe python ./make_viz/00_viz.py {} $OUT_FOLDER $VTYPE $NSTEPS $MESHING_DIR $CPOS_FILE ::: $FOLDERS

python ./make_viz/01_merge_viz.py $OUT_FOLDER $MERGE_FOLDER $VID_FOLDER 