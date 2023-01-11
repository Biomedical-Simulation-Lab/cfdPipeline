module load CCEnv StdEnv/2020 gcc/9.3.0 vtk/9.0.1 python/3.7.7
module load openmpi/4.0.3

source ~/.virtualenvs/vtk9/bin/activate

SEQ="seq_006"

MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_analysis"
PROJ_FOLDER="/scratch/s/steinman/macdo708/neuromorph_analysis/$SEQ"
RESULTS_FOLDER="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/results"
MESHING_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/meshing"

cd $MAIN_DIR

# SURF_OUT_FOLDER=$PROJ_FOLDER/data/neck_surfs
CSV_OUT_FOLDER=$PROJ_FOLDER/data/ici_data

python ./make_hemo_features/01_merge_hemo_features.py $CSV_OUT_FOLDER
