# Use mehdi's post processing script to generate SPI, OSI, etc.

SCRIPT="/home/s/steinman/mnajafi/.conda/envs/.test/Post/hemodynamic_indices_all_spi_run_me.py"

SEQ="seq_006"

# MAIN_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/"
PROJ_DIR="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ"
RESULTS_DIR=$PROJ_DIR/results

cd $PROJ_DIR

RESULTS=$(find $RESULTS_DIR -maxdepth 1 -mindepth 1 -type d -name 'art*')

for d in $RESULTS; do 
    mkdir -p $d/logs
    python $SCRIPT $d --redo;
done