RESULTS_FOLDER="/scratch/s/steinman/macdo708/surge_cfd/results"
FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )

for i in $FOLDERS; do 
    echo $(basename $i); ls $i/wss_files/* | wc -l;
done