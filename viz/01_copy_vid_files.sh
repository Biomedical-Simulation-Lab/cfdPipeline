# Copy vid files from viz to vids

viz="/scratch/s/steinman/macdo708/pereira_rupture_study_analysis/viz"
vids="/scratch/s/steinman/macdo708/pereira_rupture_study_analysis/vids"

mkdir -p $vids

for dtype in "u_mag" "qcriterion_nd"; do
    rsync -avz $viz/$dtype/*/*.mp4 $vids/$dtype/
done