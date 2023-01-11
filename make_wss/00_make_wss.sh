cd "/scratch/s/steinman/macdo708/neuromorph_analysis/make_wss/submit"

SEQ="seq_005"
RESULTS_FOLDER="/scratch/s/steinman/macdo708/neuromorph_meshing/$SEQ/results"
FOLDERS=$(find "$RESULTS_FOLDER" -mindepth 1 -maxdepth 1 -type d )

for i in $FOLDERS; do 
    python ../mehdi_wss_tools/wss_run_me.py $i -c; 
done


# To fix!

# COMPUTING

# Done fix
# art_twh1323_1739_mesh_surge_smooth_cl_pr_I1_FC_MCA_10_Q746_Per951_Newt370_ts9600_cy2_uO1
# art_twh1365_2497_mesh_standard_smooth_cl_pr_I3_FC_MCA_10_Q564_Per951_Newt370_ts9600_cy2_uO1
# art_twh1043_5139_mesh_standard_smooth_cl_pr_I5_FC_MCA_10_Q796_Per951_Newt370_ts9600_cy2_uO1
# art_twh1411_3348_mesh_surge_smooth_cl_pr_I3_FC_MCA_10_Q481_Per951_Newt370_ts9600_cy2_uO1
