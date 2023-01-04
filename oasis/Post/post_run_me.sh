#!/bin/bash
# ///////////////////////////////////////////////////////////////
#  Job to Post-Processing
#  Author: Mehdi Najafi, mnuoft@gmail.com
#  Date: 2017-07-01
#
#  This library is not intended for distributions in any form and
#  distribution of it is not allowed in any form.
# ///////////////////////////////////////////////////////////////

#SBATCH --partition=debug
#SBATCH --time=00:15:00
#SBATCH --mail-type=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name postproc
#SBATCH --output=hpclog/postproc_stdout_%j.txt
#SBATCH --error=hpclog/postproc_stderr_%j.txt

cd $SLURM_SUBMIT_DIR

sleep 15

Solver=~/.conda/envs/.test/

for var in "$@"
do
   echo "Sending post-processing jobs for:" $var
   out=$(ssh nia-login03 "cd $SLURM_SUBMIT_DIR; python $Solver/Post/volumetric_indices_run_me.py $var");
   out=$(ssh nia-login03 "cd $SLURM_SUBMIT_DIR; python $Solver/Post/wss_run_me.py $var");
   echo $out
   wss_jid=${out#*Submitted batch job};

   echo $wss_jid
   if [ ! -z $wss_jid ]
   then
     sleep 5
     out=$(ssh nia-login03 "cd $SLURM_SUBMIT_DIR; python $Solver/Post/hemodynamic_indices_run_me.py $var -p $wss_jid");
     echo $out
   elif [[ $out == *"No files left to compute WSS."* ]]
   then
     out=$(ssh nia-login03 "cd $SLURM_SUBMIT_DIR; python $Solver/Post/hemodynamic_indices_run_me.py $var");
   fi
done
