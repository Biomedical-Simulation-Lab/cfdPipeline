#!/bin/bash
#///////////////////////////////////////////////////////////////
#// CFD job script wrapper
#// Copyright (C) 2018 Mehdi Najafi (mnuoft@gmail.com)
#// Distribution of this file is not allowed in any form.
#///////////////////////////////////////////////////////////////


#casename="meshname"
#cycles=2
#timesteps_per_cycle=9600
#uOrder=1
#save_frequency=10
#estimated_required_time="24:00:00"
##partition="debug"
#partition="compute"
#post_processing_time_minutes=60

#num_cores=80


if [ $debug == "on" ]; then 
  partition="debug"
else
  partition="compute"
fi

if [[ ! -z jobname ]]; then
  jobname="${casename}_C${cycles}_TS${timesteps_per_cycle}"
fi

force="No";
clean="No";
restart_no_given=-1
for var in "$@"; do
  if [ "$var" == "force" ]; then force="Yes"; fi
  if [ "$var" == "clean" ]; then clean="Yes"; fi
  if [ $var -eq $var 2> /dev/null ]; then restart_no_given=$var; fi
done

sbatch --time=$estimated_required_time --job-name=$jobname --nodes=1 --ntasks-per-node=$num_cores --ntasks-per-node=$num_cores --partition=$partition  << EOST
#!/bin/bash
#///////////////////////////////////////////////////////////////
#// FEniCS job script
#// Copyright (C) 2018 Mehdi Najafi (mnuoft@gmail.com)
#// Distribution of this file is not allowed in any form.
#///////////////////////////////////////////////////////////////

#SBATCH --mail-type=NONE
#SBATCH --output=hpclog/art_%x_stdout_%j.txt

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

FEniCS_HOME=/home/s/steinman/mnajafi/.conda/envs/fenics201810
export DIJITSO_CACHE_DIR=\$SCRATCH/.local/dijitso-cache/
export INSTANT_CACHE_DIR=\$SCRATCH/.local/instant-cache/
export INSTANT_ERROR_DIR=\$SCRATCH/.local/instant-error/

cd \$SLURM_SUBMIT_DIR

#module purge
module load anaconda3/5.1.0 intel/2018.2 intelmpi/2018.2 boost/1.66.0 petsc/3.8.4 hdf5-mpi/1.8.20 netcdf-mpi/4.6.1 trilinos/12.12.1 fenics/2018.1.0
source activate \$FEniCS_HOME

casename_full=\$(python \$FEniCS_HOME/.test/Solver/naming.py $casename $uOrder $cycles $timesteps_per_cycle)
results_folder=./results/\$casename_full
uco=$(( $timesteps_per_cycle / $save_frequency))
#echo \$casename_full \$results_folder \$uco
simdone="0"
restart_no=0
if [ $clean == "Yes" ]; then
  echo "(!) Removing any previous results and logs."
  rm -r \$results_folder 2>/dev/null
  rm  ./logs/\${casename_full}_restart_* 2>/dev/null
fi

if [ $restart_no_given -lt 0 ]; then
  # Here we need to check if the simulation is done or not
  # find the last checkpoint folder
  #restart_no=0
  for l in \$(ls -1 \$results_folder/data/ 2>/dev/null); do if [[ \$l -gt \$restart_no ]]; then restart_no=\$l; fi; done
  uc=\$(ls -1 \$results_folder/*_up.h5 2>/dev/null | wc -l)
  #echo \$uco  " " \$uc
  if [ \$uco -eq \$uc ]; then
    simdone="1"
  fi
else
  restart_no=$restart_no_given
fi

# overwrites any previous attempts
log_file=./logs/\${casename_full}_restart_\${restart_no}

printf "~%.0s" {1..80};echo
echo "$(date)"
printf "~%.0s" {1..80};echo
echo "Aneurysm CFD v.1.2 - by Mehdi Najafi, Copyright (c) 2018 - BioMedical Simulation Lab - University of Toronto"
echo "Case Name: " $casename
echo "Number of Cycles: " $cycles
echo "Number of Time Steps per Cycle: " $timesteps_per_cycle
echo "Number of Time Steps to skip to output: " $save_frequency
echo "Case Full Name: " \$casename_full
echo "Results folder: " \$results_folder

if [ \$simdone == "0" ]; then
  if [ \$restart_no -gt 0 ]; then
    if [ -f \$results_folder/data/\$restart_no/incomplete ] || [ $force == "Yes" ]; then
      echo "Try#: " \$(( \$restart_no+1 ))
      echo "Log file: " \$log_file
      restart_folder=\${results_folder}/data/\${restart_no}/Checkpoint
      echo "restart_folder: " \$restart_folder
      echo \$(date) > \$log_file
      mpirun -n $num_cores python \$FEniCS_HOME/.test/Solver/ns uOrder=$uOrder timesteps=$timesteps_per_cycle cycles=$cycles save_frequency=$save_frequency \
             mesh_name=$casename restart_folder=\$restart_folder &>> \$log_file
    else
      echo "<!> Something went wrong! You should inspect what happened at Checkpoint#\${restart_no} in order to avoid an infinite loop at this point."
      echo "Here it is: \$restart_folder"
      echo "At the previous try the solver sould have crashed!"
      exit
    fi
  else
    echo "Try#: " \$(( \$restart_no+1 ))
    echo "Log file: " \$log_file
    echo \$(date) > \$log_file
    mpirun -n $num_cores python \$FEniCS_HOME/.test/Solver/ns uOrder=$uOrder timesteps=$timesteps_per_cycle cycles=$cycles save_frequency=$save_frequency \
           mesh_name=$casename &>> \$log_file
  fi
  echo \$(date) >> \$log_file
  echo "hpclog/art_\${SLURM_JOB_NAME}_stdout_\${SLURM_JOB_ID}.txt" >> \$log_file
  sleep 30
fi

if [ -f \$log_file ]; then
  # find the last checkpoint folder
  for l in \$(ls -1 \$results_folder/data/ 2>/dev/null); do if [[ \$l -gt \$m ]]; then m=\$l; fi; done; restart_no=\$m
  if [ -f \$results_folder/data/\$restart_no/incomplete ]; then
    echo "Submitting a job to resume previous CFD the simulation for another try: #" \$(( \$restart_no+1 ))
    ssh nia-login03 "cd \$SLURM_SUBMIT_DIR; bash $0 \$restart_no"
    exit
  fi
fi

if [ -f \$results_folder/data/\$restart_no/complete ]; then
    # simulation is done
    echo "Simulation is finished."
    echo "Sleeping for 15 seconds to let the I/O settle down."
    sleep 15
    uco=\$(( $timesteps_per_cycle / $save_frequency))
    uc=\$(ls -1 \$results_folder/*_up.h5 2>/dev/null | wc -l)
    if [ \$uco -ne \$uc ]; then
      echo "<!> No enough outputs found! Expecting " \$uco "files but found " \$uco " files!"
      echo "What to do: inspect everything first. ONLY, if there was an IO issue:"
      echo "1) Try removing file:" \$results_folder/data/\$restart_no/complete
      echo "2) Re-run this script: bash \$0"
    else
      echo "Submitting post-processing jobs for:" \$casename " stored at " \$results_folder
      wc=\$(ls -1 \$results_folder/wss_files/*_wss.h5 2>/dev/null | wc -l)
      if [ \$uco -ne \$wc ]; then
        echo "Submitting wss calculations job for:" \$casename " stored at " \$results_folder
        out=\$((ssh nia-login01 "cd \$SLURM_SUBMIT_DIR; python \$FEniCS_HOME/.test/Solver/Post/wss_run_me.py \$results_folder -t $post_processing_time_minutes")  2>&1)
        echo "\$out"
        wss_jid="\${out#*Submitted batch job}";wss_jid=\${wss_jid##* };
        # echo "<\$wss_jid>  $0 $1 \$SLURM_SUBMIT_DIR"
        if [ ! -z \$wss_jid ]; then
          #sleep 5
          out=\$((ssh nia-login01 "cd \$SLURM_SUBMIT_DIR; printf '#!/bin/bash\nssh nia-login01 \"cd \$SLURM_SUBMIT_DIR; bash $0\"' | sbatch --dependency=afterany:\$wss_jid --time=00:15:00 --mail-type=NONE --nodes=1 --ntasks-per-node=1 --job-name ${jobname}_1 --output=hpclog/art_%x_pstdout_%j.txt") 2>&1)
          echo "\$out"
        fi
      else
        if [ ! -f \$results_folder/*_hemodynamics_w.tec ] || [ "$force" == "Yes" ]; then
          echo "Submitting hemodynamic calculations job for:" $casename " stored at " \$results_folder
          out=\$((ssh nia-login01 "cd \$SLURM_SUBMIT_DIR; python \$FEniCS_HOME/.test/Solver/Post/hemodynamic_indices_run_me.py \$results_folder -t $post_processing_time_minutes") 2>&1)
          echo "\$out"
        else
          echo "Hemodynamics file were previously generated. Nothing to do!"
        fi
      fi
    fi
fi
EOST

