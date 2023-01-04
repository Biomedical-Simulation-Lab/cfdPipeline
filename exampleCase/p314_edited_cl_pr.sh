#!/bin/bash
# CFD job script wrapper
casename="p314_edited_cl_pr"
cycles=2
timesteps_per_cycle=9600
uOrder=1
save_frequency=10
estimated_required_time="24:00:00"
post_processing_time_minutes=180
num_cores=80
~/../mnajafi/solver.bin $0 "$@"
