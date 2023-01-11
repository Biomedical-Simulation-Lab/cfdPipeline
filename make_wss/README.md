# Making WSS

Originally, just wanted to inspect/debug Mehdi's WSS code. Had to copy a few files over and change a few things to get it working.
- In line 48 of job_utils, had to remove a `/` before `data/`. 
- For line 85 and 86, added a check for file exists: `if os.path.exists(mesh_filename[:-3]+'_normals.mat'):`

Location of original files:
`~/../mnajafi/.conda/envs/.test/Post/wss_run_me.py`
`~/../mnajafi/.conda/envs/.test/Post/hemodynamic_indices_run_me.py`

So now, run this `wss_run_me.py` file with your results folder.