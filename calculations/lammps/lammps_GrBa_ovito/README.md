execution order:

NEED TO CHANGE

2.nitish_ovito_gr.py
	change num_frames to how many frames you want to use. will take last num_frames to calculate
	change cutoff_distance to distance where g(r) will be measured till
	change num_bins for resolution of plot (more bins = more resolution) 

descriptions
1.last_n_frames.sh
	takes i the "all_lammps.xyz" file and outputs "last_{n}.xyz" file
	"all_lammps.xyz" in lammps xyz format for dump of all timesteps. 
	this will take the last n steps and put that into a file

2.nitish_ovito_gr.py
	uses ovito to input the last_{n}.xyz file and calculates g(r) from it
