# WHAM for CPMD
### A tool for obtaining free enregy profile from CPMD umbrella sampling

How it works
======
1. Openmp (omp) parallelized
2. Number of windows and directory index are read from the input file.
3. Details on each windows are read directly from the CPMD input files.
4. TRAJECTORIES are moved to debug folder for parallelization.
5. Number of skips is read form input.
6. Histograms are written in files under debug.

Todo list
======
1. Apply symmetry smartly or by input keywords.
2. Properly handle trajectory files when the program exits incorrectly.
3. Reformat the input file.
