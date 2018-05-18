# WHAM for CPMD
### A tool for obtaining free enregy profile from CPMD umbrella sampling

How it works
======
1. Openmp (omp) parallelized
2. Number of windows and directory index are read from the input file.
3. Details on each windows are read directly from the CPMD input files.
4. TRAJECTORIES are moved to debug folder for parallelization.
5. Number of skips is read form input.
6. Symmetry of the system is read from the input file.
7. Histograms are written in files under debug.
8. Trajectory files will be moved back when the program faces error termination.

Todo list
======
1. Reformat the input file.
2. Split the global module to multiple modules
3. Determine if nsteps array is still useful or not.
