This is my Final Project for Math 540, built for the Colorado School of Mines HPC MIO. This code implments the hierarchical domain partitioning of the multi-resolution approximation (MRA) for spatial data sets. Following are instructions to run the code. Before following the instructions, please log on to MIO and copy the files into a directory.

1) You will need to load the the required modules by typing the following commands:
module load PrgEnv/intel/latest
module load PrgEnv/mpi/intel/intel/5.0.3

2) With the modules loaded, you can type "make" to run the make file to compile the code. Doing so will produce an executable, main_exe. An outside library, Eigen, is required for compilation. This directory in included.

3) Within, userInput.txt, you may change model setup parameters should you so choose. The default settings are to run the small satellite data (n \approx 100,000). This data set is small enough that you can test it directly on the MIO log-in nodes without issue. If you test any of the other data sets, with other parameters, you will need to change the model parameters in userInput.txt and submit a batch job via the SLURM scheduler. Please see my report for more comprehensive instructions on how the userInput.txt file is structured and restrictions of the number of cores used, and how they relate to the number of levels at which regions are split across workers. In summary, make sure that the number of cores you request are a power of J, implemented to be either 2 or 4 (but typically chosen to be 2). Moreover, you may not request more cores than the number of regions at which regions are split across workers. The last parameter that you can set within userInput.txt is the level one less than the level at which regions are split across workers.

4) The model output are progress indicators and timing results.

5) If you wish you explore the codebase itself, these are the files that may be of interest: main.cpp, functions.cpp, functions.h, and buildStructure.cpp.
main.cpp: main program of the code that call the major model procedures
functions.cpp: file where many of the smaller model functions are implemented
functions.h: header file for functions.cpp and buildStructure.cpp
buildStructure.cpp: function responsible for building the multi-resolution structure
