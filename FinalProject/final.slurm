#!/bin/bash -x
#SBATCH --job-name="Final Project"
#comment	="Building multi-resolution structure: ..."
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=64
#SBATCH --export=ALL
#SBATCH --time=01:00:00
#SBATCH -o final.out
#SBATCH -e final.err

# Go to the directory from which the job was launched
cd $SLURM_SUBMIT_DIR

# Create a short JOBID based on the one provided from the scheduler
JOBID='echo $SLURM_JOBID'

# Save a copy of our environment and script
cat $0 > script.$JOBID
printenv> env.$JOBID

# Run the job
# The echo will go into the standard output for this job
# The standard output file will end up in the directory
# from which the job was launched

echo "running job"
srun main_exe
echo "job has finished"
