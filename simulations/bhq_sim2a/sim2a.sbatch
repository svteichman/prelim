#!/bin/bash
#SBATCH --job-name sim2a      # Set a name for your job.
#SBATCH --partition medium     # Slurm partition to use
#SBATCH --ntasks 1          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-12:00       # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=2000     # Memory limit for each tasks (in MB)
#SBATCH -o outfile_%j.out    # File to which STDOUT will be written
#SBATCH -e errfile_%j.err    # File to which STDERR will be written
#SBATCH --mail-type=FAIL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=teichs@uw.edu # Email to which notifications will be sent
module load R
Rscript BHq_knockoffs_sim2a.R $SLURM_ARRAY_TASK_ID