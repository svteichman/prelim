#!/bin/bash
#SBATCH --job-name combine      # Set a name for your job.
#SBATCH --partition short     # Slurm partition to use
#SBATCH --ntasks 1          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-12:00       # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=1000     # Memory limit for each tasks (in MB)
#SBATCH -o outfile_%j.out    # File to which STDOUT will be written
#SBATCH -e errfile_%j.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=teichs@uw.edu # Email to which notifications will be sent
module load R
Rscript combine_bhq.R 