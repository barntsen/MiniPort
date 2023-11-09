#!/bin/bash -l
#SBATCH --job-name=gpuport   # Job name
#SBATCH --output=log.txt     # Name of stdout output file
#SBATCH --error=log.txt      # Name of stderr error file
#SBATCH --partition=small-g  # Partition (queue) name
#SBATCH --nodes=1            # Total number of nodes 
#SBATCH --ntasks-per-node=1  # 1 MPI ranks per node, 1 total 
#SBATCH --gpus-per-node=1    # Allocate one gpu per MPI rank
#SBATCH --time=00:30:00      # Run time (d-hh:mm:ss)
###SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --account=project_465000260  # Project for billing
###SBATCH --mail-user=username@domain.com

./run.sh
