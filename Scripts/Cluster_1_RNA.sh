#!/bin/bash
#SBATCH --mail-user=user@mail.org
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err

module load VSEARCH

vsearch --cluster_size $1 --consout - --id 0.99 --minseqlength 90
