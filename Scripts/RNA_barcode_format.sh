#!/bin/bash
#SBATCH --mail-user=user@mail.org
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err

grep 'CGAAAGAATGCAGATTTTCGTCAAGAC' sample_R1_001.fastq > sample_RNA.txt

sed 's/^/>\n/' sample_RNA.txt > sample_RNA_form.fa