#!/bin/bash
#SBATCH --mail-user=user@mail.org
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err

module load VSEARCH

cat list of slurm out files > sample_clust_cat.fa

awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' sample_clust_cat.fa > sample_clust_cat_linear.fa

sed 's/seqs/size/g' sample_clust_cat_linear.fa > sample_clust_form.fa

vsearch --cluster_size sample_clust_form.fa --sizein --sizeout --consout sample_final.clust.fa --centroids sample_final.cent.fa --id 0.99 --minseqlength 90

awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' sample_final.clust.fa > sample_final_lin.fa

cat sample_final_lin.fa | paste - - > sample_clust_table.txt