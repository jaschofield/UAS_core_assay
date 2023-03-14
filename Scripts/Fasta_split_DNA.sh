#!/bin/bash
#SBATCH --mail-user=user@mail.org
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err

sed 's/"//g' sample_DNA.txt > sample_DNAa.txt
sed 's/[0-9]*//g' sample_DNAa.txt > sample_DNAb.txt
sed 's/ //g' sample_DNAb.txt > sample_DNAc.txt
sed '1d' sample_DNAc.txt > sample_DNAd.txt
sed 's/^/>\n/' sample_DNAd.txt > sample_form.fa

rm sample_DNAa.txt
rm sample_DNAb.txt
rm sample_DNAc.txt
rm sample_DNAd.txt


split -l 5000000 sample_form.fa sample_split_

