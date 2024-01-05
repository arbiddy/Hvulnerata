#!/bin/bash

conda install -c bioconda stacks=2.53

#SBATCH -J stacks ##this is the output name for err and out
#SBATCH -N 1
#SBATCH --ntasks-per-node=4 
#SBATCH -o %x.%j.out ##output file
#SBATCH -e %x.%j.err ## err file
#SBATCH -p long ## partition used
#SBATCH --export=ALL

# process_radtags

mkdir samples

process_radtags -P -p ./raw/ -o ./samples/ -c -q -r --renz_1 nlaIII --renz_2 mluCI -i gzfastq -y fastq


# ustacks

mkdir ustacks

for X in $(cat sample_list)
do 
ustacks -f ./samples/${X}-01.1.fq.gz -o ./ustacks/ -i 1 -m 5 -M 4 -p 16
done

# cstacks

cstacks -P ./ustacks/ -M ./samples/popmap.txt -n 4 -p 15 

# sstacks

sstacks -P ./ustacks/ -M ./samples/popmap.txt -p 8

# tsv2bam

mkdir filtered

tsv2bam -P ./ustacks/ -M ./samples/popmap.txt -R ./filtered -t 8

# gstacks

gstacks -P ./ustacks/ -M ./samples/popmap.txt -t 8

# populations 

mkdir vulnerata

populations -P ./ustacks1/ -O ./vulnerata/ --popmap ./samples/popmap.txt -p 3 -r 0.5 -R 0.5 -f p_value -t 18 --fstats --vcf --plink –structure --genepop --write-single-snp
