###Project overview 
The goal of this project is to assemble 6 algal genomes using short read fastq sequences, use Prodigal to predict proteins, use HMM search to search for specific giant viral markers (Linux portion), then make a graph using some of the outputs from HMM search (R Portion). 

### Environmental setup 
interact -A introtogds -p normal_q -t 1:00:00

#FastQ reads from SeqCoast, paired end reads
7637_001_S111_R1_001.fastq.gz
7637_001_S111_R2_001.fastq.gz
7637_002_S112_R1_001.fastq.gz
7637_002_S112_R2_001.fastq.gz
7637_003_S112_R1_001.fastq.gz
7637_003_S112_R2_001.fastq.gz
7637_004_S112_R1_001.fastq.gz
7637_004_S112_R2_001.fastq.gz
7637_005_S112_R1_001.fastq.gz
7637_005_S112_R2_001.fastq.gz
7637_006_S112_R1_001.fastq.gz
7637_006_S112_R2_001.fastq.gz

#Run SPAdes to assemble genome, run spades.sh

#!/bin/bash

#SBATCH --mail-user=rcesur94@vt.edu
#SBATCH --mail-type=all
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=6G
#SBATCH -t 60:00:00
#SBATCH -A aylwardlab
#SBATCH --partition=normal_q
module load gcc/8.1.0

/projects/Aylward_Lab/Robin/Sequencing_Data/7637_Illumina/
spades.py -1 7637_001_S111_R1_001.fastq.gz -2 7637_001_S111_R2_001.fastq.gz -t 16 -o 001_assembly 
spades.py -1 7637_002_S112_R1_001.fastq.gz -2 7637_002_S112_R2_001.fastq.gz -t 16 -o 002_assembly 
spades.py -1 7637_003_S113_R1_001.fastq.gz -2 7637_003_S113_R2_001.fastq.gz -t 16 -o 003_assembly 
spades.py -1 7637_004_S114_R1_001.fastq.gz -2 7637_004_S114_R2_001.fastq.gz -t 16 -o 004_assembly 
spades.py -1 7637_005_S115_R1_001.fastq.gz -2 7637_005_S115_R2_001.fastq.gz -t 16 -o 005_assembly 
spades.py -1 7637_006_S116_R1_001.fastq.gz -2 7637_006_S116_R2_001.fastq.gz -t 16 -o 006_assembly

#Run Prodigal to predict proteins 
This uses the nucleic acid fasta file, converts to proteins, then back to nucleic acid. 

prodigal -i scaffolds.fasta -a 001scaffolds.faa -d 001scaffolds.fna
prodigal -i scaffolds.fasta -a 002scaffolds.faa -d 002scaffolds.fna
prodigal -i scaffolds.fasta -a 003scaffolds.faa -d 003scaffolds.fna
prodigal -i scaffolds.fasta -a 004scaffolds.faa -d 004scaffolds.fna
prodigal -i scaffolds.fasta -a 005scaffolds.faa -d 005scaffolds.fna
prodigal -i scaffolds.fasta -a 006scaffolds.faa -d 006scaffolds.fna

#HMMER Search 
This searches for sequence homologs and makes sequence alignments using probabilistic models called profile hidden Markov models. I searched using a database we use in the lab, NCLDV_markers.hmm

hmmsearch --tblout 001hmm.hmmout -E 1e-5 NCLDV_markers.hmm 001scaffolds.faa
hmmsearch --tblout 002hmm.hmmout -E 1e-5 NCLDV_markers.hmm 002scaffolds.faa
hmmsearch --tblout 003hmm.hmmout -E 1e-5 NCLDV_markers.hmm 003scaffolds.faa
hmmsearch --tblout 004hmm.hmmout -E 1e-5 NCLDV_markers.hmm 004scaffolds.faa
hmmsearch --tblout 005hmm.hmmout -E 1e-5 NCLDV_markers.hmm 005scaffolds.faa
hmmsearch --tblout 006hmm.hmmout -E 1e-5 NCLDV_markers.hmm 006scaffolds.faa

# You can view the hmm.out files to see what markers it found. The sample with the most giant viral marker genesn



### R script 
#load required libraries
library(ggplot2)

#Read the table into R 
df <- read.table("marker_hits.tsv",
                 header = TRUE,
                 sep = "\t",
stringsAsFactors = FALSE)

#check table to see this is what you think it should look like 
head(df)

#extract contig length 
df$Contig_length <- as.numeric(
    sub(".*length_([0-9]+)_.*", "\\1", df$Contig)
 )

#Convert E-Value to base log 10
df$log10_Evalue <- -log10(as.numeric(df$E_value))

#Plot graph
ggplot(df, aes(x = Contig_length, y = log10_Evalue, color = Hit)) +
   geom_point(size = 2, alpha = 0.8) +
   theme_minimal() +
   labs(
         x = "Contig length (bp)",
        y = expression(-log[10](E-value)),
title = "HMM Marker vs Viral Contig Length"
    )

<img width="1103" height="373" alt="image" src="https://github.com/user-attachments/assets/a6d83475-16f1-4a06-a748-16ca97ffc9c6" />

