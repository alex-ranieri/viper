#!/bin/bash
FOLDER=${PWD##*/}
THREADS=$1
SAMPLES=$2

# Update pangolin
pangolin --update-data

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt

# Init assemblies in folder
cat Lista.txt | xargs -P ${SAMPLES} -I {} sh -c "bash $PIPELINE/SARS-CoV-2/Pipeline_Illumina_v7_bowtie2_ref_iVar_SNP.sh {} ${FOLDER} ${THREADS}"

# Join all fasta and statistics files
cat *.fasta > All_Fastas__${FOLDER}.fas
cat *.Statistics > All_Statistics__${FOLDER}.tsv 

# Run nextClade 
bash $PIPELINE/SARS-CoV-2/nextClade_COV.sh All_Statistics__${FOLDER}.tsv All_Fastas__${FOLDER}.fas ${FOLDER} ${THREADS}
