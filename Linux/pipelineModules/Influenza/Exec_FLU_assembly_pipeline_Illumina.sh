#!/bin/bash
FOLDER=${PWD##*/}
THREADS=$1
SAMPLES=$2

# Create list of samples based on fastq files
L=$(ls *.fastq.gz); for f in $L; do g=${f%.*}; echo ${g%_L00*};  done | uniq > Lista.txt


# Init assemblies in folder
cat Lista.txt | xargs -P ${SAMPLES} -I {} sh -c "bash $PIPELINE/Influenza/Influenza_assembly_v3.sh {} ${FOLDER} ${THREADS}"

# Join all fasta and statistics files
#cat *.fasta > All_Fastas__${FOLDER}.fas
cat *_complete.Statistics | sort -ru > All_Statistics__${FOLDER}.tsv

# Process CeVIVAS output
assemlby_path=$(realpath .)
while read value; do read_path=$(realpath $value*R1*); echo -e $value '\t' $read_path; done < Lista.txt > read_path.tsv

python $PIPELINE/Influenza/write_flu_CeVIVAS_output_v2.py All_Statistics__${FOLDER}.tsv $FOLDER

