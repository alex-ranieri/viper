S=$1; # Arquivo All Stats
F=$2; # Arquivo All.fasta
H=$3;  # data
P=$4; # Number of processors to use

# Sort All Stats file for latter
sort -ru ${S} > ${S:0:-4}_sorted.tsv

# Run NextClade 

mkdir -p ${H}__nextClade # Create NextClade output folder
nextclade run  --input-dataset $PIPELINE/SARS-CoV-2/nextstrain_files/ --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S --output-all ${H}__nextClade --output-basename nextclade --jobs ${P} ${F}

# Join NextClade results with Sorted All Stats using python
python $PIPELINE/SARS-CoV-2/process_AllStats_NextClade.py ${S:0:-4}_sorted.tsv ${H}__nextClade/nextclade.tsv

# Cleaning intermediate files
rm ${S:0:-4}_sorted.tsv
