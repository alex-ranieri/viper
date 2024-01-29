G=$1;  # nome base da amostra
H=$2;  # data
F=Indels__${G}__${H};
THREADS=$3

# Go to this folder.
s=$(pwd); cd $s


# Create folder (if it doesn't exist) for each genome to be processed.
mkdir -p ${F};


# Unzip sequences.
cat ${G}*R1*.fastq.gz  > ${F}/${F}_R1.fq.gz
cat ${G}*R2*.fastq.gz  > ${F}/${F}_R2.fq.gz



# Enter Analysis Folder.
cd ${F}



# Limpar reads.
trimmomatic PE -phred33 ${F}_R1.fq.gz ${F}_R2.fq.gz  ${F}_R1_paired.fq.gz  ${F}_R1_unpaired.fq.gz  ${F}_R2_paired.fq.gz  ${F}_R2_unpaired.fq.gz  ILLUMINACLIP:$PIPELINE/SARS-CoV-2/primer_and_adapter_colection.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36  -threads ${THREADS}

# Juntar reads unpaired

cat ${F}_R1_unpaired.fq.gz ${F}_R2_unpaired.fq.gz > ${F}_unpaired.fq.gz
rm -rf ${F}_R1_unpaired.fq.gz ${F}_R2_unpaired.fq.gz

## 1a etapa:

bowtie2 -p ${THREADS} -x $PIPELINE/SARS-CoV-2/Reference_Genomes/Ref_Wuhan -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz -U ${F}_unpaired.fq.gz | samtools view -S -b -q 10 > ${F}_Wuhan_mapped.bam
samtools sort ${F}_Wuhan_mapped.bam -o ${F}_Wuhan_mapped.sorted.bam
samtools index ${F}_Wuhan_mapped.sorted.bam


## 2a etapa:

pilon --genome $PIPELINE/SARS-CoV-2/Reference_Genomes/Ref_Wuhan.fasta --frags ${F}_Wuhan_mapped.sorted.bam --minqual 20 --minmq 10 --output ${F}.Pilon --fix "gaps,indels"  --threads ${THREADS} --mindepth 10;

## 3a etapa:

bowtie2-build ${F}.Pilon.fasta ${F}.Pilon.fasta
bowtie2 -p ${THREADS} -x ${F}.Pilon.fasta -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz -U ${F}_unpaired.fq.gz | samtools view -S -b -q 10 > ${F}.Pilon.bam
samtools sort ${F}.Pilon.bam -o ${F}.Pilon.sorted.bam
samtools index ${F}.Pilon.sorted.bam

## 4a etapa:

samtools mpileup -aa -A -d 0 -Q 0 ${F}.Pilon.sorted.bam | ivar consensus -q 10 -p ${F}_ivar -i NC_045512.2_pilon

python $PIPELINE/SARS-CoV-2/substitute_degenarate_bases.py ${F}_ivar.fa ${F}_ivar_clean.fa
rm ${F}_ivar.fa

## Ajuste de remapeamento
bowtie2-build ${F}_ivar_clean.fa ${F}_ivar_clean.fa
bowtie2 -p ${THREADS} -x ${F}_ivar_clean.fa  -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz -U ${F}_unpaired.fq.gz | samtools view -S -b -q 10 > ${F}.Pilon2.bam
samtools sort ${F}.Pilon2.bam -o ${F}.Pilon2.sorted.bam
samtools index ${F}.Pilon2.sorted.bam

# SNP calling

samtools mpileup -aa -A -d 0 -B -Q 0 ${F}.Pilon2.sorted.bam | ivar variants -p ${F}_iVar_variants -q 20 -t 0.03 -m 100 -r ${F}_ivar_clean.fa -g $PIPELINE/SARS-CoV-2/genemap.gff
grep 'TRUE' ${F}_iVar_variants.tsv | cut -f 1-4 | sort | uniq | grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.SNPsCount

grep 'TRUE' ${F}_iVar_variants.tsv | grep 'ORF1a' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN' | wc -l > ${F}.ORF1aSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'ORF1b' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.ORF1bSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'ORF3a' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.ORF3aSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'Spike' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN' | wc -l > ${F}.SpikeSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'Envelope' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.EnvelopeSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'Membrane_Glycoprotein' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.Membrane_GlycoproteinSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'ORF6' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.ORF6SNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'ORF7a' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.ORF7aSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'ORF7b' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.ORF7bSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'Nucleocapsid' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.NucleocapsidSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'ORF8' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.ORF8SNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep 'ORF9b' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.ORF9bSNPsCount
grep 'TRUE' ${F}_iVar_variants.tsv | grep -P '\tNA\t' | cut -f 1-4 |grep -v '+' | grep -v '-' | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.nonCodingSNPsCount


mv ${F}_ivar_clean.fa Genoma_${F}.fasta


# CHECK KEY PANGOLIN REGIONS
samtools depth -b $PIPELINE/SARS-CoV-2/v3_check_S.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.S_MeanDepth;
samtools depth -b $PIPELINE/SARS-CoV-2/v3_check_ORF1ab.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.ORF1ab_MeanDepth;
samtools depth -b $PIPELINE/SARS-CoV-2/v3_check_N.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.N_MeanDepth;

# Check Spike coverage
total=$(seqtk comp -r $PIPELINE/SARS-CoV-2/v3_check_S.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($3-$2)}' | cut -f2)
bases=$(seqtk comp -r $PIPELINE/SARS-CoV-2/v3_check_S.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($7+$4+$5+$6)}' | cut -f2)
echo "scale=4; ($bases / $total) * 100" | bc > ${F}.S_Coverage;

# Check coverage for non-overlaping regions of amplicons from different PCR
samtools depth -b $PIPELINE/SARS-CoV-2/pcr_mix1.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.PCR1_MeanDepth;
samtools depth -b $PIPELINE/SARS-CoV-2/pcr_mix2.bed ${F}.Pilon2.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > ${F}.PCR2_MeanDepth;
samtools depth -b $PIPELINE/SARS-CoV-2/pcr_mix1.bed ${F}.Pilon2.sorted.bam | awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}' > ${F}.PCR1_MedianDepth;
samtools depth -b $PIPELINE/SARS-CoV-2/pcr_mix2.bed ${F}.Pilon2.sorted.bam | awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}' > ${F}.PCR2_MedianDepth;

# PANGOLIN
pangolin -t ${THREADS} Genoma_${F}.fasta --outfile Pangolin__${F}.tsv
sed -i "s/,/\t/g" Pangolin__${F}.tsv


# Obtain assembly statistics.
#samtools view -c ${F}_hmpox_mapped.sorted.bam > ${F}.ReadCount;
echo $(zcat ${F}_R1.fq.gz ${F}_R2.fq.gz | wc -l)/4|bc > ${F}.ReadCount;
samtools view -c -F 260 ${F}.Pilon2.sorted.bam > ${F}.ReadsMapped;
x=$(cat ${F}.ReadsMapped); y=$(cat ${F}.ReadCount); python -c "print(round(float(${x}/${y}*100), 2))" > ${F}.PercentMapped;
samtools depth -a  ${F}.Pilon2.sorted.bam  |  awk '{sum+=$3} END {print sum/NR}' > ${F}.MeanDepth;
samtools depth -a  ${F}.Pilon2.sorted.bam  |  awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}' > ${F}.MedianDepth;
samtools depth -a  ${F}.Pilon2.sorted.bam  |  awk '{print $3 >= 10}' | grep '1' | wc -l > ${F}.Depth10;
samtools depth -a  ${F}.Pilon2.sorted.bam  |  awk '{print $3 >= 25}' | grep '1' | wc -l > ${F}.Depth25;
blastn -query Genoma_${F}.fasta -db $PIPELINE/SARS-CoV-2/Reference_Genomes/Ref_Wuhan.fasta -outfmt 6 > ${F}.blastn;
cat ${F}.blastn | awk '{x = $8-$7; print x < 0 ? -x+1 : x+1}' | awk '{sum+=$1} END {coverage = sum/29903 * 100"%"; printf "%0.2f\n", coverage}' > ${F}.CoverageBlastn;
seqtk comp Genoma_${F}.fasta | awk '{x+=$9}END{print x}' > ${F}.CountNs;
cat Pangolin__${F}.tsv | awk NR==2 | awk '{print $2}' > ${F}.Pangolin;
cut -f 5  Pangolin__${F}.tsv | tail -n+2 > ${F}.scorpio;

# Print statistics to .Statistics file.
ls Genoma_${F}.fasta > ${F}.GenomeName;

# Rename genome fasta header.
K=Genoma_${F}.fasta # New header
sed -i "1s/^.*$/>${K}/" Genoma_${F}.fasta;


printf "Genome\tN_Reads\tReads_mapped\tPercent_mapped\tMean_depth\tMedian_depth\tNpos_Depth>=10\tNpos_Depth>=25\tCoverage\tNumber_of_Ns\tPangolin_lineage\tscorpio_call\tS_coverage\tS_Mean_depth\tORF1ab_mean_depth\tN_mean_depth\tSNPs_count\tORF1aSNPsCount\tORF1bSNPsCount\tORF3aSNPsCount\tSpikeSNPsCount\tEnvelopeSNPsCount\tMembrane_GlycoproteinSNPsCount\tORF6SNPsCount\tORF7aSNPsCount\tORF7bSNPsCount\tNucleocapsidSNPsCount\tORF8SNPsCount\tORF9bSNPsCount\tnonCodingSNPsCount\tPCR1_mean_depth\tPCR2_mean_depth\tPCR1_median_depth\tPCR2_median_depth\n" > ${F}.Statistics;
paste -d "\t" *GenomeNam* *ReadC* *ReadsM* *Percent* *.MeanDepth *.MedianDepth *Depth10 *Depth25 *CoverageBlastn *CountN* *.Pangolin *.scorpio *.S_Coverage *.S_MeanDepth *.ORF1ab_MeanDepth *.N_MeanDepth *.SNPsCount *.ORF1aSNPsCount *.ORF1bSNPsCount *.ORF3aSNPsCount *.SpikeSNPsCount *.EnvelopeSNPsCount *.Membrane_GlycoproteinSNPsCount *.ORF6SNPsCount *.ORF7aSNPsCount *.ORF7bSNPsCount *.NucleocapsidSNPsCount *.ORF8SNPsCount *.ORF9bSNPsCount *.nonCodingSNPsCount *.PCR1_MeanDepth *.PCR2_MeanDepth *.PCR1_MedianDepth *.PCR2_MedianDepth >> ${F}.Statistics;

# Copy Genome and Statistics to folder up.
cp Genoma_${F}.fasta ../;
cp ${F}.Statistics ../;

rm *.gz
rm  ${F}.Pilon.bam ${F}.Pilon2.bam ${F}.Pilon.sorted.bam ${F}_Wuhan_mapped.sorted.bam ${F}_Wuhan_mapped.bam

# Back to folder up.
cd ..;

