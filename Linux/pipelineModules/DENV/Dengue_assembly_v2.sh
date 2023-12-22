G=$1;  # nome base da amostra
H=$2;  # data
F=DENV__${G}__${H};
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

# Excluir primers e adaptadores

cutadapt -b file:$PIPELINE/DENV/primer_and_adapter_colection_dengue.fasta -B file:$PIPELINE/DENV/primer_and_adapter_colection_dengue.fasta -j ${THREADS} -o ${F}_R1_cutadapt.fq.gz -p ${F}_R2_cutadapt.fq.gz ${F}_R1.fq.gz ${F}_R2.fq.gz

# Limpar reads.
trimmomatic PE -phred33 ${F}_R1_cutadapt.fq.gz ${F}_R2_cutadapt.fq.gz ${F}_R1_paired.fq.gz  ${F}_R1_unpaired.fq.gz  ${F}_R2_paired.fq.gz  ${F}_R2_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:35 TOPHRED33 -threads ${THREADS}

rm -rf ${F}_R1_cutadapt.fq.gz ${F}_R2_cutadapt.fq.gz

## 1a etapa - mapear as reads em cada sorotipo:

sorotipo=(1 2 3 4)

for i in "${sorotipo[@]}" 
do
    bwa mem -t ${THREADS} $PIPELINE/DENV/Reference_Genomes/dengue_${i}.fasta ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz > ${F}_dengue_${i}_mapped.sam
    samtools view -S -b -q 30 ${F}_dengue_${i}_mapped.sam > ${F}_dengue_${i}_mapped.bam
    samtools sort ${F}_dengue_${i}_mapped.bam -o ${F}_dengue_${i}_mapped.sorted.bam
    samtools index ${F}_dengue_${i}_mapped.sorted.bam
    rm ${F}_dengue_${i}_mapped.sam
done

## 2a etapa - calcular a porcentagem de reads mapeadas e decidir qual sorotipo usar com base na maior quantidade de reads mapeadas
for i in "${sorotipo[@]}" 
do
    samtools view -c -F 260 ${F}_dengue_${i}_mapped.sorted.bam > ${F}_${i}.ReadsMapped;
done 

awk '{print $0"\t"FILENAME}' *.ReadsMapped > stats_ReadsMapped.tsv
selected_serotype=$(sort -k1 -rn stats_ReadsMapped.tsv | head -n1 | cut -f 2 | tail -c 14 | head -c 1)


## 3a etapa - usar pilon no sorotipo selecionado

pilon --genome $PIPELINE/DENV/Reference_Genomes/dengue_${selected_serotype}.fasta --frags ${F}_dengue_${selected_serotype}_mapped.sorted.bam --output ${F}_${selected_serotype}.Pilon --fix "gaps,indels"  --threads ${THREADS} --mindepth 5 --minmq 30;


## 4a etapa - remapear em cima do fasta pilon

bwa index ${F}_${selected_serotype}.Pilon.fasta
bwa mem -t ${THREADS} ${F}_${selected_serotype}.Pilon.fasta ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz > ${F}_dengue_${selected_serotype}_Pilon_mapped.sam
samtools view -S -b -q 30 ${F}_dengue_${selected_serotype}_Pilon_mapped.sam > ${F}_dengue_${selected_serotype}_Pilon_mapped.bam
samtools sort ${F}_dengue_${selected_serotype}_Pilon_mapped.bam -o ${F}_dengue_${selected_serotype}_Pilon_mapped.sorted.bam
samtools index ${F}_dengue_${selected_serotype}_Pilon_mapped.sorted.bam
rm ${F}_dengue_${selected_serotype}_Pilon_mapped.sam

## 5a etapa - Pegar consenso com iVar
samtools mpileup -aa -A -d 0 -Q 0 ${F}_dengue_${selected_serotype}_Pilon_mapped.sorted.bam | ivar consensus -q 20 -m 5 -p ${F}_ivar -i ${F}_pilon

## 6a etapa - Ajuste de remapeamento
bwa index ${F}_ivar.fa
bwa mem ${F}_ivar.fa  ${F}_R1_paired.fq.gz  ${F}_R2_paired.fq.gz  -t ${THREADS} | samtools view -S -b -q 30 > ${F}.Pilon2.bam
samtools sort -o ${F}.Pilon2.sorted.bam ${F}.Pilon2.bam
samtools index ${F}.Pilon2.sorted.bam

## 7a etapa -  SNP calling
samtools mpileup -aa -A -d 0 -B -Q 0 ${F}.Pilon2.sorted.bam | ivar variants -p ${F}_iVar_variants -q 20 -t 0.25 -m 5 -r ${F}_ivar.fa
grep 'TRUE' ${F}_iVar_variants.tsv | grep -vP 'N\t' | grep -vP '\tN'  | wc -l > ${F}.SNPsCount

python $PIPELINE/DENV/substitute_degenarate_bases.py ${F}_ivar.fa Genoma_${F}.fasta
rm -rf ${F}_ivar.fa

# Obtain assembly statistics.
#samtools view -c ${F}.Pilon2.sorted.bam > ${F}.ReadCount;
echo $(zcat ${F}_R1.fq.gz ${F}_R1.fq.gz | wc -l)/4|bc > ${F}.ReadCount;
samtools view -c -F 260 ${F}.Pilon2.sorted.bam > ${F}.ReadsMappedFinal;
x=$(cat ${F}.ReadsMappedFinal); y=$(cat ${F}.ReadCount); python -c "print(round(float(${x}/${y}*100), 2))" > ${F}.PercentMapped;
samtools depth -a  ${F}.Pilon2.sorted.bam  |  awk '{sum+=$3} END {print sum/NR}' > ${F}.MeanDepth;
samtools depth -a  ${F}.Pilon2.sorted.bam  |  awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}' > ${F}.MedianDepth;
samtools depth -a  ${F}.Pilon2.sorted.bam  |  awk '{print $3 >= 10}' | grep '1' | wc -l > ${F}.Depth10;
samtools depth -a  ${F}.Pilon2.sorted.bam  |  awk '{print $3 >= 25}' | grep '1' | wc -l > ${F}.Depth25;
bases=$(seqtk comp Genoma_${F}.fasta | awk '{print $1 "\t" ($3+$4+$5+$6)}' | cut -f 2)
total=$(seqtk comp Genoma_${F}.fasta | cut -f2)
echo "scale=4; ($bases / $total) * 100" | bc > ${F}.CoverageBlastn
if [[ ${selected_serotype} -eq 1 ]]
then
    #cat ${F}.blastn | awk '{x = $8-$7; print x < 0 ? -x+1 : x+1}' | awk '{sum+=$1} END {coverage = sum/10735 * 100"%"; printf "%0.2f\n", coverage}' > ${F}.CoverageBlastn;
    #check gene coverage - use sed with " to create bed in order to avoid changing fasta header after mapping
    for gene in E NS1 NS3 NS5
    do
        cat $PIPELINE/DENV/gene_check/denv1_${gene}.bed | sed "s/gene/${F}_pilon/g" > denv1_${gene}.bed
        samtools depth -a  ${F}.Pilon2.sorted.bam -b denv1_${gene}.bed | awk '{sum+=$3} END {print sum/NR}' > ${F}.Coverage_${gene}_MeanDepth
        bases=$(seqtk comp -r denv1_${gene}.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($7+$4+$5+$6)}' | cut -f 2)
        total=$(seqtk comp -r denv1_${gene}.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($3-$2)}' | cut -f2)
        echo "scale=4; ($bases / $total) * 100" | bc > ${F}.Coverage_gene_${gene}
    done
elif [[ ${selected_serotype} -eq 2 ]]
then
    for gene in E NS1 NS3 NS5
    do
        cat $PIPELINE/DENV/gene_check/denv2_${gene}.bed | sed "s/gene/${F}_pilon/g" > denv2_${gene}.bed
        samtools depth -a  ${F}.Pilon2.sorted.bam -b denv2_${gene}.bed | awk '{sum+=$3} END {print sum/NR}' > ${F}.Coverage_${gene}_MeanDepth
        bases=$(seqtk comp -r denv2_${gene}.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($7+$4+$5+$6)}' | cut -f 2)
        total=$(seqtk comp -r denv2_${gene}.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($3-$2)}' | cut -f2)
        echo "scale=4; ($bases / $total) * 100" | bc > ${F}.Coverage_gene_${gene}
    done
elif [[ ${selected_serotype} -eq 3 ]]
then
    for gene in E NS1 NS3 NS5
    do
        cat $PIPELINE/DENV/gene_check/denv3_${gene}.bed | sed "s/gene/${F}_pilon/g" > denv3_${gene}.bed
        samtools depth -a  ${F}.Pilon2.sorted.bam -b denv3_${gene}.bed | awk '{sum+=$3} END {print sum/NR}' > ${F}.Coverage_${gene}_MeanDepth
        bases=$(seqtk comp -r denv3_${gene}.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($7+$4+$5+$6)}' | cut -f 2)
        total=$(seqtk comp -r denv3_${gene}.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($3-$2)}' | cut -f2)
        echo "scale=4; ($bases / $total) * 100" | bc > ${F}.Coverage_gene_${gene}
    done
else
    for gene in E NS1 NS3 NS5
    do
        cat $PIPELINE/DENV/gene_check/denv4_${gene}.bed | sed "s/gene/${F}_pilon/g" > denv4_${gene}.bed
        samtools depth -a  ${F}.Pilon2.sorted.bam -b denv4_${gene}.bed | awk '{sum+=$3} END {print sum/NR}' > ${F}.Coverage_${gene}_MeanDepth
        bases=$(seqtk comp -r denv4_${gene}.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($7+$4+$5+$6)}' | cut -f 2)
        total=$(seqtk comp -r denv4_${gene}.bed Genoma_${F}.fasta | awk '{print $1 "\t" ($3-$2)}' | cut -f2)
        echo "scale=4; ($bases / $total) * 100" | bc > ${F}.Coverage_gene_${gene}
    done
fi
seqtk comp Genoma_${F}.fasta | awk '{x+=$9}END{print x}' > ${F}.CountNs;
echo $selected_serotype > ${F}.SeroType 

nextclade run -D $PIPELINE/DENV/nextclade_files/denv${selected_serotype} -j 3 -t nextclade.tsv Genoma_${F}.fasta
csvcut -t -c clade nextclade.tsv | tail -n+2 > ${F}.Genotype

# Print statistics to .Statistics file.
ls Genoma_${F}.fasta > ${F}.GenomeName;
printf "Genome\tN_Reads\tReads_mapped\tPercent_mapped\tMean_depth\tMedian_depth\tNpos_Depth>=10\tNpos_Depth>=25\tCoverage\tE_Coverage\tE_Depth\tNS1_Coverage\tNS1_Depth\tNS3_Coverage\tNS3_Depth\tNS5_Coverage\tNS5_Depth\tNumber_of_Ns\tSNPs\tSerotype\tGenotype\n" > ${F}.Statistics;
paste -d "\t" *GenomeNam* *ReadC* *ReadsMappedFinal *Percent* *.MeanDepth *.MedianDepth *Depth10 *Depth25 *CoverageBlastn *.Coverage_gene_E *.Coverage_E_MeanDepth *.Coverage_gene_NS1 *.Coverage_NS1_MeanDepth *.Coverage_gene_NS3 *.Coverage_NS3_MeanDepth *.Coverage_gene_NS5 *.Coverage_NS5_MeanDepth *CountN* *.SNPsCount *.SeroType *.Genotype >> ${F}.Statistics;


# Copy Genome and Statistics to folder up.
cp Genoma_${F}.fasta ../;
cp ${F}.Statistics ../;

rm *.gz
rm  *_Pilon_mapped.sorted.bam ${F}.Pilon2.bam

# Back to folder up.
cd ..;
