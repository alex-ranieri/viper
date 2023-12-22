G=$1;  # nome base da amostra
H=$2;  # data
F=FLU__${G}__${H};
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
echo $F > sample_name.txt

# Excluir primers e adaptadores

cutadapt -b file:$PIPELINE/Influenza/primer_and_adapter_colection.fasta -B file:$PIPELINE/Influenza/primer_and_adapter_colection.fasta -j ${THREADS} -o ${F}_R1_cutadapt.fq.gz -p ${F}_R2_cutadapt.fq.gz ${F}_R1.fq.gz ${F}_R2.fq.gz

# Limpar reads.
trimmomatic PE -phred33 ${F}_R1_cutadapt.fq.gz ${F}_R2_cutadapt.fq.gz ${F}_R1_paired.fq.gz  ${F}_R1_unpaired.fq.gz  ${F}_R2_paired.fq.gz  ${F}_R2_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:35 TOPHRED33 -threads ${THREADS}

rm -rf ${F}_R1_cutadapt.fq.gz ${F}_R2_cutadapt.fq.gz

# Join unpaired reads

cat ${F}_R1_unpaired.fq.gz ${F}_R2_unpaired.fq.gz > ${F}_unpaired.fq.gz
rm -rf ${F}_R1_unpaired.fq.gz ${F}_R2_unpaired.fq.gz


# Identify with vapor

# Identify best reference with vapor for each segment
cat $PIPELINE/Influenza/segment_list.txt | xargs -P ${THREADS} -I {} sh -c 'N=$(cat sample_name.txt); python "$PIPELINE"/Influenza/software/vapor/vapor.py -fa $PIPELINE/Influenza/segments_database/{}.fasta -fq "$N"_R1_paired.fq.gz "$N"_R2_paired.fq.gz | cut -f 6 | sed "s/>//g"> {}_result.txt'

# Get influenza identification

if [ -s segment_7_result.txt ]; then
    cut -f 1 -d '_' segment_7_result.txt > A_or_B_genotype.txt
else
    echo "not_detected" > A_or_B_genotype.txt
fi

if [ -s segment_4_result.txt ]; then
    cut -f 1 -d '_' segment_4_result.txt > H_genotype.txt
else
    echo "not_detected" > H_genotype.txt
    echo "not_detected" > clade.txt
fi

if [ -s segment_6_result.txt ]; then
    cut -f 1 -d '_' segment_6_result.txt > N_genotype.txt
else
    echo "not_detected" > N_genotype.txt
fi


# Get ref from DB

cat $PIPELINE/Influenza/segment_list.txt | xargs -P ${THREADS} -I {} sh -c "seqtk subseq $PIPELINE/Influenza/segments_database/{}.fasta {}_result.txt > {}_ref.fasta"

# Map and assemble
mkdir segments
for i in 1 2 3 4 5 6 7 8; do
    if [ -s segment_${i}_result.txt ]; then # Check if segment was detected by vapor
        bowtie2-build segment_${i}_ref.fasta segment_${i}_ref.fasta
        bowtie2 --very-sensitive -p ${THREADS} -x segment_${i}_ref.fasta -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz | samtools view -S -b > segment_${i}_mapped.bam
        samtools sort  segment_${i}_mapped.bam -o segment_${i}_mapped.sorted.bam
        rm -rf segment_${i}_mapped.bam
        samtools index segment_${i}_mapped.sorted.bam
        # Recuperando fastq
        samtools view  -u -f 1 -F 12 segment_${i}_mapped.sorted.bam > segment_${i}_mapped.paired.bam
        rm -rf segment_${i}_mapped.sorted.bam
        samtools sort  -n -o segment_${i}_mapped.sorted.paired.bam segment_${i}_mapped.paired.bam
        rm -rf segment_${i}_mapped.paired.bam
        bedtools bamtofastq -i segment_${i}_mapped.sorted.paired.bam -fq segment_${i}_R1.fq -fq2 segment_${i}_R2.fq
        gzip segment_${i}_R1.fq segment_${i}_R2.fq
        spades.py -1 segment_${i}_R1.fq.gz -2 segment_${i}_R2.fq.gz -o SPAdes_segment${i} -t ${THREADS} --only-assembler --careful --cov-cutoff 10.0
        if test -f SPAdes_segment${i}/scaffolds.fasta; then # Test if there were enough reads for assembly
            minimap2 -a segment_${i}_ref.fasta SPAdes_segment${i}/scaffolds.fasta | samtools view -Sb > minimap_segment_${i}.bam
            samtools sort -o minimap_segment_${i}.sorted.bam minimap_segment_${i}.bam
            samtools index minimap_segment_${i}.sorted.bam
            rm -rf minimap_segment_${i}.bam
            samtools mpileup -aa -A -d 1000 -C 50 -Q 0 minimap_segment_${i}.sorted.bam | ivar consensus -p minimap_segment_${i}_preconsensus -m 1 -q 0 -t 0
            python $PIPELINE/Influenza/substitute_degenarate_bases.py minimap_segment_${i}_preconsensus.fa minimap_segment_${i}_preconsensus_fixed.fa
            rm -rf minimap_segment_${i}_preconsensus.fa
            K=segment${i}_${F}.fasta # New header
            mv minimap_segment_${i}_preconsensus_fixed.fa segments/segment_preconsensus_${i}_${F}.fasta
            bowtie2-build segments/segment_preconsensus_${i}_${F}.fasta segments/segment_preconsensus_${i}_${F}.fasta
            bowtie2 --very-sensitive -p ${THREADS} -x segments/segment_preconsensus_${i}_${F}.fasta -1 segment_${i}_R1.fq.gz -2 segment_${i}_R2.fq.gz | samtools view -S -b > segments/segment_preconsensus_${i}_${F}.bam
            samtools sort  segments/segment_preconsensus_${i}_${F}.bam -o segments/segment_preconsensus_${i}_${F}.sorted.bam
            rm -rf segments/segment_preconsensus_${i}_${F}.bam segments/segment_preconsensus_${i}_${F}.fasta
            samtools mpileup -aa -A -d 0 -Q 0 segments/segment_preconsensus_${i}_${F}.sorted.bam | ivar consensus -p segment_${i}_${F}_ivar -i ${K}
            python $PIPELINE/Influenza/substitute_degenarate_bases.py segment_${i}_${F}_ivar.fa segments/segment_${i}_${F}.fasta
            rm -rf segment_${i}_${F}_ivar.fa
            rm -rf segments/segment_preconsensus_${i}_${F}.sorted.bam
            if [ $(cat H_genotype.txt) == "H1" ] && [ $i -eq 4 ]; then
                nextclade run -D $PIPELINE/Influenza/nextclade_files/H1 -j ${THREADS} -t nextclade.tsv segments/segment_${i}_${F}.fasta
                csvcut -t -c clade nextclade.tsv | tail -n+2 > clade.txt
            elif [ $(cat H_genotype.txt) == "H3" ] && [ $i -eq 4 ]; then
                nextclade run -D $PIPELINE/Influenza/nextclade_files/H3 -j ${THREADS} -t nextclade.tsv segments/segment_${i}_${F}.fasta
                csvcut -t -c clade nextclade.tsv | tail -n+2 > clade.txt
            elif [ $(cat H_genotype.txt) == "Victoria" ] && [ $i -eq 4 ]; then
                nextclade run -D $PIPELINE/Influenza/nextclade_files/Vic -j ${THREADS} -t nextclade.tsv segments/segment_${i}_${F}.fasta
                csvcut -t -c clade nextclade.tsv | tail -n+2 > clade.txt
            elif [ $(cat H_genotype.txt) == "Yamagata" ] && [ $i -eq 4 ]; then
                nextclade run -D $PIPELINE/Influenza/nextclade_files/Yam -j ${THREADS} -t nextclade.tsv segments/segment_${i}_${F}.fasta
                csvcut -t -c clade nextclade.tsv | tail -n+2 > clade.txt
            elif [ $i -eq 4 ]; then
                echo "nextclade_not_available" > clade.txt
            fi
            bowtie2-build segments/segment_${i}_${F}.fasta segments/segment_${i}_${F}.fasta
            bowtie2 --very-sensitive -p ${THREADS} -x segments/segment_${i}_${F}.fasta -1 segment_${i}_R1.fq.gz -2 segment_${i}_R2.fq.gz | samtools view -S -b > segments/segment_${i}_${F}.bam
            samtools sort  segments/segment_${i}_${F}.bam -o segments/segment_${i}_${F}.sorted.bam
            rm -rf segments/segment_${i}_${F}.bam
            samtools mpileup -aa -A -d 0 -B -Q 0 segments/segment_${i}_${F}.sorted.bam | ivar variants -p segments/segment_${i}_${F}_iVar_variants -t 0.25 -m 10 -r segments/segment_${i}_${F}.fasta
            grep 'TRUE' segments/segment_${i}_${F}_iVar_variants.tsv | grep -vP 'N\t' | grep -vP '\tN' | wc -l > segments/segment_${i}.SNPsCount
            samtools depth -a  segments/segment_${i}_${F}.sorted.bam |  awk '{sum+=$3} END {print sum/NR}' > segments/segment_${i}.MeanDepth
            samtools depth -a  segments/segment_${i}_${F}.sorted.bam  |  awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}' > segments/segment_${i}.MedianDepth
            samtools depth -a  segments/segment_${i}_${F}.sorted.bam|  awk '{print $3 >= 10}' | grep '1' | wc -l > segments/segment_${i}.Depth10
            samtools depth -a  segments/segment_${i}_${F}.sorted.bam |  awk '{print $3 >= 25}' | grep '1' | wc -l > segments/segment_${i}.Depth25
            seqtk comp segments/segment_${i}_${F}.fasta | awk '{x+=$9}END{print x}' > segments/segment_${i}.CountNs
            total_segment=$(seqtk comp  segments/segment_${i}_${F}.fasta | awk '{print $1 "\t" ($3+$4+$5+$6)}' |  cut -f2)
            total_ref=$(seqtk comp segment_${i}_ref.fasta | cut -f2)
            echo "scale=2; ($total_segment / $total_ref) * 100" | bc > segments/segment_${i}.coverage
            printf "segment_${i}_Mean_depth\tsegment_${i}_Median_depth\tsegment_${i}_Npos_Depth>=10\tsegment_${i}_Npos_Depth>=25\tsegment_${i}_Coverage\tsegment_${i}_Number_of_Ns\tsegment_${i}_SNPs\n" > segments/segment_${i}.Statistics
            paste -d "\t" segments/segment_${i}.MeanDepth segments/segment_${i}.MedianDepth segments/segment_${i}.Depth10 segments/segment_${i}.Depth25 segments/segment_${i}.coverage segments/segment_${i}.CountNs segments/segment_${i}.SNPsCount >> segments/segment_${i}.Statistics
            rm -rf segment_${i}_R1.fq.gz segment_${i}_R2.fq.gz
        else # Could not assemble the segment
            printf "segment_${i}_Mean_depth\tsegment_${i}_Median_depth\tsegment_${i}_Npos_Depth>=10\tsegment_${i}_Npos_Depth>=25\tsegment_${i}_Coverage\tsegment_${i}_Number_of_Ns\tsegment_${i}_SNPs\n" > segments/segment_${i}.Statistics
            printf "could_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled\tcould_not_be_assembled" >> segments/segment_${i}.Statistics
            if [ $i -eq 4 ]; then
                echo "not_detected" > clade.txt
            fi
        fi
    else # segment was not detected by vapor
        printf "segment_${i}_Mean_depth\tsegment_${i}_Median_depth\tsegment_${i}_Npos_Depth>=10\tsegment_${i}_Npos_Depth>=25\tsegment_${i}_Coverage\tsegment_${i}_Number_of_Ns\tsegment_${i}_SNPs\n" > segments/segment_${i}.Statistics
        printf "not_detected\tnot_detected\tnot_detected\tnot_detected\tnot_detected\tnot_detected\tnot_detected" >> segments/segment_${i}.Statistics
    fi
done


paste -d "\t" segments/*.Statistics > segments.Statistics

## Recuperar estatísticas gerais para o genoma inteiro
cat segments/*.fasta > Genoma_${F}.fasta

if [ -s Genoma_${F}.fasta ]; then
    bowtie2-build Genoma_${F}.fasta Genoma_${F}.fasta
    bowtie2 --very-sensitive -p ${THREADS} -x Genoma_${F}.fasta -1 ${F}_R1_paired.fq.gz -2 ${F}_R2_paired.fq.gz | samtools view -S -b  > ${F}.bam
    samtools sort  ${F}.bam -o ${F}.sorted.bam
    samtools index ${F}.sorted.bam
    rm -rf ${F}.bam
    samtools view -c -F 260 ${F}.sorted.bam > ${F}.ReadsMappedFinal;
    echo $(zcat ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz | wc -l)/4|bc > ${F}.ReadCount;
    x=$(cat ${F}.ReadsMappedFinal); y=$(cat ${F}.ReadCount); python -c "print(round(float(${x}/${y}*100), 2))" > ${F}.PercentMapped;
    grep -c '>' Genoma_${F}.fasta > ${F}.SegmentsAssembled
    paste H_genotype.txt N_genotype.txt | sed 's/\t//g' > subtype.txt
else
    echo $(zcat ${F}_R1_paired.fq.gz ${F}_R2_paired.fq.gz | wc -l)/4|bc > ${F}.ReadCount
    echo "0" > ${F}.ReadsMappedFinal;
    echo "0" > ${F}.PercentMapped;
    echo "0" > ${F}.SegmentsAssembled
    echo "not_detected" > subtype.txt
fi
ls Genoma_${F}.fasta > ${F}.GenomeName;
printf "Genome\tN_Reads\tReads_assembled\tPercent_assembled\tSegments_Assembled\tType\tH_genotype\tN_genotype\tSubtype\tClade\n" > Genome_basic.Statistics
paste -d "\t" ${F}.GenomeName ${F}.ReadCount ${F}.ReadsMappedFinal ${F}.PercentMapped ${F}.SegmentsAssembled A_or_B_genotype.txt H_genotype.txt N_genotype.txt subtype.txt clade.txt >> Genome_basic.Statistics

paste -d "\t" Genome_basic.Statistics segments.Statistics > ${F}_complete.Statistics

# Copy Genome and Statistics to folder up.
cp Genoma_${F}.fasta ../;
cp ${F}_complete.Statistics ../;
rm -rf *.fq.gz
# Back to folder up.
find . -type d | xargs -i chmod g+w {}
cd ..

