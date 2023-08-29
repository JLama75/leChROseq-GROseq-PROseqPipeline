#!/usr/bin/env bash
#PRO-seq data from Abuhashem et al. 2022 (https://genesdev.cshlp.org/content/36/13-14/770)
#Only taking mouse Embryonic stem cells (mESCs) that are untreated controls 
#the ESCs were derived from male embryo

# Declare the array #list alternative in bash
#cbsuGpu01 only for downloading rest of the analysis done in cbsudanko machine

scratch = /workdir/jl3285/mESCs

declare -a accessions=(
    "SRR18010270"
    "SRR18010271"
    "SRR18010272"
    "SRR18010275"
    "SRR18010276"
    "SRR18010277"
    "SRR18010278"
    "SRR18010279"
    "SRR18010280"
)

count=0
for i in "${accessions[@]}"; do
    echo "Downloading ${i}"
    prefetch ${i}
    fasterq-dump ${i}/${i}.sra
    cp ${i}/${i}_1.fastq ${i}/${i}_2.fastq -t ${scratch}
    gzip ${i}_1.fastq ${i}_2.fastq  #SRR1482439_1.fastq  SRR1482439_2.fastq 
    count=$((count + 1))
    rm -r ${i}
done

echo "downloaded ${count} files"
echo "done!"

#!/bin/bash
# Merge the technical replicates for both biological replicates
echo "Merging the technical replicates for both biological replicates..."

zcat SRR18010270_1.fastq.gz SRR18010271_1.fastq.gz SRR18010272_1.fastq.gz > Abood123_R1_spikeIn.fastq
zcat SRR18010270_2.fastq.gz SRR18010271_2.fastq.gz SRR18010272_2.fastq.gz > Abood123_R2_spikeIn.fastq

zcat SRR18010275_1.fastq.gz SRR18010276_1.fastq.gz SRR18010277_1.fastq.gz > Abood19.20.21_R1_spikeIn.fastq
zcat SRR18010275_2.fastq.gz SRR18010276_2.fastq.gz SRR18010277_2.fastq.gz > Abood19.20.21_R2_spikeIn.fastq

zcat SRR18010278_1.fastq.gz SRR18010279_1.fastq.gz SRR18010280_1.fastq.gz > Abood4_R1_spikeIn.fastq
zcat SRR18010278_2.fastq.gz SRR18010279_2.fastq.gz SRR18010280_2.fastq.gz > Abood4_R2_spikeIn.fastq

echo "Merging completed!"

gzip Abood*


###################################################
#cbsudanko machine

scratch=/workdir/jl3285/mESCs/PRO-seq_Abuhassem2022
cd ${scratch}

mkdir -p Mapped_dir
mkdir -p FPKM

for file in *1.fastq.gz #Abood123_spikeIn_1.fastq, Abood19.20.21_spikeIn_1.fastq Abood4_spikeIn_1.fastq

do

NAMES=$(echo "$file" | cut -d _ -f 1  | sort | uniq)  
echo processing $file  #Abood123_spikeIn_1.fastq, Abood19.20.21_spikeIn_1.fastq Abood4_spikeIn_1.fastq

echo processed the names into "$NAMES"... #Abood123, Abood19.20.21, Abood4
#Since the fastq files have both dm3 spikein and mouse data, we will perform “competitive alignment” 
# which results in the creation of BAM files containing a mix of chromosomes

export mouse_genome="/local/storage/data/short_read_index/mm10/bwa.indexes_mm10.rRNA.dm3/mm10.dm3" #indices with both mm10 and dm3
export mouse_chinfo="/local/storage/data/short_read_index/mm10/bwa.indexes_mm10.rRNA.dm3/mm10.dm3.chromInfo"

#PREFIX="Abood123"

export PATH=/programs/kentUtils/bin:$PATH
export PATH=/programs/seqtk:$PATH
export PATH=/programs/prinseq-lite-0.20.4:$PATH
export PATH=/programs/proseq2.0:$PATH
export PATH=/programs/cutadapt-4.1/bin:$PATH
export PYTHONPATH=/programs/cutadapt-4.1/lib/python3.9/site-packages:/programs/cutadapt-4.1/lib64/python3.9/site-packages

#does the demultiplexing throug princelite and adapter cutting along with competitive alignment with mouse and drosophila genome. 
#The fastq files has drosophila DNA as a spike in.

#changing the name of the files
#mv ${NAMES}_R1_spikeIn.fastq ${NAMES}_spikeIn_R1.fastq
#mv ${NAMES}_R2_spikeIn.fastq ${NAMES}_spikeIn_R2.fastq
#mv ${NAMES}_spikeIn_1.fastq.gz ${NAMES}_spikeIn_R1.fastq.gz
#mv ${NAMES}_spikeIn_2.fastq.gz ${NAMES}_spikeIn_R2.fastq.gz

echo using proseq2.0 to align ${NAMES}_spikeIn_R1.fastq.gz ${NAMES}_spikeIn_R2.fastq.gz
bash proseq2.0.bsh --RNA5=R2_5prime -PE --UMI1=6 -i ${mouse_genome} -c ${mouse_chinfo} -I ${NAMES}_spikeIn -T ./Mapped_dir/ -O ./Mapped_dir/ --thread=15
#output: Abood123_spikeIn_dedup_QC_end.sort.bam 

#PROseq sorts by name
echo indexing the bam files #sorting by chromosome 
samtools sort -@ 10 Mapped_dir/${NAMES}_spikeIn_dedup_QC_end.sort.bam -o Mapped_dir/${NAMES}_QC_end.sorted.bam
samtools index -@ 8 Mapped_dir/${NAMES}_QC_end.sorted.bam

echo -e "\n filtering quality reads, removing rRNA and converting bam to bed..." 
echo -e "\n reporting only single nucleotide position..." 
        
bedtools bamtobed -i Mapped_dir/${NAMES}_QC_end.sorted.bam | awk 'BEGIN{{OFS="\t"}} ($5 > 20){{print $0}}' | grep -v "rRNA" | grep -v "chrM" | awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' | gzip > FPKM/${NAMES}.bed.gz

# Removing the dm3 portion 
zcat FPKM/${NAMES}.bed.gz | grep -v "dm3" | gzip > FPKM/${NAMES}.dm3removed.bed.gz

rm Mapped_dir/${NAMES}_QC_end.sorted.bam 
gzip Mapped_dir/${NAMES}_spikeIn_dedup_QC_end.sort.bam

echo done!

done

##############################################################################################################################
cd /workdir/jl3285/mESCs/FPKM
#To convert gtf to B12 format that is need for FPKM calculation downstream
source /home/jl3285/miniconda3/bin/activate
#conda create -n bedparse -c bioconda bedparse
conda activate bedparse #downloading the package called bedparse

gzip -d gencode.vM25.annotation.gtf.gz
bedparse gtf2bed gencode.vM25.annotation.gtf > gencode.vM25.annotation.bed
gzip gencode.vM25.annotation.gtf

##############################################################################################################################

export mm10_chinfo=/local/storage/data/mm10/mm10.chromInfo
echo $mm10_chinfo

for file in *dm3removed.bed.gz #Abood123.dm3removed.bed.gz
do

NAMES=$(echo "$file" | cut -d . -f 1  | sort | uniq)  #Abood123
echo processing $file  #Abood123.dm3removed.bed.gz

echo processed the names into "$NAMES" #Abood123

echo converting bed to bam...

gzip -d ${NAMES}.dm3removed.bed.gz
bedtools bedtobam -i ${NAMES}.dm3removed.bed -g ${mm10_chinfo} > ${NAMES}.filtered.bam 

echo sorting and creating index files from bam files
samtools sort ${NAMES}.filtered.bam -o ${NAMES}.sorted.filtered.bam
samtools index ${NAMES}.sorted.filtered.bam

rm ${NAMES}.filtered.bam 

Rscript getRefGenes.R

echo calculating the FPKM ...
echo activating RSeQC ...

export PYTHONPATH=/programs/RSeQC-5.0.1/lib64/python3.9/site-packages:/programs/RSeQC-5.0.1/lib/python3.9/site-packages
export PATH=/programs/RSeQC-5.0.1/bin:$PATH

#FPKM_count.py -i ${NAMES}.sorted.filtered.bam  -r Ref_bodies.bed –d ‘1++,1–-,2+-,2-+’ -o ${NAMES}
FPKM_count.py -i ${NAMES}.sorted.filtered.bam  -r Ref_bodies.bed -o ${NAMES}
#IF YOU DON'T KNOW THE STRANDED NESS RUN 
infer_experiment.py -i  ${NAMES}.sorted.filtered.bam -r Ref_bodies.bed > ${NAMES}.infer.txt #unstranded #single end: Paired end information might have been lost when converting bed to bam format
    
echo done

done

echo removing files
rm *bam 
