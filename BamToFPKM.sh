#!/bin/bash
#Converting the bam files from leChRO-seq to FPKM values for three spermatocytes stages: Lepto/Zygotene, Pachytene and Zygotene.
#Refgenome used: mm10
#To convert gtf to B12 format that is need for FPKM calculation downstream

source /home/jl3285/miniconda3/bin/activate
conda create -n bedparse -c bioconda bedparse
conda activate bedparse #downloading the package called bedparse

gzip -d gencode.vM25.annotation.gtf.gz
bedparse gtf2bed gencode.vM25.annotation.gtf > gencode.vM25.annotation.bed
gzip gencode.vM25.annotation.gtf

#check to see if all the genes have 6 columns
awk '{if ($6=="+"|| $6=="-") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' gencode.vM31.annotation.bed| wc -l #1122 rows have gene name missing
awk '{if ($6=="+"|| $6=="-") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' gencode.vM25.annotation.bed| wc -l #5 rows have gene name missing
awk '$6!="+" || $6!="-" {print}' gencode.vM25.annotation.bed > gencode.vM25.annotation.filtered.bed


##############################################################################################################################
#Using bam files to generate FPKM counts
source /home/jl3285/miniconda3/bin/activate 

cd /local/storage/projects/sexChromosome_dREG_ATAC/ChRO-seq/bamFiles

cp *.bam.gz -t /workdir/jl3285/ChRO_seq/mm10geneAnnotation/FPKM

scratch=/workdir/jl3285/ChRO_seq/mm10geneAnnotation/FPKM

cd $scratch
echo $scratch

gzip -d *bam.gz
export mm10_chinfo=/local/storage/data/mm10/mm10.chromInfo

echo $mm10_chinfo

for file in *.bam #ChRO_P_R2.sorted.bam.gz
do

NAMES=$(echo "$file" | cut -d . -f 1 | cut -d _ -f 1-3 | sort | uniq)   
echo processing $file #ChRO_P_R2.sorted.bam.gz

echo processed the names into "$NAMES" #ChRO_D_R1_minus

echo converting bam to bed and filter sorting... 

bedtools bamtobed -i ${NAMES}.sorted.bam | awk 'BEGIN{OFS="\t"} ($5 > 20){print $0}' | grep -v "rRNA" | awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' > ${NAMES}.sorted.bed

echo converting bed to bam...

bedtools bedtobam -i ${NAMES}.sorted.bed -g ${mm10_chinfo} > ${NAMES}.filtered.bam 

echo sorting and creating index files from bam files
samtools sort ${NAMES}.filtered.bam -o ${NAMES}.sorted.filtered.bam
samtools index ${NAMES}.sorted.filtered.bam

Rscript getRefGenes.R

echo calculating the FPKM ...
echo activating RSeQC ...

#conda activate rseqc

export PYTHONPATH=/programs/RSeQC-5.0.1/lib64/python3.9/site-packages:/programs/RSeQC-5.0.1/lib/python3.9/site-packages

export PATH=/programs/RSeQC-5.0.1/bin:$PATH

FPKM_count.py -i ${NAMES}.sorted.filtered.bam  -r Ref_bodies.bed –d ‘1++,1–,2+-,2-+’ -o ${NAMES}.FPKM 

#cat ${NAMES}.FPKM | awk print'{$1, $2, $3, $4, $5, $6, $7, $8, $9*NORMRAD }'
echo done

done

echo removing files
rm *bam  *sorted.bed 
