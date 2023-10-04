
#!/usr/bin/env bash
#GRO-seq data from Adelman's lab (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402150/) Williams et al. 2015
#The mouse Embryonic Stem Cells (mESCs) are derived from the female embryo.
# Declare the array #list alternative in bash
#Ilumina MiSeq: 2-3 billion reads

declare -a accessions=(
    "SRR1482440"
    "SRR1482439"
    "SRR1482438"
    "SRR1482437"
    "SRR1482436"
    "SRR1482441"
    "SRR1482442"
    "SRR1482443"
    "SRR1482444"
)

count=0
for i in "${accessions[@]}"; do
    echo "Downloading ${i}"
    prefetch ${i}
    fasterq-dump ${i}/${i}.sra
    gzip ${i}_1.fastq ${i}_2.fastq  #SRR1482439_1.fastq  SRR1482439_2.fastq
    count=$((count + 1))
    rm -r ${i}
done

echo "downloaded ${count} files"
echo "done!"

# Merge the technical replicates for both biological replicates
echo "Merging the technical replicates for both biological replicates..."

zcat SRR1482440_1.fastq.gz SRR1482439_1.fastq.gz SRR1482438_1.fastq.gz SRR1482437_1.fastq.gz SRR1482436_1.fastq.gz > adelman_Rep1_R1.fastq
zcat SRR1482440_2.fastq.gz SRR1482439_2.fastq.gz SRR1482438_2.fastq.gz SRR1482437_2.fastq.gz SRR1482436_2.fastq.gz > adelman_Rep1_R2.fastq

zcat SRR1482441_1.fastq.gz SRR1482442_1.fastq.gz SRR1482443_1.fastq.gz SRR1482444_1.fastq.gz > adelman_Rep2_R1.fastq
zcat SRR1482441_2.fastq.gz SRR1482442_2.fastq.gz SRR1482443_2.fastq.gz SRR1482444_2.fastq.gz > adelman_Rep2_R2.fastq

echo "Merging completed!"

#Cutadapt...
###################################################
export PATH=/programs/cutadapt-4.1/bin:$PATH
export PYTHONPATH=/programs/cutadapt-4.1/lib/python3.9/site-packages:/programs/cutadapt-4.1/lib64/python3.9/site-packages

#Cut the adapters
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAAC -e 0.10 --minimum-length=15 --cores=5 -o adelman_Rep1_R1.no-A.fastq -p adelman_Rep1_R2.no-A.fastq adelman_Rep1_R1.fastq adelman_Rep1_R2.fastq
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAAC -e 0.10 --minimum-length=15 --cores=5 -o adelman_Rep2_R1.no-A.fastq -p adelman_Rep2_R2.no-A.fastq adelman_Rep2_R1.fastq adelman_Rep2_R2.fastq


#Cut the polyA tails
cutadapt -a AAAAAAAAAAAAAAAAAAAA -G TTTTTTTTTTTTTTTTTTTT -e 0.10 --minimum-length=15 --cores=5 -o adelman_Rep1_noAdapters_R1.fastq -p adelman_Rep1_noAdapters_R2.fastq adelman_Rep1_R1.no-A.fastq adelman_Rep1_R2.no-A.fastq
cutadapt -a AAAAAAAAAAAAAAAAAAAA -G TTTTTTTTTTTTTTTTTTTT -e 0.10 --minimum-length=15 --cores=5 -o adelman_Rep2_noAdapters_R1.fastq -p adelman_Rep2_noAdapters_R2.fastq adelman_Rep2_R1.no-A.fastq adelman_Rep2_R2.no-A.fastq

#Mapping to the mm10 reference genome using proseq2.0 from Charles Danko's lab. (https://github.com/Danko-Lab/proseq2.0.git)

mkdir -p Mapped_dir
mkdir -p FPKM

export mouse_genome="/local/storage/data/short_read_index/mm10/bwa.rRNA-0.7.8-r455/mm10.rRNA.fa.gz"
export mouse_chinfo="/local/storage/data/mm10/mm10.chromInfo"
#PREFIX="adelman_Rep1_noAdapters"

export PATH=/programs/kentUtils/bin:$PATH
export PATH=/programs/seqtk:$PATH
export PATH=/programs/prinseq-lite-0.20.4:$PATH
export PATH=/programs/proseq2.0:$PATH
export PATH=/programs/cutadapt-4.1/bin:$PATH
export PYTHONPATH=/programs/cutadapt-4.1/lib/python3.9/site-packages:/programs/cutadapt-4.1/lib64/python3.9/site-packages

echo using proseq2.0 to align adelman_Rep1_noAdapters_R1.fastq adelman_Rep1_noAdapters_R2.fastq
bash proseq2.0.bsh --RNA5=R1_5prime -PE --ADAPT1=GATCGTCGGACTGTAGAACTCTGAAC --ADAPT2=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --map5=FALSE  -i ${mouse_genome} -c ${mouse_chinfo} -I adelman_Rep1_noAdapters  -T ./Mapped_dir/ -O ./Mapped_dir/ --thread=12

#PREFIX="adelman_Rep2_noAdapters"
echo using proseq2.0 to align adelman_Rep2_noAdapters_R1.fastq adelman_Rep2_noAdapters_R2.fastq
bash proseq2.0.bsh --RNA5=R1_5prime -PE --ADAPT1=GATCGTCGGACTGTAGAACTCTGAAC --ADAPT2=TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --map5=FALSE  -i ${mouse_genome} -c ${mouse_chinfo} -I adelman_Rep2_noAdapters  -T ./Mapped_dir/ -O ./Mapped_dir/ --thread=12

cd /workdir/jl3285/GRO-seq_mESCs

echo indexing the bam files
gzip -d  Mapped_dir/adelman_Rep1_noAdapters_QC_end.sort.bam.gz
gzip -d  Mapped_dir/adelman_Rep2_noAdapters_QC_end.sort.bam.gz

samtools sort -@ 10 Mapped_dir/adelman_Rep1_noAdapters_QC_end.sort.bam -o Mapped_dir/adelman_Rep1_noAdapters_QC_end.sorted.bam
samtools sort -@ 10 Mapped_dir/adelman_Rep2_noAdapters_QC_end.sort.bam -o Mapped_dir/adelman_Rep2_noAdapters_QC_end.sorted.bam

samtools index -@ 8 Mapped_dir/adelman_Rep1_noAdapters_QC_end.sorted.bam
samtools index -@ 8 Mapped_dir/adelman_Rep2_noAdapters_QC_end.sorted.bam

echo -e "\n filtering quality reads, removing rRNA, chromosome M and converting bam to bed..." 
        
bedtools bamtobed -i Mapped_dir/adelman_Rep1_noAdapters_QC_end.sorted.bam | awk 'BEGIN{{OFS="\t"}} ($5 > 20){{print $0}}' | grep -v "rRNA" | grep -v "chrM" > FPKM/adelman_Rep1.bed
bedtools bamtobed -i Mapped_dir/adelman_Rep2_noAdapters_QC_end.sorted.bam | awk 'BEGIN{{OFS="\t"}} ($5 > 20){{print $0}}' | grep -v "rRNA" | grep -v "chrM" > FPKM/adelman_Rep2.bed

gzip Mapped_dir/adelman_Rep1_noAdapters_QC_end.sort.bam
gzip Mapped_dir/adelman_Rep2_noAdapters_QC_end.sort.bam

rm Mapped_dir/adelman_Rep*_noAdapters_QC_end.sorted.bam 

#####################################

cd FPKM
#scratch=/workdir/jl3285/GRO-seq_mESCs/FPKM
#echo $scratch

#Reporting only single nucleotide position information of where the read maps to/RNApolII is bound
cat adelman_Rep1.bed | awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' > adelman_Rep1.sorted.bed
cat adelman_Rep2.bed | awk 'BEGIN{OFS="\t"} ($6 == "+") {print $1,$2,$2+1,$4,$5,$6}; ($6 == "-") {print $1,$3-1,$3,$4,$5,$6}' > adelman_Rep2.sorted.bed


##############################################################################################################################

#To convert gtf to B12 format that is needed for FPKM calculation downstream
source /home/jl3285/miniconda3/bin/activate
conda activate bedparse #downloading the package called bedparse

gzip -d gencode.vM25.annotation.gtf.gz
bedparse gtf2bed gencode.vM25.annotation.gtf > gencode.vM25.annotation.bed
gzip gencode.vM25.annotation.gtf

#check to see if all the genes have 6 columns
awk '{if ($6=="+"|| $6=="-") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' gencode.vM25.annotation.bed| wc -l #no rows have gene name missing
awk '$6!="+" || $6!="-" {print}' gencode.vM25.annotation.bed > gencode.vM25.annotation.filtered.bed
 
 #done!!
##############################################################################################################################
#Converting bed to bam format as FPKM_count.py code uses bam as input to calculate FPKM values
export mm10_chinfo=/local/storage/data/mm10/mm10.chromInfo
echo $mm10_chinfo

for file in adelman_Rep*.sorted.bed #adelman_Rep1.sorted.bed
do

NAMES=$(echo "$file" | cut -d . -f 1  | sort | uniq)   
echo processing $file #adelman_Rep1.sorted.bed

echo processed the names into "$NAMES" #adelman_Rep1

echo converting bed to bam...

bedtools bedtobam -i ${NAMES}.sorted.bed -g ${mm10_chinfo} > ${NAMES}.filtered.bam 

echo sorting and creating index files from bam files
samtools sort ${NAMES}.filtered.bam -o ${NAMES}.sorted.filtered.bam
samtools index ${NAMES}.sorted.filtered.bam

Rscript getRefGenes.R #This generates a text file with all gene bodies and their positions needed as input for FPKM_count.py program

echo calculating the FPKM ...
echo activating RSeQC ... #FPKM_count.py code is in RSeQC package

export PYTHONPATH=/programs/RSeQC-5.0.1/lib64/python3.9/site-packages:/programs/RSeQC-5.0.1/lib/python3.9/site-packages
export PATH=/programs/RSeQC-5.0.1/bin:$PATH

FPKM_count.py -i ${NAMES}.sorted.filtered.bam  -r Ref_bodies.bed –d -o ${NAMES}
#This will give you the Excel file with gene bodies and their FPKM values. 

#IF YOU DON'T KNOW THE STRANDED NESS RUN 
infer_experiment.py -i  ${NAMES}.sorted.filtered.bam -r Ref_bodies.bed > ${NAMES}.infer.txt #unstranded #single-end: paired end information is possibly lost when converting bed to bam.

echo done

done

echo removing files
rm *bam 
