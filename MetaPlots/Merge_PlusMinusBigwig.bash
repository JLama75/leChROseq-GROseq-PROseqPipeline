
# Merge the plus and minus strand of bigwig to generate the metaplot

scratch=/workdir/jl3285/Adriana_metaplot/ChroSeq_plusMinus

final=/local/storage/projects/prophaseI/Jyoti/merged

cp *rad.bw $scratch
export mm10_chinfo=/local/storage/data/mm10/mm10.chromInfo
echo $mm10_chinfo

cd $scratch
echo $scratch

source /home/jl3285/miniconda3/bin/activate bigWigToBedGraph
echo bigWigToBedGraph activated

for file in `cat FileNames.txt` #FileNames.txt has stored all the names of required files. 
do

NAMES=$(echo "$file" | cut -d _ -f 1-2 | sort | uniq)   
echo processing $file # reading all the names of files one by one (eg:LZ_merged_all_minus.rpkm.rad.bw)
echo processed the names into "$NAMES" #LZ_merged (taking part of the name we want)

echo converting bigwig to bedGraph...
bigWigToBedGraph ${NAMES}_all_plus.rpkm.rad.bw ${NAMES}.tmp.plus.bedGraph  #check wc -l of tmp file
bigWigToBedGraph ${NAMES}_all_minus.rpkm.rad.bw ${NAMES}.tmp.minus.bedGraph #check wc -l of tmp file

echo merging the plus and minus strand data using bedtools unionbedg...
bedtools unionbedg -i ${NAMES}.tmp.plus.bedGraph ${NAMES}.tmp.minus.bedGraph > tmp.merge.bedGraph #compare the #check wc -l of merge vs tmp file
echo -e "line count of ${NAMES}.tmp.merge.bedGraph is:\n"
wc -l tmp.merge.bedGraph >> out

echo deleting the minus sign from read count obtained from minus strand...
echo adding the count reads of plus and minus strand...

Rscript GetAbsoluteValue_Sum.R # check the output. Its wc -l should be equal to merge file 

echo -e "line count of ${NAMES}.Output.bedGraph is:\n"
wc -l Output.bedGraph >> out

cat Output.bedGraph | awk '{print $2, $3, $4, $5}' | tr -d \" > ${NAMES}_sum.bedGraph # check the sum.bedGraph file   #LZ_merged_sum.bedGraph

done


source /home/jl3285/miniconda3/bin/activate bedGraphToBigWig
echo bedGraphToBigWig activated


for file in *sum.bedGraph
do

NAMES=$(echo "$file" | cut -d _ -f 1-2 | sort | uniq)  

echo converting ${NAMES}_sum.bedGraph bedGraph to bigwig...
bedGraphToBigWig ${NAMES}_sum.bedGraph ${mm10_chinfo} ${NAMES}.ChroSeq.bw 

cp ${NAMES}.ChroSeq.bw $final

done

echo done copying!
