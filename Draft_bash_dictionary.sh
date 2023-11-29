#!/bin/bash
declare -A dict #declaring dictionary

while IFS=$'\t' read -r stage norm 	 # key=stage; value=norm ; stage= #ChRO_D_R1  value=###
do

dict[${stage}]=${norm}		#dictionary[$key]=$value
echo within while loop processing $norm 

done < norm.files

echo outside processing ${dict[ChRO_D_R1]} and ${dict[ChRO_LZ_R1]}

for file in *_plus_bedGraph.gz
do

NAMES=$(echo "$file" | cut -d _ -f 1-3 | sort | uniq) #ChRO_D_R1_minus_bedGraph.gz
echo processing the names of the $file into $NAMES #ChRO_D_R1
#echo associative array ${dict[*]}
echo reading dictionary for $NAMES : ${dict[${NAMES}]}

done
