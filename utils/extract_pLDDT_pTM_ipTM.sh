#!/bin/bash

AF2_models_folder=$1
output_file=$2

echo "ID1 ID2 rank pLDDT pTM ipTM" > $output_file

grep "rank_00" $AF2_models_folder/*/log.txt | \
	grep -oE '[^/]+__vs__[^/ ]+|rank_00[0-9.]+|pLDDT=[0-9.]+|pTM=[0-9.]+|ipTM=[0-9.]+' | \
	awk 'ORS=NR%5?FS:RS' | \
	awk '{gsub("__vs__"," "); gsub("rank_00",""); gsub("pLDDT=",""); gsub(" pTM="," "); gsub("ipTM=",""); print}' >> $output_file
