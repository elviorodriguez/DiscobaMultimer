#!/bin/bash

AF2_folder=$1
ID_pair_list=$2

while read line; do
	IFS=$'\t' read -r -a IDs_array <<< "$line"
	paired_name=`printf '%s__vs__' "${IDs_array[@]}" | sed 's/__vs__$//'`
	expected_AF2_dirname=$paired_name
	ls -d $AF2_folder/$expected_AF2_dirname
	if [ $? -ne 0 ]; then
		reversed_IDs_array=("${IDs_array[@]:1}" "${IDs_array[0]}")
		reversed_paired_name=`printf '%s__vs__' "${reversed_IDs_array[@]}" | sed 's/__vs__$//'`
		expected_reversed_AF2_dirname=$reversed_paired_name
		ls -d $AF2_folder/$expected_reversed_AF2_dirname
		if [ $? -eq 0 ]; then
			echo "IDs_backwards: $line"
			echo "Rearange them if you're going to work with this AF2 folder model."
		fi
	fi

done < "$ID_pair_list"

