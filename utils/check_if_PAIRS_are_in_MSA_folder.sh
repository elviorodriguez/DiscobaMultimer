#!/bin/bash

msa_folder=$1
ID_list=$2

while read line; do
	IFS=$'\t' read -r -a IDs_array <<< "$line"
	paired_name=`printf '%s__vs__' "${IDs_array[@]}" | sed 's/__vs__$//'`
	expected_a3m_filename=$paired_name.a3m
	ls $msa_folder/$expected_a3m_filename
done < "$ID_list"

