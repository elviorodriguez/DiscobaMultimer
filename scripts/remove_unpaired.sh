#!/bin/bash

# This script removes the unpaired part of a paired+unpaired a3m file
# The a3m file needs to have the CARDINALITY as first line
# (eg. #232,162	1,1)

# USAGE:
# ./remove_unpaired.sh path/to/file.a3m
# 
# OUTPUT:
# generates "path/to/file_paired.a3m" with only paired sequences

# Check for possitional argument
if [ -z "$1" ]; then
	echo "USAGE: $0 path/to/file.a3m"
	echo "OUTPUT: path/to/file_paired.a3m"
	exit 1
fi

input_a3m=$1

# Test if a3m file starts with cardinality -----------------------------------
# Regular expression of the cardinality
card_regex="^#[0-9]\{1,5\},[0-9]\{1,5\}[[:space:]][0-9]\{1,2\},[0-9]\{1,2\}$"
# Use head to extract the first line of the file and pipe it to grep
if head -1 $input_a3m | grep -q $card_regex; then
	echo "a3m test PASS: a3m file contains cardinality"
else
	echo "a3m test NOT PASSED: a3m file does not contain cardinality (e.g. #123,523	1,1)"
	echo "USAGE: $0 path/to/file.a3m"
	echo "OUTPUT: path/to/file_paired.a3m"
	exit 1
fi

# Basename
input_a3m_dirname=`dirname $input_a3m`
input_a3m_basename=`basename $input_a3m`
input_a3m_name=${input_a3m_basename%.a3m}

# Protein IDs
protein1_ID=`echo $input_a3m_name | awk -F '__vs__' '{print $1}'`
protein2_ID=`echo $input_a3m_name | awk -F '__vs__' '{print $2}'`

# Extract the length of the first (L1) and second (L2) sequences
first_seq_length=`head -1 $input_a3m | cut -d "#" -f 2 | cut -d "," -f 1`
second_seq_length=`head -1 $input_a3m | cut -d "#" -f 2 | cut -d "," -f 2 | cut -d $'\t' -f 1`
first_seq_length=$((first_seq_length))		# convert to number
second_seq_length=$((second_seq_length))	# convert to number

# Original length of the MSA
MSA_length=`grep -c '^>' $input_a3m`
MSA_length=$((MSA_length))			# convert to number

# Output some data
echo "-------------------------------------"
echo "Input a3m file dirname: $input_a3m_dirname"
echo "Input a3m file basename: $input_a3m_basename"
echo "Input a3m file name: $input_a3m_name"
echo "Protein 1 ID: $protein1_ID"
echo "Protein 2 ID: $protein2_ID"
echo "Protein 1 length (L1): $first_seq_length"
echo "Protein 2 length (L2): $second_seq_length"
echo "Original number of sequences: $MSA_length"

# Remove unpaired ----------------------------------------

# Progress
echo "Starting to process a3m file"

# Output file creation
output_a3m="${input_a3m_dirname}/${input_a3m_name}_paired.a3m"
head -1 $input_a3m > $output_a3m

# Open the file for reading
exec 3< $input_a3m

# Read and discard the first line
read -r <&3

# Loop over the file in steps of two lines (seq by seq)
while read -r HEADER && read -r SEQUENCE; do
	#echo "Line 1: $HEADER"		# Debug
	#echo "Line 2: $SEQUENCE"	# Debug
	if [[ $SEQUENCE =~ ^-{$first_seq_length}|-{$second_seq_length}$ ]]; then
		echo "Unpaired sequence detected. Removing it."
	else
		echo $HEADER >> $output_a3m
		echo $SEQUENCE >> $output_a3m
	fi
done <&3

# Close the file
exec 3<&-

# Output info
final_MSA_length=`grep -c '^>' $output_a3m`
echo "FINISH"
echo "Output file: $output_a3m"
echo "Remaining paired sequences: $final_MSA_length"
