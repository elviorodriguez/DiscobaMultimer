#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 <ColabFold_MSA.a3m> <DiscobaDB_MSA.a3m> <output_name>"
    exit 1
}

# Check if exactly three arguments are provided
if [ "$#" -ne 3 ]; then
    usage
fi

ColabFold_MSA=$1
DiscobaDB_MSA=$2
output_name=$3

# Check if the input files have .a3m extension
if [[ ! "$ColabFold_MSA" == *.a3m || ! "$DiscobaDB_MSA" == *.a3m ]]; then
    echo "Error: Both input files must have the .a3m extension"
    usage
fi

# Check if the input files exist
if [ ! -f "$ColabFold_MSA" ]; then
    echo "Error: $ColabFold_MSA does not exist"
    exit 1
fi

if [ ! -f "$DiscobaDB_MSA" ]; then
    echo "Error: $DiscobaDB_MSA does not exist"
    exit 1
fi

# Use the discoba paired+unpaired MSA
cat "$ColabFold_MSA" > "${output_name}.a3m"
echo "" >> "${output_name}.a3m"
tail -n +2 "$DiscobaDB_MSA" >> "${output_name}.a3m"

echo "Output written to ${output_name}.a3m"

