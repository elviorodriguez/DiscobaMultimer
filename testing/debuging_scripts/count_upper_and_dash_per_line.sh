#!/bin/bash

# Check if a file is provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# Check if the provided argument is a file
if [ ! -f "$1" ]; then
    echo "Error: $1 is not a valid file"
    exit 1
fi

# Read the file line by line
while IFS= read -r line
do
    # Count uppercase letters and hyphens in the line
    count=$(echo "$line" | grep -o '[A-Z-]' | wc -l)
    # Print the line and the count
    echo "$line"
    echo "COUNT: $count"
done < "$1"

