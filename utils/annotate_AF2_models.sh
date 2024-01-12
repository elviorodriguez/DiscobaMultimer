#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
	echo "Usage: $0 <annotations_file> <AF2_folder_path>"
 	echo ""
  	echo "Used to rename all the filesystem replacing the IDs for matching protein symbols."
   	echo ""
   	echo "Parameters:"
  	echo "  - annotations_file: TSV file containg the IDs to protein symbols correspondance."
   	echo "                      Take a look at utils/annotations.txt file to see its format."
   	echo "  - AF2_folder_path : Path to AF2 folder created by DiscobaMultimer that you want "
        echo "                      to rename the filesystem."
	exit 1
fi

# Assign the arguments to variables
annotations_file="$1"
folder_path="$2"

# Check if the provided folder path exists
if [ ! -d "$folder_path" ]; then
	echo "Error: Folder '$folder_path' not found."
	exit 1
fi

# Check if the provided annotations file exists
if [ ! -f "$annotations_file" ]; then
	echo "Error: File '$annotations_file' not found."
	exit 1
fi

# Generate file to iterate over
ls -1 $folder_path > AF2_dirnames.tmp

# Rename directories based on Symbol value
while read -r AF2_prediction; do
	new_name=""
	# Separate in individual IDs
	IDs=($(echo "$AF2_prediction" | tr "__vs__" " "))
	for ID in "${IDs[@]}"; do
		annotation=$(grep "$ID" $annotations_file)
		symbol=${annotation##*$'\t'}
		if [[ "$symbol" == "N/A" ]]; then
			new_name=${new_name}__vs__${ID}
		else
			new_name=${new_name}__vs__${symbol}
		fi
	done
	new_name=${new_name#*__vs__}
	if [[ "$new_name" == "$AF2_prediction" ]]; then
		echo $AF2_prediction IDs has no annotations.
	else
		ls -1 $folder_path/$AF2_prediction > AF2_filenames.tmp
		while read -r AF2_filename; do
			new_filename=${AF2_filename//$AF2_prediction/$new_name}
			mv $folder_path/$AF2_prediction/$AF2_filename $folder_path/$AF2_prediction/$new_filename
		done < AF2_filenames.tmp
		rm AF2_filenames.tmp
		mv "$folder_path/$AF2_prediction" "$folder_path/$new_name"
	fi

done < AF2_dirnames.tmp

# Remove temporary file
rm AF2_dirnames.tmp

