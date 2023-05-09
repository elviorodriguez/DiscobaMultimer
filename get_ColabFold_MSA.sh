#!/bin/bash -e

# This script runs MMseqs2 to get Discoba MSA 

# DEPENDENCIES -------------------------------------------------------------
# colabfold_batch
# --------------------------------------------------------------------------

WD=`pwd`

usage() {
	echo "USAGE: $0 <database.fasta> <protein_ID1> [<protein_IDn>]"
	echo "OPTIONS:"
	echo "	- At least one protein ID must be provided."
	echo "OUTPUT: stored in ./colabfold_MSA directory"
	echo "	- protein_ID[__vs__protein_IDn].a3m"
	echo "NOTES:"
	echo "	- MMseqs2 and DiscobaDB must be installed in advance"
	echo "	  by running install_MMseqs2_and_DiscobaDB.sh"
	exit 1
}

# Check input ---------------------------------------------------------
if [ $# -lt 2 ]; then
	echo "ERROR: 2 or more arguments are needed"
	usage
elif [ ! -f $1 ]; then
	echo "ERROR: the database argument $1 is not a file"
	usage
fi

# Assign positional arguments to variables
database=$1
IDs_number=$(($# - 1))
IDs_array=()

echo "Get MSA using ColabFold custom MSA server"
echo "Number of IDs: $IDs_number"
# Assing each ID to a variable
for i in $(seq 2 $#); do
	# Create a variable with a dynamic name
	var_name="ID_$((i - 1))"
	# Get the value of the corresponding positional argument
	var_value=${!i}
	# Assign the value to the variable
	eval "$var_name=$var_value"
	echo "	${var_name}: ${var_value}"
	IDs_array+=(${var_value})
done



# Aux Fx ---------------------------------------------------------------
# Searches and finds $ID in $database. Outputs ${ID}.fasta in the WD
get_sequence() {
	database=$1
	ID=$2
	# Search for the ID in the fasta database and extract the sequence
	find_flag=0
	while read line; do
		if [[ "$line" != ">"* ]]; then
			continue
		elif [[ "$line" == ">${ID}"* ]]; then
			find_flag=1
			echo "$line" > ${ID}.fasta
			while read next_line; do
				if [[ "$next_line" == ">"* ]]; then
					break
				fi
				echo "$next_line" >> ${ID}.fasta
			done
			break
		fi
	done < $database
	if [[ $find_flag -eq 1 ]]; then
		echo "INFO: ID $ID found in input database ($database)" 
	else
		echo "ERROR: ID $ID does NOT FOUND in database ($database)"
		exit 1
	fi
}

# Generate the Discoba MSA ---------------------------------------------------

# Output directory
output_dir=colabfold_MSA
if [ ! -d "$output_dir" ]; then
	mkdir colabfold_MSA
fi

# Output a3m file
output_a3m=""
for ((i=0; i<${#IDs_array[@]}; i++)); do
	output_a3m="$output_a3m${IDs_array[i]}"
	if (( i < ${#IDs_array[@]} - 1 )); then
		output_a3m="$output_a3m""__vs__"
	fi
done
output_a3m_file=${output_a3m}.a3m


# If the protein was not previously queried
if [ ! -f "./colabfold_MSA/$output_a3m_file" ]; then
	
	[ ! -d fasta_tmp/ ] && mkdir fasta_tmp

	# Recover the sequences. If not, exit the program
	for ID in ${IDs_array[@]}; do
		# This generates protein_ID.fasta in the WD
		get_sequence $database $ID || exit 1
		mv $ID.fasta fasta_tmp/
	done

	# If it is a monomer just use the single sequence
	if [ "${#IDs_array[@]}" == "1" ]; then
		echo "WARNING: Is monomer"
		mv fasta_tmp/${IDs_array[i]}.fasta >> tmp_MSA2/${output_a3m}.fasta
		
	# If it is a complex combine the sequences into a query.fasta
	else
		[ ! -d tmp_MSA2 ] && mkdir tmp_MSA2
		# Add the header
		echo ">$output_a3m" > tmp_MSA2/${output_a3m}.fasta
		# Add the combined sequence
		for ((i=0; i<${#IDs_array[@]}; i++)); do
			grep -v "^>" fasta_tmp/${IDs_array[i]}.fasta >> tmp_MSA2/${output_a3m}.fasta
			if [ "$((i + 1))" -ne "${#IDs_array[@]}" ]; then
				echo ":" >> tmp_MSA2/${output_a3m}.fasta
			fi
		done
	fi

	colabfold_batch --num-recycle 0 --num-models 1 --num-seeds -1 tmp_MSA2/${output_a3m}.fasta tmp_MSA || echo "DONE: IndexError at this point is expected."

	# Move results and remove temporary files
	mv tmp_MSA/$output_a3m_file $output_dir
#	mv tmp_MSA/${output_a3m}_coverage.png $output_plot_dir
#	rm ${output_a3m}.fasta
	rm -r tmp_MSA
	rm -r fasta_tmp/
	rm -r tmp_MSA2/

else
	echo "WARNING: $output_a3m MSA generated beforehand. The search was not performed."
	exit 0
fi


cd $WD
