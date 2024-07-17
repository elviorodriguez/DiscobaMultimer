#!/bin/bash -e

# This script generates a paired+unpaired a3m file from using DiscobaDB

# DEPENDENCIES -------------------------------------------------------------
# Path to run_MMseqs2_to_get_DiscobaMSA.sh
RUN_MMSEQS_DISCOBA=$DiscobaMultimerPath/scripts/run_MMseqs2_to_get_DiscobaMSA_3.0.sh
# Path to perform_pairing.py
PAIRING=$DiscobaMultimerPath/scripts/perform_pairing_general_solution.py
GREEDY_PAIRING=$DiscobaMultimerPath/scripts/perform_greedy_pairing.py

# SUB-DEPENDENCIES ---------------------------------------------------------
# mmseqs in the PATH
command -v mmseqs >/dev/null 2>&1 || { echo >&2 "ERROR: MMseqs2 not found. Install it using installation script. Aborting."; exit 1; }
# DiscobaDB as env variable (Installation folder)
: "${DiscobaDB:? ERROR: DiscobaDB not found. Install it using installation script. Aborting.}"
: "${DiscobaMultimerPath:? ERROR: DiscobaMultimerPath not found. Install it using installation script. Aborting.}"
# Path to reformat_mmseqs_alignment.py
REFORMAT=$DiscobaMultimerPath/scripts/reformat_mmseq_table.py
# --------------------------------------------------------------------------

usage() {
	echo "USAGE: $0 [-greedy] <database> <protein_ID_1> <protein_ID_2> [<protein_ID_n>]"
	echo "  -greedy			: Use greedy pairing method"
	echo "  database		: protein database in fasta containing the queries"
	echo "			  (typically the proteome of an organism)"
	echo "  protein_ID_1	: ID of the first protein. Must be contained in the DB"
	echo "  protein_ID_2	: ID of the second protein. Must be contained in the DB"
	echo "  protein_ID_n	: ID of the Nth protein (optional, as many as needed)"
	echo "OUTPUT:"
	echo "  - a3m file named <protein_ID_1>__vs__<protein_ID_2>__vs__[etc].a3m"
	echo "		Has cardinality as first line (e.g. #421,512	1,1)"
	echo "	- output will be stored in discoba_paired_unpaired/ directory"
	echo "NOTES:"
	echo "	- If the search of any of the proteins was already done, the program"
	echo "	  will skip this part and continue with the pairing."
	exit 1
}

# Check if -greedy option is passed
greedy_mode=false
if [ "$1" == "-greedy" ]; then
	greedy_mode=true
	shift
fi

# Check input ---------------------------------------------------------
if [ $# -lt 3 ]; then
	echo "ERROR: 3 or more arguments are needed"
	usage
elif [ ! -f $1 ]; then
	echo "ERROR: the database argument $1 is not a file"
	usage
fi

# Assign positional arguments to variables
database=$1
IDs_number=$(($# - 1))
IDs_array=()

echo "Number of IDs: $IDs_number"
# Assign each ID to a variable
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

# Generate the monomer a3m files with run_MMseqs2_to_get_DiscobaMSA.sh --------

# Output files locations
a3m_files=()

echo "Generating isolated Discoba alignments"
for i in $(seq 1 $IDs_number); do
	ID_num="ID_$i"
	ID_value=${!ID_num}
	echo "Parsing $ID_num: $ID_value..."
	$RUN_MMSEQS_DISCOBA $database $ID_value

	# Append file location
	a3m_files+=(discoba_mmseqs_alignments/${ID_value}/${ID_value}.a3m)
done

# Generate paired+unpaired a3m file -----------------------------------

# Output a3m file
printf -v output_a3m '%s__vs__' "${IDs_array[@]}"
output_a3m=${output_a3m%"__vs__"}.a3m

# Output directory to store the a3m MSAs
output_dir=discoba_paired_unpaired

# If the output directory does not exist, create it
if [ ! -d $output_dir ]; then
	mkdir $output_dir
fi

# Check if paired+unpaired MSA was not already created
if [ ! -f $output_dir/$output_a3m ]; then
	echo "Generating paired+unpaired alignment..."

	# Choose the pairing method based on the greedy_mode flag
	if [ "$greedy_mode" = true ]; then
		# perform greedy pairing
		python3 $GREEDY_PAIRING $output_a3m ${a3m_files[@]}
	else
		# perform the standard pairing
		python3 $PAIRING $output_a3m ${a3m_files[@]}
	fi

	# move the results to the output dir
	mv $output_a3m $output_dir

	echo "Paired+unpaired a3m MSA generated successfully"
else
	echo "WARNING: $output_a3m paired+unpaired MSA generated beforehand"
	echo "WARNING: The alignment was not performed"
	exit 0
fi

echo "Output file: $output_dir/$output_a3m"
