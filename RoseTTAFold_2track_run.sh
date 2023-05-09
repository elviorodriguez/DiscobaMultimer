#!/bin/bash -e

# This script runs RoseTTAFold 2-track network from RoseTTAFold env

# DEPENDENCIES -------------------------------------------------------
# Path to RoseTTAFold 2-track
rosetta_2track=~/RoseTTAFold/network_2track/predict_msa.py
# Path to remove_paired.sh
script_remove_p=$DiscobaMultimerPath/scripts/remove_paired.sh
# Path to remove_unpaired.sh
script_remove_u=$DiscobaMultimerPath/scripts/remove_unpaired.sh
# Path to plot_coev_RoseTTAFold_2track_run.py
script_plot=$DiscobaMultimerPath/scripts/RoseTTAFold_2track_plot_coev.py
# --------------------------------------------------------------------

usage() {
  echo "Usage: $0 [-p|-u] [-t <int>] -f <path/to/file.a3m>" 1>&2
  echo "Mandatory argument:"
  echo "-f path : path to the a3m file produced with MMseqs2 (ColabFold)"
  echo "	 the name must be formatted like this: ID1__vs__ID2.a3m"
  echo "Optional:"
  echo "-p      : remove the paired part of the a3m file (NOT RECOMENDED)"
  echo "-u      : remove the unpaired part of the a3m file (RECOMENDED)"
  echo "-t int  : top contacts to calculate (default 30)"
  echo "-o	: order by sequence similarity the MSA befor calculations"
  echo "-s      : plot number of sequences in the output heatmap"
  echo "NOTE: if no -p/-u option is passed, the complete a3m file will be used,"
  echo "      i.e. paired+unpaired (NOT RECOMENDED, but better than -p)"
  echo "Output files: (OUTDATED)"
  echo "- path/to/2track_output-<Nseq>/file.npz"
  echo "	Coevolutionary information stored as a numpy array"
  echo "- path/to/2track_output-<Nseq>/file-coevolution.png"
  echo "	Coevolutionary information from file.npz represented as heatmap"
  echo "- path/to/2track_output-<Nseq>/file.contacts"
  echo "	tsv file with pairs of residues predicted to be in contact"
  exit 1
}

# top_contacts by default
top_contacts=30

# Check arguments and options --------------------------------------------------
while getopts "put:sf:o" opt; do
  case ${opt} in
      	p) p_flag=1;;
	u) u_flag=1;;
	t) t_flag=1 ; top_contacts=$OPTARG;;
	f) f_flag=1 ; input_a3m=$OPTARG;;
	s) s_flag=1;;
	o) o_flag=1;;
	\?) usage;;
	*) usage;;
  esac
done

# Check if option -f was passed
if [[ $f_flag -ne 1 ]]; then
	echo "Error: No -f option was passed" ; usage
fi

# Check if both -u and -p options were passed
if [[ $u_flag -eq 1 && $p_flag -eq 1 ]]; then
  echo "Error: -u and -p cannot be passed together" ; usage
fi

# Test if a3m file starts with cardinality 
card_regex="^#[0-9]\{1,5\},[0-9]\{1,5\}[[:space:]][0-9]\{1,2\},[0-9]\{1,2\}$"
if head -1 $input_a3m | grep -q $card_regex; then
	echo "a3m test PASSED: a3m file contains cardinality"
else
	echo "a3m test NOT PASSED: a3m file does not contain cardinality (e.g. #123,523	1,1)"
	usage
fi

# Aux function -----------------------------------------------------------------
# Switchs the positions between the first and the second protein in a paired alignment
# Generates the reversed paired a3m (ID2__vs__ID1.a3m) 
switch_a3m() {

	# Function input
	a3m_file=$1		# Original a3m file
	no_card_a3m=$2		# No cardinality, no unpaired a3m file ($temp_a3m)


	# Protein pair info
	a3m_basename=`basename ${a3m_file} .a3m`
	ID1=${a3m_basename%__vs__*}
	ID2=${a3m_basename#*__vs__}
	L1=$(( $(head -1 $a3m_file | cut -d "#" -f 2 | cut -d "," -f 1) ))
	L2=$(( $(head -1 $a3m_file | cut -d "#" -f 2 | cut -d "," -f 2 | cut -d $'\t' -f 1) ))
	combined_L=$(( L1 + L2 ))
	switched_a3m_basename=${ID2}__vs__${ID1}_paired

	# Debug
#	echo "$a3m_file"
#	echo "$a3m_basename"
#	echo "ID1: $ID1"
#	echo "ID2: $ID2"
#	echo "L1: $L1"
#	echo "L2: $L2"
#	echo "comb_L: $combined_L"

	# Check if dimer
	IDs_count=$(( $(echo "$a3m_basename" | grep -o "__vs__" | wc -l) + 1 ))
	if [ $IDs_count -ne 2 ]; then
		echo "ERROR: Proteins in a3m file not equal to 2"
		echo "ERROR: a3m file $a3m_file contains $IDs_count"
		echo "ERROR: Aborting"
		exit 1
	fi

	# Definition of counting residues
	AA_upper=("A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "X" "-")

	while read line; do
		if [[ "$line" == ">"* ]]; then
			echo "$line" >> $switched_a3m_basename.a3m
		else
			# Iterate over each character of the sequence
			count=0
			for (( i=0; i<${#line}; i++ )); do

				# If possition correspond to L1, make the switch
				if [[ "$count" == "$L1" ]]; then
					prot_1=${line:0:i}
					prot_2=${line:i:$combined_L}

#					printf -- $prot_2 >> $switched_a3m_basename.a3m
#					printf -- $prot_1 >> $switched_a3m_basename.a3m
#					echo >> $switched_a3m_basename.a3m

					echo "${prot_2}${prot_1}" >> $switched_a3m_basename.a3m
					break
				fi

				# Get the current character at the index i
				char=${line:i:1}
				# Add one if it is a gap or Uppercase letter
				if [[ " ${AA_upper[*]} " =~ " $char " ]]; then
#				if [[ "$char" == "-" || "$char" =~ [[:upper:]] ]]; then
					count=$((count + 1))
				fi

			done
		fi
	done < "$no_card_a3m"
}

# Sorts the sequences in the a3m file by decreasing similarity to the queries
sort_a3m(){

	unsorted_a3m=$1
	sorted_a3m=$2

	head -n 2 $unsorted_a3m > temp_first
	rest_seqs_fasta=$(tail -n +3 $unsorted_a3m)

	HEADERS=()
	SEQUENCES=()
	SCORES=()

	while read -r line1 && read -r line2; do

		# Header and sequence
		header=$line1
		sequence=$line2
		echo ">subject" > temp_second
		echo "$sequence" >> temp_second
		
		# Calculate global alignment
		needle -asequence temp_first -bsequence temp_second -gapopen 10.0 -gapextend 0.5 -outfile temp_global &> /dev/null

		# Extract score
		score=$(grep "Score:" temp_global)
		score=${score##*" "}

		HEADERS+=($header)
		SEQUENCES+=($sequence)
		SCORES+=($score)

	done <<< "$rest_seqs_fasta"

	# Combine the arrays into a single file
	paste <(printf '%s\n' "${HEADERS[@]}") \
	      <(printf '%s\n' "${SEQUENCES[@]}") \
	      <(printf '%s\n' "${SCORES[@]}") \
	      > temp.txt

	# Sort the file based on the third column (the scores)
	sort -rn -k3 temp.txt > sorted.txt

	# Extract the sorted arrays
	HEADERS_sorted=($(cut -f1 sorted.txt))
	SEQUENCES_sorted=($(cut -f2 sorted.txt))
	SCORES_sorted=($(cut -f3 sorted.txt))

	# Remove the temporary files
	rm temp.txt sorted.txt temp_first temp_second temp_global

	# Write sorted a3m
	head -n 2 $unsorted_a3m > $sorted_a3m
	for i in "${!SCORES_sorted[@]}"; do
		if [ "${HEADERS_sorted[$i]}" == "^>"* ]; then
			echo "${HEADERS_sorted[$i]}" >> $sorted_a3m
		else
			echo ">${HEADERS_sorted[$i]}" >> $sorted_a3m
		fi
		echo "${SEQUENCES_sorted[$i]}" >> $sorted_a3m
	done
}

# Checks if any sequence has troubles and removes it
format_a3m() {
	
	# function input
	a3m_F=$1
	real_combined_L=$2

	AA_upper=("A" "C" "D" "E" "F" "G" "H" "I" "K" "L" "M" "N" "P" "Q" "R" "S" "T" "V" "W" "Y" "X" "-")
	
	# line1==header line2==sequence 
	while read -r line1 && read -r line2; do
		
		computed_comb_L=0
		for (( i=0; i<${#line2}; i++ )); do
			char="${line2:i:1}"
			if [[ " ${AA_upper[*]} " =~ " $char " ]]; then
				computed_comb_L=$(( computed_comb_L + 1 ))
			fi
		done

		if [ "$computed_comb_L" != "$real_combined_L" ] ; then
			echo "-----------------sequence has different length:"
			echo "Computed: $computed_comb_L"
			echo "Real    : $real_combined_L"
			echo "Sequence: $line1"
			echo "Header  : $line2"
			echo "Removing it"
		else
			echo "$line1" >> temporal_a3m
			echo "$line2" >> temporal_a3m
		fi
	done < "$a3m_F"

	mv temporal_a3m $a3m_F
}



add_time() {
  # Get the current time in the format "YYYY-MM-DD HH:MM:SS"
  current_time=$(date '+%Y-%m-%d %H:%M:%S')
  
  # Add the current time to the beginning of the input string
  output_string="${current_time} $1"
  
  # Echo the output string
  echo "$output_string"
}


#-------------------------------------------------------------------------------

# Usefull variables with names
input_a3m_dirname=`dirname $input_a3m`
input_a3m_basename=`basename $input_a3m`
input_a3m_name=${input_a3m_basename%.a3m}
output_npz=$input_a3m_dirname/${input_a3m_name}.npz
protein1_ID=`echo $input_a3m_name | awk -F '__vs__' '{print $1}'`
protein2_ID=`echo $input_a3m_name | awk -F '__vs__' '{print $2}'`
output_npz_switched=$input_a3m_dirname/${protein2_ID}__vs__${protein1_ID}.npz

# Extract the length of the first (L1) and second (L2) sequences
first_seq_length=`head -1 $input_a3m | cut -d "#" -f 2 | cut -d "," -f 1`
second_seq_length=`head -1 $input_a3m | cut -d "#" -f 2 | cut -d "," -f 2 | cut -d $'\t' -f 1`
first_seq_length=$((first_seq_length))				# convert to number
second_seq_length=$((second_seq_length))			# convert to number
combined_seq_length=$((first_seq_length + second_seq_length))	# convert to number

# Original length of the MSA
MSA_length=`grep -c '^>' $input_a3m`
MSA_length=$((MSA_length))			# convert to number


# Output some data
echo "------------- RoseTTAFold 2-track --------------"
echo "Input a3m file dirname: $input_a3m_dirname"
echo "Input a3m file basename: $input_a3m_basename"
echo "Input a3m file name: $input_a3m_name"
echo "Protein 1 ID: $protein1_ID"
echo "Protein 2 ID: $protein2_ID"
echo "Protein 1 length (L1): $first_seq_length"
echo "Protein 2 length (L2): $second_seq_length"
echo "Original number of sequences: $MSA_length"

#############################################################################
########################## Preprocess the MSA ###############################
#############################################################################

# Check if a3m file reduction was called
if [[ $p_flag -eq 1 ]]; then
	echo "WARNING: -p option was passed. This option is NOT RECOMMENDED."
	echo "WARNING: paired sequences in the a3m file will not be used."

	# Remove paired sequences
	$script_remove_p $input_a3m > /dev/null

	# Location of the unpaired a3m file
	temp_a3m=$input_a3m_dirname/${input_a3m_name}_unpaired.a3m

	# Remove cardinality
	tail -n +2 $temp_a3m > tmp.a3m
	mv tmp.a3m $temp_a3m
	
	remaining_seqs=`grep -c "^>" $temp_a3m`
	echo "Remainig sequences: $remaining_seqs"

elif [[ $u_flag -eq 1 ]]; then
	add_time "INFO: -u option was passed (RECOMMENDED OPTION)"
	add_time "INFO: unpaired sequences in the a3m file will not be used."
	
	# Remove unpaired sequences
	$script_remove_u $input_a3m > /dev/null
	
	# Location of the paired a3m files (stored in current WD)
	temp_a3m=$input_a3m_dirname/${input_a3m_name}_paired.a3m	# Direct a3m
	switched_a3m=${protein2_ID}__vs__${protein1_ID}_paired.a3m	# Reversed a3m

	# Remove cardinality
	tail -n +2 $temp_a3m > tmp.a3m
	mv tmp.a3m $temp_a3m

	remaining_seqs=`grep -c "^>" $temp_a3m`
	echo "Remainig sequences: $remaining_seqs"
	
	# Sometimes the MSA is too big and needs to be reduced using hhblits diversity filter
	if [ "$remaining_seqs" -gt "50000" ]; then
		add_time "Remaining sequences are bigger than 50000"
		add_time "Applying hhfilter to reduce the MSA length"
		echo "hhfilter options: "
		hhfilter -i $temp_a3m -o temp_filtered.a3m -v 2 -diff 1000 -qid 0.0,0.2,0.4,0.6,0.8,1.0 -qsc 0 -id 95
		mv temp_filtered.a3m $temp_a3m
		remaining_seqs=`grep -c "^>" $temp_a3m`
		add_time "Remaining filtered sequences: $remaining_seqs"

	elif [ "$remaining_seqs" -gt "30000" ]; then
		add_time "Remaining sequences are bigger than 30000"
		add_time "Applying hhfilter to reduce the MSA length"
		echo "hhfilter options: "
		hhfilter -i $temp_a3m -o temp_filtered.a3m -v 2 -diff 1500 -qid 0.0,0.2,0.4,0.6,0.8,1.0 -qsc 0 -id 95
		mv temp_filtered.a3m $temp_a3m
		remaining_seqs=`grep -c "^>" $temp_a3m`
		add_time "Remaining filtered sequences: $remaining_seqs"

	elif [ "$remaining_seqs" -gt "10000" ]; then
		add_time "Remaining sequences are bigger than 10000"
		add_time "Applying hhfilter to reduce the MSA length"
		echo "hhfilter options: "
		hhfilter -i $temp_a3m -o temp_filtered.a3m -v 2 -diff 2000 -qid 0.0,0.2,0.4,0.6,0.8,1.0 -qsc 0 -id 95
		mv temp_filtered.a3m $temp_a3m
		remaining_seqs=`grep -c "^>" $temp_a3m`
		add_time "Remaining filtered sequences: $remaining_seqs"
	fi
	
	# Sort by similarity to query (TO DO) before RF2-track calculations	
	if [[ $o_flag -eq 1 ]]; then
		add_time "Reordering sequences by similarity to the query..."
		sort_a3m $temp_a3m temp_sorted_direct
		mv temp_sorted_direct $temp_a3m
		add_time "Reordering paired MSA complete"
	fi	

	# Switch a3m
	add_time "Switching place of proteins..."
	switch_a3m $input_a3m $temp_a3m
	add_time "Switching complete"

else
	echo "WARNING: no option was passed. This is NOT RECOMMENDED."
	echo "WARNING: both paired+unpaired sequences in the a3m file will be used"
	
	# Location of the paired a3m files (stored in current WD)
	temp_a3m=$input_a3m_dirname/${input_a3m_name}_paired.a3m	# Direct a3m
	switched_a3m=${protein2_ID}__vs__${protein1_ID}_paired.a3m	# Reversed a3m

	# Remove cardinality and store it in a temp file
	tail -n +2 $input_a3m > $temp_a3m

	remaining_seqs=`grep -c "^>" $temp_a3m`

	# Sometimes the MSA is too big and needs to be reduced using hhblits diversity filter
	if [ "$remaining_seqs" -gt "50000" ]; then
		add_time "Remaining sequences are bigger than 50000"
		add_time "Applying hhfilter to reduce the MSA length"
		echo "hhfilter options: "
		hhfilter -i $temp_a3m -o temp_filtered.a3m -v 2 -diff 500 -qid 0.0,0.2,0.4,0.6,0.8,1.0 -qsc 0 -id 95
		mv temp_filtered.a3m $temp_a3m
		remaining_seqs=`grep -c "^>" $temp_a3m`
		add_time "Remaining filtered sequences: $remaining_seqs"

	elif [ "$remaining_seqs" -gt "30000" ]; then
		add_time "Remaining sequences are bigger than 30000"
		add_time "Applying hhfilter to reduce the MSA length"
		echo "hhfilter options: "
		hhfilter -i $temp_a3m -o temp_filtered.a3m -v 2 -diff 700 -qid 0.0,0.2,0.4,0.6,0.8,1.0 -qsc 0 -id 95
		mv temp_filtered.a3m $temp_a3m
		remaining_seqs=`grep -c "^>" $temp_a3m`
		add_time "Remaining filtered sequences: $remaining_seqs"

	elif [ "$remaining_seqs" -gt "10000" ]; then
		add_time "Remaining sequences are bigger than 10000"
		add_time "Applying hhfilter to reduce the MSA length"
		echo "hhfilter options: "
		hhfilter -i $temp_a3m -o temp_filtered.a3m -v 2 -diff 900 -qid 0.0,0.2,0.4,0.6,0.8,1.0 -qsc 0 -id 95
		mv temp_filtered.a3m $temp_a3m
		remaining_seqs=`grep -c "^>" $temp_a3m`
		add_time "Remaining filtered sequences: $remaining_seqs"

	fi
		
	
	# Sort by similarity to query (TO DO) before RF2-track calculations	
	if [[ $o_flag -eq 1 ]]; then
		add_time "Reordering sequences by similarity to the query..."
		sort_a3m $temp_a3m temp_sorted_direct
		mv temp_sorted_direct $temp_a3m
		add_time "Reordering paired MSA complete"
	fi	

	# Switch a3m
	add_time "Switching place of proteins..."
	switch_a3m $input_a3m $temp_a3m
	add_time "Switching complete"

fi


#############################################################################
################### Obtain coevolution with RF2-track #######################
#############################################################################

# format a3m files (removes inconsistencies)
add_time "Formating a3m before process..."
format_a3m $temp_a3m $combined_seq_length
format_a3m $switched_a3m $combined_seq_length

add_time "Obtaining co-evolutionary information with GPU..."

# Needed to run RoseTTAFold 2-track
. ~/anaconda3/etc/profile.d/conda.sh	# source the conda.sh script (necessary)
conda init bash > /dev/null
conda activate RoseTTAFold

run_failed=false

# Run the 2-track prediction

# Try with GPU
python $rosetta_2track -msa $temp_a3m -npz $output_npz -L1 $first_seq_length || run_failed=true
python $rosetta_2track -msa $switched_a3m -npz $output_npz_switched -L1 $second_seq_length || run_failed=true

if [ "$run_failed" == "true" ]; then
	# Use CPU
	add_time "GPU failed, trying with CPU..."
	python $rosetta_2track -msa $temp_a3m -npz $output_npz -L1 $first_seq_length --cpu
	python $rosetta_2track -msa $switched_a3m -npz $output_npz_switched -L1 $second_seq_length --cpu
fi

# Check if RoseTTAFold 2-track succeded
if [ "$?" != "0" ]; then
	add_time "RoseTTAFold 2-track failed. Potential cause: no enough VRAM"
else
	add_time "Output files: $output_npz"
	add_time "Output files: $output_npz_switched"
fi

conda deactivate

rm $temp_a3m
rm $switched_a3m


#############################################################################
##################### Produce plot and contact heatmap ######################
#############################################################################

# Convert to integer
remaining_seqs=$((remaining_seqs))

# Plot contact map and 
python $script_plot $output_npz $output_npz_switched $top_contacts $remaining_seqs

# Move everything to output directory
output_dir=./RoseTTAFold_2track_results/${input_a3m_name}
mkdir -p $output_dir
mv $output_npz $output_dir
mv $output_npz_switched $output_dir
mv ./*${input_a3m_name}-coevolution.png $output_dir
mv ./*${input_a3m_name}.metrics $output_dir
mv ./*${input_a3m_name}.contacts $output_dir
mv ./*${input_a3m_name}.chimeraX $output_dir

