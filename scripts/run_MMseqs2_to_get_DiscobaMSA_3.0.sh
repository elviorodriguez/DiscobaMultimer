#!/bin/bash -e

# This script runs MMseqs2 to get Discoba MSA 

# DEPENDENCIES -------------------------------------------------------------
# mmseqs in the PATH
command -v mmseqs >/dev/null 2>&1 || { echo >&2 "ERROR: MMseqs2 not found. Install it using installation script. Aborting."; exit 1; }
# DiscobaDB as env variable (Installation folder)
: "${DiscobaDB:? ERROR: DiscobaDB not found. Install it using installation script. Aborting.}"
: "${DiscobaMultimerPath:? ERROR: DiscobaMultimerPath not found. Install it using installation script. Aborting.}"
DISCOBA_tmp=$DiscobaMultimerPath/discoba/tmp4
# Path to reformat_mmseqs_alignment.py
REFORMAT=$DiscobaMultimerPath/scripts/reformat_mmseq_table_2.0.py
# --------------------------------------------------------------------------

usage() {
	echo "USAGE: $0 <database.fasta> <protein_ID>"
	echo "OPTIONS:"
	echo "	- To be implemented..."
	echo "OUTPUT: stored in ./mmseqs_alignments/<protein_ID> directory"
	echo "	- protein_ID.a3m"
	echo "	- other files produced by MMseqs2"
	echo "NOTES:"
	echo "	- MMseqs2 and DiscobaDB must be installed in advance"
	echo "	  by running install_MMseqs2_and_DiscobaDB.sh"
	exit 1
}

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

# Check input ---------------------------------------------------------------

if [ $# -ne 2 ]; then
	usage
elif [ ! -f $1 ]; then
	echo "ERROR: $1 is not a file"
	usage
elif ! head -n 1 "$1" | grep -q "^>"; then
	echo "ERROR: $1 does not look like a sequence database"
	usage
fi

database=$1
ID=$2
WD=`pwd`

# Generate the Discoba MSA ---------------------------------------------------

# Output directory
output_dir=discoba_mmseqs_alignments/${ID}
output_file=discoba_mmseqs_alignments/${ID}/${ID}.a3m

# If the protein was not previously queried
if [ ! -f "$output_file" ]; then

	# recover the sequence. If not found, exit the program
	get_sequence $database $ID || exit 1

	# create directory to store results
	mkdir -p $output_dir
	rm -f $output_dir/*

	# Move the fasta file with the sequence to the output_dir
	mv ${ID}.fasta $output_dir/query.fasta

	# search DiscobaDB and produce Discoba MSA
	cd $output_dir
		echo "Converting query.fasta to database..."
		mmseqs createdb query.fasta queryDB -v 2
		echo "Searching Discoba database..."
		mmseqs search queryDB $DiscobaDB resultDB $DISCOBA_tmp --remove-tmp-files 1 -v 2
		echo "Aligning hits..."
		mmseqs align queryDB $DiscobaDB resultDB alignDB -a -v 2
		echo "Formatting result..."
		mmseqs convertalis queryDB $DiscobaDB alignDB queryDB.tab --format-output target,qlen,qstart,qend,tstart,tend,tseq,cigar,taln -v 2
		echo "Reformatting mmseqs table to a3m..."
		# Reformat the mmseqs table to a3m
		python $REFORMAT $ID
		# Add query to the start of the a3m file
		head -1 query.fasta > a3m.tmp
		echo `tail -n +2 query.fasta | tr -d '\n'` >> a3m.tmp
		cat ${ID}.a3m >> a3m.tmp
		mv a3m.tmp ${ID}.a3m

		echo "Results in $output_dir"
else
	echo "WARNING: $ID Discoba MSA generated beforehand. The search was not performed"
	exit 0
fi


cd $WD
