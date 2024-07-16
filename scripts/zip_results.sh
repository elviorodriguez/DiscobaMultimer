#!/bin/bash

usage() {
	echo "USAGE:"
	echo " $0 [OPTIONS] <IDs_table.txt> <project_name>"
	echo "OPTIONS:"
	echo " -m : include all MSA types"
	echo " -p : include MSA plots"
	echo " -r : include RF2-track results"
	echo " -a : include AF2-multimer predictions"
	echo "NOTES:"
	echo " - Must be run inside the dir in which discoba-multimer_batch.sh"
	echo "   was executed"
	exit 1
}

usage_error(){
	echo "ERROR: $1"
	usage
}

retrieve_msa=false
retrieve_plot=false
retrieve_rf=false
retrieve_af2=false
while getopts "mpra" opt; do
	case ${opt} in
		m) retrieve_msa=true ;;
		p) retrieve_plot=true ;;
		r) retrieve_rf=true ;;
		a) retrieve_af2=true ;;
		\?) echo "ERROR: invalid option" ; usage ;;
	esac
done
shift $((OPTIND -1))

[ $# -ne 2 ] && usage_error "Two arguments are needed. $# were given."
IDs_table=$1
project_name=$2
[ ! -f $IDs_table ] && usage_error "IDs_list is not a file"
[ -f ${project_name}.zip ] && usage_error "${project_name}.zip already exists"
[ -d ${project_name} ] && usage_error "${project_name} already exists"

while read line; do
	
	# Read IDs from the line
	IFS=$'\t' read -r -a IDs_array <<< "$line"
	paired_name=`printf '%s__vs__' "${IDs_array[@]}" | sed 's/__vs__$//'`
	IDs_number=${#IDs_array[@]}

	echo "$paired_name ----------------------------------------"
	
	mono_discoba_msa_found=false
	paired_discoba_msa_found=false
	colabfold_msa_found=false
	merged_msa_found=false
	if [ "$retrieve_msa" == "true" ]; then
		echo "Retrieving MSA..."
	fi
	
	plot_found=false
	if [ "$retrieve_plot" == "true" ]; then
		echo "Retrieving plot..."
	fi
	
	rf_found=false
	if [ "$retrieve_rf" == "true" ]; then
		echo "Retrieving RF2-track results..."
	fi
	
	af2_found=false
	if [ "$retrieve_af2" == "true" ]; then
		echo "Retrieving AF2-multimer results..."
	fi

done < $IDs_table
echo "Zipping results..."

zip -r $project_name.zip $project_name


echo "Finished"
