#!/bin/bash -e

# DEPENDENCIES -------------------------------------------------------------
GetDiscobaMSA=$DiscobaMultimerPath/scripts/get_Discoba_MSA.sh
GetColabFoldMSA=$DiscobaMultimerPath/scripts/get_ColabFold_MSA.sh
GENERATE_PLOT=$DiscobaMultimerPath/scripts/plot_msa.py
RF2track_run=$DiscobaMultimerPath/scripts/RoseTTAFold_2track_run.sh

# SUB-DEPENDENCIES ---------------------------------------------------------
# Path to run_MMseqs2_to_get_DiscobaMSA.sh
RUN_MMSEQS_DISCOBA=$DiscobaMultimerPath/scripts/run_MMseqs2_to_get_DiscobaMSA_3.0.sh
# Path to perform_pairing.py
PAIRING=$DiscobaMultimerPath/scripts/perform_pairing_general_solution.py
# mmseqs in the PATH
command -v mmseqs >/dev/null 2>&1 || { echo >&2 "ERROR: MMseqs2 not found. Install it using installation script. Aborting."; exit 1; }
# DiscobaDB as env variable (Installation folder)
: "${DiscobaDB:? ERROR: DiscobaDB not found. Install it using installation script. Aborting.}"
: "${DiscobaMultimerPath:? ERROR: DiscobaMultimerPath not found. Install it using installation script. Aborting.}"
# Path to reformat_mmseqs_alignment.py
REFORMAT=$DiscobaMultimerPath/scripts/reformat_mmseq_table.py
# --------------------------------------------------------------------------

WD=`pwd`

usage() {
	echo "USAGE: $0 [OPTIONS] <database.fasta> <IDs_table.txt>"
	echo " database.fasta	: protein database (only ID headers)"
	echo " IDs_table.txt	: list of IDs to generate the paired+unpaired MSAs (.tsv)"
	echo "OPTIONS:"
	echo " -h 			: full help message"
	echo " -m			: computes MSAs"
	echo " -p			: performs MSA plots"
	echo " -r			: runs RoseTTAFold 2-track"
	echo " -a			: runs AlphaFold2-multimer"
	echo " -g <integer>		: GPUs to use. If surpass system GPUs, all available GPUs will be used"
	echo "ADVANCED OPTIONS:"
	echo " *-d			: use default pre-filter (see help)"
	echo " *-f <method_metric>	: use custom pre-filter given by <method_metric> (see help)"
	echo " *-t <value>		: threshold to surpass (see help)"
	echo " -c <AF2.conf>		: path to custom AF2 configuration file"
	echo " -i <MSAs_path>		: path to custom MSAs (a3m format with cardinality)"
	echo " -s <MIN_MAX>		: minimum and maximum sizes to compute with RF and AF2"
	echo "" 
	echo "*NOT IMPLEMENTED"
	echo ""
	echo "Any bugs? ---> rodriguezaraya@ibr-conicet.gov.ar"
	exit 1
}
help_msg() {
	echo "USAGE: $0 [OPTIONS] <database.fasta> <IDs_table.txt>"
	echo " database.fasta	: protein database (only ID headers)"
	echo " IDs_table.txt	: list of IDs to generate the paired+unpaired MSAs"
	echo "		  Each line must contain at least 2 IDs separated by tabs"
	echo "		  E.g.: Q6GZX3	Q6GZX4"
	echo "		  	Q197F8	Q197F7	Q6GZX2"
	echo "		  	Q6GZX2	Q197F8"
	echo "OPTIONS:"
	echo " -h (help)	: Displays help message"
	echo ""
	echo " -m (MSA)	: Creates MSA using Discoba-Multimer"
	echo ""
	echo " -p (plot)	: Creates a folder with MSA plot representations"
	echo "		  Plots will be stored in a separate ./msa_plots/ folder"
	echo ""
	echo " -r (rosetta)	: Perform RoseTTAFold 2-track calculations"
	echo "		  Results in ./RF2-track/ID1__vs__ID2/"
	echo ""
	echo " *-f (pre-filter): Use <method_metric> for prefilter *NOT_IMPLEMENTED"
	echo "		  if not passed, all AF2 pairs will be predicted"
	echo "		  -f <method_metric> (only compatible with -ra)"
	echo "		  available methods"
	echo "			dir  : direct (ID1__vs__ID2)"
	echo "			rev  : reversed (ID2__vs__ID2)"
	echo "			max  : max contact prob. of both dir and rev methods"
	echo "			min  : min contact prob. of both dir and rev methods"
	echo "			mean : mean contact prob. between dir and rev methods"
	echo "		  available metrics"
	echo "		  	max  :	maximum contact probability"
	echo "		  	summ :  summ of contact probabilities"
	echo "			norm :  normalized summ of contact probability"
	echo "		  examples: -f min_max  -v 0.6"
	echo "		            -f mean_summ  -v 20.0"
	echo "		  benchmark studies showed that the best combination is"
	echo "		  ???_??? combined with a threshold of ??? (-d option)"
	echo ""
	echo " *-d (default)	: Use default prefilter method, metric and threshold"
	echo ""
	echo " *-t (threshold) : Threshold to surpass in order to compute AF2 model"
	echo "		  -t <value> (only compatible with -f)"
	echo "		  values ranges (depends on metric)"
	echo "			max  :	[0.0, 1.0]"
	echo "			summ :	[0.0, inf]"
	echo "			norm :	[0.0, 1.0]"
	echo "		  examples: -f min_max  -v 0.6"
	echo "		            -f mean_summ  -v 20.0"
	echo "		  benchmark studies showed that the best combination is"
	echo "		  ???_??? combined with a threshold of ??? (-d option)"
	echo ""
	echo " -a (alphafold)	: Perform AF2 calculations."
	echo "		  Results in ./AF2/ID1__vs__ID2/ folder"
	echo "		  If no AF2.conf file is passed, AF2 will be executed with"
	echo "		  default parameters (5 models, 3 recycles, no_tol, no_relax)."
	echo ""
	echo " -c (AF_conf)	: config file for running AF2 with custom parameters"
	echo "		  -c <AF2.conf> (only compatible with -a option)"
	echo "		  Look at AF2.config sample file in scripts folder for format"
	echo ""
	echo " -i (import)	: imports merged alignments from previous calculations"
	echo "		  -i path/merged_MSA/ storing ID1__vs__ID2.a3m precomputed MSAs"
	echo "		  Incompatible with -m"
	echo ""
	echo " -s (sizes)	: defines a range of combined lengths to be processed by"
	echo "		  RF2-track and AF2-multimer. Bigger or lower sizes will not be"
	echo "		  processed"
	echo "		  -s <min>_<max> (eg 0_1000 or 600_1400)"
	echo ""
	echo "OUTPUTs:"
	echo "(-m):"
	echo "	./discoba_mmseqs_alignments/"
	echo "		content: monomers MSAs"
	echo "		formats: <protein_ID>.a3m"
	echo "	./discoba_paired_unpaired/"
	echo "		content: paired discobaDB alignments"
	echo "		format: <protein_ID1>__vs__<protein_ID2>__vs__[etc].a3m"
	echo "	./colabfold_MSA/"
	echo "		content: paired colabfoldDB paired+unpaired alignments"
	echo "		format: <protein_ID1>__vs__<protein_ID2>__vs__[etc].a3m"
	echo "	./merged_MSA/"
	echo "		content: combined colabfoldDB+discobaDB paired+unpaired MSA"
	echo "		format: <protein_ID1>__vs__<protein_ID2>__vs__[etc].a3m"
	echo "(-p):"
	echo "	./msa_plots/"
	echo "		content: msa paired+unpaired plots representations"
	echo "		format: <protein_ID1>__vs__<protein_ID2>__vs__[etc].png"
	echo "(-r):"
	echo "	./RoseTTAFold_2track_results/"
	echo "		subfolders named ID1__vs__ID2/ with pair results"
	echo "		only heterodimers will be parsed"
	echo "		content: predicted contact maps (matrix and plots)"
	echo "		         predicted contact pairs and probability (.contacts)"
	echo "		         metrics for PPI prediction (.metrics file)"
	echo "(-a):"
	echo "	./AF2/"
	echo "		subfolders named ID1__vs__ID2/ with pair results"
	echo "		content: predicted structures"
	echo "		         PAE, pLDDT and MSA plots"
	echo ""
	echo "NOTES:"
	echo "	- If the search of any of the proteins was already done,"
	echo "	  the program will skip this part and continue with the pairing."
	echo "	  The same will happen for the paired+upaired alignments."
	echo "	- Lines in IDs_table with only one ID will be dropped, along"
	echo "	  with empty ones"
	echo "	- The program can be used for homodimers, heterotrimers/tetramers/etc"
	echo "	  but RF2-track will only be conducted for heterodimers"
	echo "	- It does not work for monomers"
	echo ""
	echo "Doubts, bugs, comments ---> rodriguezaraya@ibr-conicet.gov.ar"
	exit 1
}

# Check input ---------------------------------------------------------

# Options
make_MSA=false
make_plot=false
rosettafold=false
rosettafold_tag=""
alphafold=false
alphafold_tag=""
import_MSA=false
import_MSA_tag=""
sizes=false
sizes_tag=""
AF2_conf=false
AF2_conf_tag=""
AF2_conf_file=""
default_filter=false
custom_filter=false
custom_threshold=false
use_multiple_GPUs=false
while getopts "hmpag:c:ri:s:" opt; do
  case ${opt} in
    h)
      help_msg
      ;;
    m)
      make_MSA=true
      ;;
    p)
      make_plot=true
      ;;
    r)
      rosettafold=true
      rosettafold_tag="-r"
      ;;
    a)
      alphafold=true
      alphafold_tag="-a"
      ;;
    g)
      use_multiple_GPUs=true
      gpus_number=$OPTARG
      if [[ "$gpus_number" == "-"* ]]; then
      	echo "ERROR: seems that -g argument is missing"
	usage
      elif [[ ! $gpus_number =~ ^[0-9]+$ ]]; then
	echo "ERROR: seems that -g argument is not a positive integer"
	usage
      elif [[ $gpus_number -eq 0 ]]; then
	echo "ERROR: -g argument cannot be zero"
	usage
      fi
      ;;
    c)
      AF2_conf_file=$OPTARG
      if [[ "$AF2_conf_file" == "-"* ]]; then
      	echo "ERROR: seems that -c argument is missing"
	usage
      elif [ ! -f "$AF2_conf_file" ]; then
      	echo "ERROR: AF2.conf argument ($AF2_conf_file) is not a file"
	usage
      fi
      AF2_conf_tag="-c $AF2_conf_file"
      ;;
    i)
      import_MSA=true
      merged_msa_path="$OPTARG"
      if [[ "$AF2_conf_file" == "-"* ]]; then
      	echo "ERROR: seems that -i argument is missing"
	usage
      elif [ ! -d $merged_msa_path ]; then
      	echo "ERROR: MSA path (-i argument) do not exists"
	usage
      fi
      import_MSA_tag="-i $merged_msa_path"
      ;;
    s)
      sizes=true
      if [[ "$OPTARG" == "-"* ]]; then
      	echo "ERROR: seems that -s argument is missing"
	usage
      fi
      min_size=${OPTARG%_*}
      max_size=${OPTARG#*_}
      sizes_tag="-s $OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
  esac
done

# Shift the processed options
shift $((OPTIND -1))


# Check positional arguments
if [ $# -ne 2 ]; then
	echo "ERROR: two positional arguments are needed. $# were given"
	usage
fi
if [ ! -f $1 ]; then
	echo "ERROR: database argument $1 is not a file"
	usage
elif [ ! -f $2 ]; then
	echo "ERROR: IDs_table argument $2 is not a file"
	usage
fi

if [ "$import_MSA" == "true" ] && [ "$make_MSA" == "true" ]; then
	echo "ERROR: -i and -m options are incompatible"
	usage
elif [ "$alphafold" == "false" ] && [ "$AF2_conf_file" != "" ]; then
	echo "ERROR: AF2.conf file (-c $AF2_conf_file) given without calling -a"
	usage
elif [ "$make_MSA" == "false" ] && [ "$make_plot" == "true" ]; then
	echo "ERROR: -p only possible if -m was passed"
	echo "If you want to produce plots use plot_msa.py"
	echo ""
	usage
fi

if [ "$sizes" == "true" ]; then
	if  [[ ! "$min_size" =~ ^[0-9]+$ ]]; then
		echo "ERROR: minimum and maximum size format is incorrect"
		echo "	Min size: $min_size"
		echo "	Max size: $max_size"
		usage
	elif [[ ! "$max_size" =~ ^[0-9]+$ ]]; then
		echo "ERROR: minimum and maximum size format is incorrect"
		usage
	elif [ "$min_size" -gt "$max_size" ]; then
		echo "ERROR: minimum size is bigger than maximum size"
		usage
	fi
fi


add_time() {
  # Get the current time in the format "YYYY-MM-DD HH:MM:SS"
  current_time=$(date '+%Y-%m-%d %H:%M:%S')
  
  # Add the current time to the beginning of the input string
  output_string="${current_time} $1"
  
  # Echo the output string
  echo "$output_string"
}

# Assign positional arguments to variables
database_file=$1
IDs_table_file=$2

add_time "STARTING process..."
echo "---------------- Executing Discoba-Multimer ----------------"
echo "Possitional argument:"
echo " - Database: $database_file"
echo " - ID table: $IDs_table_file"
[ "$AF2_conf_file" != "" ] && echo " - AF2 conf: $AF2_conf_file"
echo "OPTIONS:"
echo " - MSA    (-m): $make_MSA"
echo " - Plot   (-p): $make_plot"
echo " - RF2    (-r): $rosettafold"
echo " - AF2    (-a): $alphafold"
echo " - Import (-i): $import_MSA" ; [ "$import_MSA" == "true" ] && echo "    path: $merged_msa_path"
echo " - Sizes  (-s): $sizes" ; [ "$sizes" == "true" ] && echo "    min: $min_size" && echo "    max: $max_size"
echo " - GPUs   (-g): $use_multiple_GPUs"; [ "$use_multiple_GPUs" == "true" ] && echo "    number: $gpus_number"

#####################################################################################
#################################### MSA module #####################################
#####################################################################################

# Get Discoba and ColabFold MSAs one at a time -----------------------------
if [ "$make_MSA" == "true" ]; then
	echo ""
	echo "---------------------------------------------------------------------------"
	[ ! -d merged_MSA ] && mkdir merged_MSA		# Output path
	add_time "Batch generation of MSA..."
	while read line; do
		# Read IDs from the line
		IFS=$'\t' read -r -a IDs_array <<< "$line"
		paired_name=`printf '%s__vs__' "${IDs_array[@]}" | sed 's/__vs__$//'`
		IDs_number=${#IDs_array[@]}
		
		# Return IFS to default value
		IFS=$' \t\n'
		
		echo ""
		add_time "Batch IDs: ${IDs_array[*]}"

		# Check if homooligomer
		is_homooligomer=true
		first_element="${IDs_array[0]}"
		for element in "${IDs_array[@]}"; do
			if [[ "$element" != "$first_element" ]]; then
				is_homooligomer=false
				break
			fi
		done

		if [ $IDs_number -lt 2 ]; then
			add_time "	WARNING: At least 2 IDs are requiered. Ignoring IDs line."
		else
			# Get MSAs
			$GetDiscobaMSA $database_file "${IDs_array[@]}"
			$GetColabFoldMSA $database_file "${IDs_array[@]}"
			
			# Merge the MSAs
			add_time "Merging ColabFold and Discoba MSAs..."
			if [ -f ./merged_MSA/$paired_name.a3m ]; then
				add_time "WARNING: ./merged_MSA/$paired_name.a3m already exists. Merging not performed."
			elif [ "$is_homooligomer" == "true" ]; then
				# Use the discoba monomer MSA
				cat ./colabfold_MSA/$paired_name.a3m > ./merged_MSA/$paired_name.a3m
				echo "" >> ./merged_MSA/$paired_name.a3m
				cat ./discoba_mmseqs_alignments/$first_element/$first_element.a3m >> ./merged_MSA/$paired_name.a3m
			else
				# Use the discoba paired+unpaired MSA
				cat ./colabfold_MSA/$paired_name.a3m > ./merged_MSA/$paired_name.a3m
				echo "" >> ./merged_MSA/$paired_name.a3m
				tail -n +2 ./discoba_paired_unpaired/$paired_name.a3m >> ./merged_MSA/$paired_name.a3m
			fi
			add_time "DONE: output in ./merged_MSA/$paired_name.a3m"
			
			# Perform plot
			if [ "$make_plot" = true ]; then
				[ ! -d "./msa_plots" ] && mkdir ./msa_plots
				output_png="${paired_name}_msa.png"
				if [ ! -f ./msa_plots/$output_png ]; then
					add_time "Generating MSA plot $output_png..."
					python $GENERATE_PLOT ./merged_MSA/$paired_name.a3m
					mv $output_png ./msa_plots
					add_time "DONE: output in ./msa_plots/$output_png"
				else
					add_time "WARNING: Plot generated beforehand. Not performed."
				fi
				
			fi
					
		fi
	done < "$IDs_table_file"
fi

# Return IFS to default value
IFS=$' \t\n'

#####################################################################################
#################### Parallel GPU processing module (recursion) #####################
#####################################################################################

if [ "$use_multiple_GPUs" == "true" ]; then

	[ ! -d ./reports ] && mkdir ./reports/
	[ ! -d ./split_IDs ] && mkdir ./split_IDs/

	# Check system GPUs number availability
	max_gpus=$(nvidia-smi --query-gpu=name --format=csv,noheader | wc -l)
	if [[ "$gpus_number" -gt "$max_gpus" ]]; then
		add_time "WARNING: GPU number provided ($gpus_number) is bigger than available GPUs ($max_gpus)"
		add_time "WARNING: GPU number to use will be set to $max_gpus"
		gpus_number=$((max_gpus))
	fi

	# Check that complexes to predict do not surpass the number of GPUs
	lines_number_in_IDs_file=$(cat $IDs_table_file | wc -l)
	if [[ "$gpus_number" -gt "$lines_number_in_IDs_file" ]]; then
		add_time "WARNING: GPU number provided ($gpus_number) is bigger than lines in IDs file ($lines_number_in_IDs_file)"
		add_time "WARNING: GPU number to use will be set to $lines_number_in_IDs_file"
		gpus_number=$((lines_number_in_IDs_file))
	fi

	# Split IDs_file for parallel process
	split --numeric-suffixes=1 -n l/$gpus_number $IDs_table_file ./split_IDs/$(basename $IDs_table_file)_temp_

	# Perform parallel processing in the background	
	date_to_add=$(date '+%Y%m%d_%H%M%S')
	for ((i=1; i<=$gpus_number; i++)); do
		add_time "Starting to process GPU number: $i"
		GPU_index=$((i-1))
		formatted_index=$(printf "%02d" "$i")
		IDs_table_file_i=./split_IDs/${IDs_table_file}_temp_${formatted_index}
		report_file_i=./reports/report_${date_to_add}_${formatted_index}.log
		CUDA_VISIBLE_DEVICES=$GPU_index $0 $rosettafold_tag $alphafold_tag $sizes_tag $AF2_conf_tag $import_MSA_tag $database_file $IDs_table_file_i > $report_file_i 2>&1 &
		add_time "PID GPU number $i: $!"
	done
	
	add_time "Parallel processes running on the background..."

	exit
fi

#####################################################################################
############################## RoseTTAFold 2-track ##################################
#####################################################################################

sort_a3m(){

	unsorted_a3m=$1
	sorted_a3m=$2
	basename_unsorted_a3m=$(basename $unsorted_a3m .a3m)

	head -n 2 $unsorted_a3m > ${basename_unsorted_a3m}temp_first
	rest_seqs_fasta=$(tail -n +3 $unsorted_a3m)

	HEADERS=()
	SEQUENCES=()
	SCORES=()

	while read -r line1 && read -r line2; do

		# Header and sequence
		header=$line1
		sequence=$line2
		echo "$header" > ${basename_unsorted_a3m}temp_second
		echo "$sequence" >> ${basename_unsorted_a3m}temp_second

		# Calculate global alignment
		needle -asequence ${basename_unsorted_a3m}temp_first -bsequence ${basename_unsorted_a3m}temp_second -gapopen 10.0 -gapextend 0.5 -outfile ${basename_unsorted_a3m}temp_global &> /dev/null

		# Extract score
		score=$(grep "Score:" ${basename_unsorted_a3m}temp_global)
		score=${score##*" "}

		HEADERS+=($header)
		SEQUENCES+=($sequence)
		SCORES+=($score)

	done <<< "$rest_seqs_fasta"

	# Combine the arrays into a single file
	paste <(printf '%s\n' "${HEADERS[@]}") \
	      <(printf '%s\n' "${SEQUENCES[@]}") \
	      <(printf '%s\n' "${SCORES[@]}") \
	      > ${basename_unsorted_a3m}temp.txt

	# Sort the file based on the third column (the scores)
	sort -rn -k3 ${basename_unsorted_a3m}temp.txt > ${basename_unsorted_a3m}sorted.txt

	# Extract the sorted arrays
	HEADERS_sorted=($(cut -f1 ${basename_unsorted_a3m}sorted.txt))
	SEQUENCES_sorted=($(cut -f2 ${basename_unsorted_a3m}sorted.txt))
	SCORES_sorted=($(cut -f3 ${basename_unsorted_a3m}sorted.txt))

	# Remove the temporary files
	rm ${basename_unsorted_a3m}temp.txt ${basename_unsorted_a3m}sorted.txt ${basename_unsorted_a3m}temp_first ${basename_unsorted_a3m}temp_second ${basename_unsorted_a3m}temp_global

	# Write sorted a3m
	head -n 2 $unsorted_a3m > $sorted_a3m
	for i in "${!SCORES_sorted[@]}"; do
		if [ "${HEADERS_sorted[$i]}" == "^>"* ]; then
			echo "${HEADERS_sorted[$i]}" >> $sorted_a3m
		else
			echo ">${HEADERS_sorted[$i]}" >> $sorted_a3m
		fi
		echo "${SEQUENCES_sorted[$i]}" >> $sorted_a3m
#		echo "${SCORES_sorted[$i]}"
	done
}



if [ "$rosettafold" == "true" ]; then
	echo ""
	echo "---------------------------------------------------------------------------"
	add_time "Batch generation of RoseTTAFold 2-track contact maps..."

	# Scan IDs_file one line at a time
	while read line; do
		

		# Read IDs from the line
		IFS=$'\t' read -r -a IDs_array <<< "$line"
		paired_name=`printf '%s__vs__' "${IDs_array[@]}" | sed 's/__vs__$//'`
		IDs_number=${#IDs_array[@]}
		
		# Progress
		echo ""
		add_time "Batch IDs: ${IDs_array[*]}"

		# Select were to find MSAs
		if [ "$import_MSA" == "true" ]; then
			input_a3m_file=$merged_msa_path/$paired_name.a3m
			if [ ! -f "$input_a3m_file" ]; then
				echo "ERROR: $input_a3m_file is missing in custom MSA database"
				usage
			fi
		else
			input_a3m_file=./merged_MSA/$paired_name.a3m
		fi
		
		# Check if prediction was already performed
		if [ -d "./RoseTTAFold_2track_results/$paired_name/" ]; then
			add_time "RF2-track prediction for $paired_name already performed"
			add_time "Ignoring line"
			continue	
		fi

		# Check if homooligomer
		is_homooligomer=true
		first_element="${IDs_array[0]}"
		for element in "${IDs_array[@]}"; do
			if [[ "$element" != "$first_element" ]]; then
				is_homooligomer=false
				break
			fi
		done
		
		# If it is a homooligomer: ignore it (homooligomers gives internal contact map)
		if [ "$is_homooligomer" == "true" ]; then
			add_time "WARNING: IDs corresponds to HOMO-oligomer"
			add_time "WARNING: Only HETERO-dimers can be parsed"
			add_time "WARNING: Ignoring IDs line. Saved in ./ignored.txt"
			echo "ignored:RF-2track	reason:is_homooligomer	msa_file:$input_a3m_file" >> ignored.txt

		# If it has more than 2 proteins: ignore it (only pairwise comparisons)
		elif [ $IDs_number -ne 2 ]; then
			add_time "WARNING: Only dimers can be parsed"
			add_time "WARNING: Instead, $IDs_number IDs are in the line"
			add_time "WARNING: Ignoring IDs line. Saved in ./ignored.txt"
			echo "ignored:RF-2track	reason:proteins>2	msa_file:$input_a3m_file" >> ignored.txt

		# If it is a heterodimer: check size and compute coevolution
		else
			# Check combined size
			L1=`head -1 $input_a3m_file | cut -d "#" -f 2 | cut -d "," -f 1`
			L2=`head -1 $input_a3m_file | cut -d "#" -f 2 | cut -d "," -f 2 | cut -d $'\t' -f 1`
			combined_L=$(( L1 + L2 ))

			# If size is outside range, ignore it
			if [ "$sizes" == "true" ]; then
				if [ "$combined_L" -lt "$min_size" ]; then
					add_time "WARNING: $paired_name is smaller than $min_size"
					add_time "WARNING: ignoring it. Saved in ./ignored.txt"
					echo "ignored:RF-2track	reason:combined_size($combined_L)<min_size($min_size)	msa_file:$input_a3m_file" >> ignored.txt
					continue
				elif [ "$combined_L" -gt "$max_size" ]; then
					add_time "WARNING: $paired_name is bigger than $max_size"
					add_time "WARNING: ignoring it. Saved in ./ignored.txt"
					echo "ignored:RF-2track	reason:combined_size($combined_L)>max_size($max_size)	msa_file:$input_a3m_file" >> ignored.txt
					continue
				fi
			fi
	
			# Extract coevolution with RF 2-track
			add_time "ROSETTAFOLD 2-track: $input_a3m_file"
			[ -f $input_a3m_file ] && echo "File exists"
			$RF2track_run -us -f $input_a3m_file
		fi
	done < "$IDs_table_file"

fi


#####################################################################################
############################## AlphaFold2-multimer ##################################
#####################################################################################

# Converts AF2.config file
parse_cfg_file() {
	cfg_file=$1
	configurations=()
	while read line; do
		# Ignore comments
		if [[ "$line" == "#"* ]]; then
			continue
		else
			configurations+=($line)
		fi
	done < $cfg_file
	echo ${configurations[@]}
}

if [ "$alphafold" == "true" ]; then
	echo ""
	echo "---------------------------------------------------------------------------"
	add_time "Batch generation of AF2 models..."
	[ ! -d ./AF2 ] && mkdir ./AF2/ || add_time "AF2 directory already exist. Continuing."

	# Scan IDs_file one line at a time
	while read line; do
		
		# Read IDs from the line
		IFS=$'\t' read -r -a IDs_array <<< "$line"
		paired_name=`printf '%s__vs__' "${IDs_array[@]}" | sed 's/__vs__$//'`
		IDs_number=${#IDs_array[@]}
		
		# Progress
		echo ""
		add_time "Batch IDs: ${IDs_array[*]}"

		# Select were to find MSAs
		if [ "$import_MSA" == "true" ]; then
			input_a3m_file=$merged_msa_path/$paired_name.a3m
			add_time "Using precomputed MSA"
			if [ ! -f "$input_a3m_file" ]; then
				add_time "ERROR: $input_a3m_file is missing in custom MSA database"
				usage
			fi
		else
			input_a3m_file=./merged_MSA/$paired_name.a3m
		fi
		
		# dir to store AF2 models
		output_dir_AF2=./AF2/$paired_name/
		
		# Check if homooligomer
		is_homooligomer=true
		first_element="${IDs_array[0]}"
		for element in "${IDs_array[@]}"; do
			if [[ "$element" != "$first_element" ]]; then
				is_homooligomer=false
				break
			fi
		done

		# Check combined size
		string=$(head -1 $input_a3m_file)
		list1=$(echo "$string" | awk -F'\t' '{print $1}')
		list2=$(echo "$string" | awk -F'\t' '{print $2}')
		list1_arr=($(echo "${list1//\#/}" | tr ',' ' '))
		list2_arr=($(echo "$list2" | tr ',' ' '))
		combined_L=0
		for i in "${!list1_arr[@]}"; do
			product=$((list1_arr[i] * list2_arr[i]))
			combined_L=$((combined_L + product))
		done

		# If size is outside range, ignore it
		if [ "$sizes" == "true" ]; then
			if [ "$combined_L" -lt "$min_size" ]; then
				add_time "WARNING: $paired_name combined length is $combined_L"
				add_time "WARNING: length smaller than $min_size. Ignoring it"
				add_time "WARNING: ignored a3m file stored in ./ignored.txt"
				echo "ignored:AF2	reason:combined_size($combined_L)<min_size($min_size)	msa_file:$input_a3m_file" >> ignored.txt
				continue
			elif [ "$combined_L" -gt "$max_size" ]; then
				add_time "WARNING: $paired_name combined length is $combined_L"
				add_time "WARNING: length bigger than $max_size. Ignoring it"
				add_time "WARNING: ignored a3m file stored in ./ignored.txt"
				echo "ignored:AF2	reason:combined_size($combined_L)>max_size($max_size)	msa_file:$input_a3m_file" >> ignored.txt
				continue
			fi
		fi

		# Return IFS to default value
		IFS=$' \t\n'


		# If AF2.conf file has been passed
		if [ "$AF2_conf_file" != "" ]; then
			
			# Parse options
			options=$(parse_cfg_file $AF2_conf_file)
			
			# Use AF2.conf file as configurations
			echo "AF2 with config file: $AF2_conf_file"
			echo "Options: $options"
			colabfold_batch $options $input_a3m_file $output_dir_AF2
			
		else 
			# Use default configuration
			echo "AF2 without config"
			echo "Options: --num-models 5 --num-recycle 6 --rank iptm --stop-at-score 80 --recycle-early-stop-tolerance 1.5 --num-relax 1 --use-gpu-relax"
			colabfold_batch --num-models 5 --num-recycle 6 --rank iptm --stop-at-score 80 --recycle-early-stop-tolerance 1.5 --num-relax 1 --use-gpu-relax $input_a3m_file ./AF2/$paired_name
		fi
	done < "$IDs_table_file"
fi

add_time "discoba-multimer-batch FINISHED"
