#!/bin/bash -e

usage(){
	echo "USAGE: $0 [OPTIONS] <database.fasta> <IDs_table.txt>"
	echo " database.fasta	: protein database (only ID headers)"
	echo " IDs_table.txt	: list of IDs to generate discoba MSAs/AF2 monomers (.tsv)"
	echo "OPTIONS:"
	echo " -m			: computes MSAs"
	echo " -p			: performs MSA plots"
	echo " -a			: runs AlphaFold2-multimer"
	echo "ADVANCED OPTIONS:"
	echo " -s <MIN_MAX>		: minimum and maximum sizes to compute AF2 monomers"
	echo " -c <AF2.conf>		: path to custom AF2 configuration file"
	echo " -i <MSAs_path>		: path to custom MSAs (a3m format with cardinality)"
	echo ""
	echo "Any bugs? ---> rodriguezaraya@ibr-conicet.gov.ar"
	exit 1
}


# DEPENDENCIES -------------------------------------------------------------

GetColabFoldMSA=$DiscobaMultimerPath/scripts/get_ColabFold_MSA.sh
	# USAGE: $GetColabFoldMSA <database.fasta> <protein_ID1>
	# OUTPUT: ./colabfold_MSA/<protein_ID>.a3m"

GetDiscobaMSA=$DiscobaMultimerPath/scripts/run_MMseqs2_to_get_DiscobaMSA_3.0.sh
	# USAGE: $GetDiscobaMSA <database.fasta> <protein_ID>
	# OUTPUT: ./mmseqs_alignments/<protein_ID>/protein_ID.a3m

GENERATE_PLOT=$DiscobaMultimerPath/scripts/plot_msa.py


# SUB-DEPENDENCIES ---------------------------------------------------------

# mmseqs in the PATH
command -v mmseqs >/dev/null 2>&1 || { echo >&2 "ERROR: MMseqs2 not found. Install it using installation script. Aborting."; exit 1; }

# DiscobaDB as env variable (Installation folder)
: "${DiscobaDB:? ERROR: DiscobaDB not found. Install it using installation script. Aborting.}"
: "${DiscobaMultimerPath:? ERROR: DiscobaMultimerPath not found. Install it using installation script. Aborting.}"

# Path to reformat_mmseqs_alignment.py
REFORMAT=$DiscobaMultimerPath/scripts/reformat_mmseq_table.py

# --------------------------------------------------------------------------

# Options
make_MSA=false
make_plot=false
alphafold=false
alphafold_tag=""
import_MSA=false
import_MSA_tag=""
sizes=false
sizes_tag=""
AF2_conf=false
AF2_conf_tag=""
AF2_conf_file=""
# while getopts "mpa:c:i:s:" opt; do
while getopts "mpac:i:s:" opt; do
  case ${opt} in
    m)
      make_MSA=true
      ;;
    p)
      make_plot=true
      ;;
    a)
      alphafold=true
      alphafold_tag="-a"
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

# Set min and max sizes (if provided)
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


# Progress
add_time "STARTING process..."
echo "---------------- Executing Discoba-Monomer ----------------"
echo "Possitional argument:"
echo " - Database: $database_file"
echo " - ID table: $IDs_table_file"
[ "$AF2_conf_file" != "" ] && echo " - AF2 conf: $AF2_conf_file"
echo "OPTIONS:"
echo " - MSA    (-m): $make_MSA"
echo " - Plot   (-p): $make_plot"
echo " - AF2    (-a): $alphafold"
echo " - Import (-i): $import_MSA" ; [ "$import_MSA" == "true" ] && echo "    path: $merged_msa_path"
echo " - Sizes  (-s): $sizes" ; [ "$sizes" == "true" ] && echo "    min: $min_size" && echo "    max: $max_size"



#####################################################################################
#################################### MSA module #####################################
#####################################################################################

# Get Discoba and ColabFold MSAs one at a time -----------------------------
if [ "$make_MSA" == "true" ]; then
	echo ""
	echo "---------------------------------------------------------------------------"

	# Create output folder if it does not exists
	[ ! -d merged_MSA ] && mkdir merged_MSA		# Output path
	add_time "Batch generation of MSA..."
	
	# Read IDs_table.txt one line at a time
	while IFS=$'\t' read -r -a IDs_array; do
		
		# If the line contains more than one ID 
		if [ "${#IDs_array[@]}" -ne 1 ]; then
			echo ""
			add_time "${IDs_array[@]} is not a monomer. Jumping to next ID."
       			continue  # Skip to the next iteration
   		fi

		protein_ID=$IDs_array
		
		echo ""
		add_time "Batch ID (MSA): ${IDs_array[*]}"

		# Get MSAs
		$GetDiscobaMSA $database_file "${IDs_array[@]}"
		$GetColabFoldMSA $database_file "${IDs_array[@]}"
			
		# Merge the MSAs
		add_time "Merging ColabFold and Discoba MSAs..."

		# If protein_ID merged MSA was already computed, do not compute it
		if [ -f ./merged_MSA/${protein_ID}.a3m ]; then
			add_time "WARNING: ./merged_MSA/${protein_ID}.a3m already exists. Merging not performed."
		# If protein_ID merged MSA does not exists, compute it
		else
			# Concatenate ColabFoldMSA with DiscobaMSA monomers
			cat ./colabfold_MSA/${protein_ID}.a3m > ./merged_MSA/${protein_ID}.a3m
			echo "" >> ./merged_MSA/${protein_ID}.a3m
			cat ./discoba_mmseqs_alignments/${protein_ID}/${protein_ID}.a3m >> ./merged_MSA/${protein_ID}.a3m
		fi
		add_time "DONE: output in ./merged_MSA/${protein_ID}.a3m"
		
		# Perform plot
		if [ "$make_plot" = true ]; then
			[ ! -d "./msa_plots" ] && mkdir ./msa_plots
			output_png="${protein_ID}_msa.png"
			if [ ! -f ./msa_plots/$output_png ]; then
				add_time "Generating MSA plot $output_png..."
				python $GENERATE_PLOT ./merged_MSA/${protein_ID}.a3m
				mv $output_png ./msa_plots
				add_time "DONE: output in ./msa_plots/$output_png"
			else
				add_time "WARNING: Plot generated beforehand. Not performed."
			fi
			
		fi
				

	done < "$IDs_table_file"
fi

# Return IFS to default value
IFS=$' \t\n'


#####################################################################################
############################### AlphaFold2-Monomer ##################################
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

	# Create output folder if it does not exists
	[ ! -d ./AF2 ] && mkdir ./AF2/ || add_time "AF2 directory already exist. Continuing."

	# Scan IDs_file one line at a time
	while IFS=$'\t' read -r -a IDs_array; do

		# If the line contains more than one ID
                if [ "${#IDs_array[@]}" -ne 1 ]; then
                        echo ""
                        add_time "${IDs_array[@]} is not a monomer. Jumping to next ID."
                        continue  # Skip to the next iteration
                fi

                protein_ID=$IDs_array

                echo ""
		add_time "Batch ID (AF2): ${protein_ID}"

		# Select were to find MSAs
		if [ "$import_MSA" == "true" ]; then
			input_a3m_file=$merged_msa_path/${protein_ID}.a3m
			add_time "Using precomputed MSA"
			if [ ! -f "$input_a3m_file" ]; then
				add_time "ERROR: $input_a3m_file is missing in custom MSA database"
				usage
			fi
		else
			input_a3m_file=./merged_MSA/${protein_ID}.a3m
		fi

		# dir to store AF2 models
		output_dir_AF2=./AF2/${protein_ID}/

		# Check protein size
		protein_sequence=$(grep -A 1 "^>" $input_a3m_file | tail -n 1)
		protein_length=$(echo -n "$protein_sequence" | wc -c)

		# If size is outside range, ignore it
		if [ "$sizes" == "true" ]; then
			if [ "$protein_length" -lt "$min_size" ]; then
				add_time "WARNING: $protein_ID length is $protein_ID"
				add_time "WARNING: length smaller than $min_size. Ignoring it"
				add_time "WARNING: ignored a3m file stored in ./ignored.txt"
				echo "ignored:AF2	reason:size($protein_length)<min_size($min_size)	msa_file:$input_a3m_file" >> ignored.txt
				continue
			elif [ "$protein_length" -gt "$max_size" ]; then
				add_time "WARNING: $protein_ID length is $protein_length"
				add_time "WARNING: length bigger than $max_size. Ignoring it"
				add_time "WARNING: ignored a3m file stored in ./ignored.txt"
				echo "ignored:AF2	reason:size($combined_L)>max_size($max_size)	msa_file:$input_a3m_file" >> ignored.txt
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
			echo "Options: --num-models 5 --num-recycle 6 --rank plddt --recycle-early-stop-tolerance 0.5 --num-relax 1 --use-gpu-relax --save-all"
			colabfold_batch --num-models 5 --num-recycle 6 --rank plddt --recycle-early-stop-tolerance 0.5 --num-relax 1 --use-gpu-relax --save-all $input_a3m_file ./AF2/$protein_ID
		fi
	done < "$IDs_table_file"
fi

add_time "discoba-monomer-batch FINISHED"
