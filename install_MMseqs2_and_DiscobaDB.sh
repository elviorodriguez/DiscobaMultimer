#!/bin/bash -e

# This script Installs MMseqs2, fetch and build Discoba database.

usage() {
	echo "Before installation, makes sure the following dependencies are met:"
	echo "	-curl                    (sudo apt-get curl)"
	echo "	-miniconda               (see intall_dependencies_aws.txt)"
	echo "	-biopython               (pip install biopython)"
	echo "	-emboss                  (sudo apt-get install emboss)"
	echo "	-hhsuite                 (sudo apt-get install hhsuite)"
	echo "	-RoseTTAFold 2-track     (see intall_dependencies_aws.txt)"
	echo "	-LocalColabFold          (see intall_dependencies_aws.txt)"
	echo " "
	echo "To perform the instalation run the script in the directory to be installed:"
	echo "$0 install"
	echo "This will create a folder called DiscobaMultimer in the working dir."
	echo "Once installed, do not move it."
	echo ""
	echo "Some environmental variables will be created/modified:"
	echo "	-DiscobaMultimerPath: path to the installation folder"
	echo "	-DiscobaDB: path to discoba database"
	echo "	-PATH: mmseqs will be added to PATH"
	exit 1
}

[ "$1" != "install" ] && usage

mkdir DiscobaMultimer
cd DiscobaMultimer

WD=`pwd`

echo "-------------------- Discoba-Multimer installation --------------------"

#Install MMseqs2, for search and alignment
echo "Verifying pre-existing mmseqs2 installation"
mmseqs --help >> /dev/null && touch MMSEQSLOC_READY
if [ ! -f MMSEQSLOC_READY ]; then
	
	if [ -d mmseqs ]; then
		echo "	Damaged installation detected. It will be removed and reinstalled"
		rm -r mmseqs
	else 
		echo "	No MMseqs2 detected"
	fi

	MMSEQS2_DOWNLOAD_PAGE=https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
	MMSEQS2_TAR=`basename $MMSEQS2_DOWNLOAD_PAGE`
	echo "Starting download from $MMSEQS2_DOWNLOAD_PAGE"
	echo "This will take some time..."
	curl $MMSEQS2_DOWNLOAD_PAGE -s -L -o $MMSEQS2_TAR
	echo "Unzipping MMseqs2"
	tar -xf $MMSEQS2_TAR
	rm $MMSEQS2_TAR

	# Add mmseqs to the PATH and autocompletion
	MMSEQS_PATH=${WD}/mmseqs/bin
	echo "" >> ~/.bashrc
	echo "# Add mmseqs to PATH" >> ~/.bashrc
	echo "export PATH=${MMSEQS_PATH}:\$PATH" >> ~/.bashrc
	echo "# Add mmseqs autocompletion" >> ~/.bash_profile
	echo "source ${WD}/mmseqs/util/bash-completion.sh" >> ~/.bash_profile
	
	touch MMSEQSLOC_READY
else
	echo "Local MMseqs2: preexisting installation found. Installation not performed."
fi
echo "Local MMseqs2: READY"

#Download the custom Discoba database and create MMseqs2 search database
echo "Verifying pre-existing Discoba database"
[ -z ${DiscobaDB} ] || touch DISCOBA_READY
if [ ! -f DISCOBA_READY ]; then
	echo "Searching for pre-existing folders"
	if [ -d discoba ]; then
		echo "Damaged installation detected. It will be removed and reinstalled"
		read -p "Do you want to continue? (y/n) " response && [ "$response" == "n" ] && { echo "Program ended."; exit 1; }
		rm -r discoba
	else
		echo "No previous database detected"
	fi

	# Original DiscobaDB:
	#DISCOBA_DOWNLOAD_PAGE=http://wheelerlab.net/discoba.fasta.gz

	# DiscobaDB with TaxIDs in the header for complex prediction:
	DISCOBA_DOWNLOAD_PAGE=https://www.dropbox.com/s/m8vrdijx0slexsn/discoba_TaxID.fasta.gz
	DISCOBA_ZIP=discoba_TaxID.fasta.gz
	DISCOBA_FILE=discoba_TaxID.fasta
	
	mkdir discoba
	cd discoba
		
		# discoba
		echo "DiscobaDB statistics"
		curl http://wheelerlab.net/discobaStats.txt
		echo "Discoba TaxIDs were added to perform MSA pairings"
		echo "Starting DiscobaDB download from $DISCOBA_DOWNLOAD_PAGE"
		echo "This will take some time..."
		curl $DISCOBA_DOWNLOAD_PAGE -s -L -o $DISCOBA_ZIP
		echo "Unzipping DiscobaDB"
		gzip -d $DISCOBA_ZIP
		
		# Create mmseqs searchable database
		echo "Generating MMseqs2 searchable DiscobaDB..."
		../mmseqs/bin/mmseqs createdb discoba_TaxID.fasta discoba
		echo "Creating indexes..."
		../mmseqs/bin/mmseqs createindex discoba tmp4 --remove-tmp-files 1
		
	cd ..
	
	# Add DiscobaDB location as env variable
	echo "# Adds DiscobaDB as env variable" >> ~/.bashrc
	echo "export DiscobaDB=${WD}/discoba/discoba" >> ~/.bashrc
	
	# Create env variable to store the installation path
	echo "# env variable to store the installation folder" >> ~/.bashrc
	echo "export DiscobaMultimerPath=${WD}" >> ~/.bashrc
	
	touch DISCOBA_READY
else
	echo "MMseqs2 DiscobaDB: preexisting installation found. Installation not performed."
fi
echo "MMseqs2 DiscobaDB: READY"

cd $WD

echo "Restart the shell to use the program or simply run these commands:"
echo "source ~/.bashrc"
echo "source ~/.bash_profile"

