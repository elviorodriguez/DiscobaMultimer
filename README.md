# DiscobaMultimer
DiscobaMultimer is designed to systematically perform high performance computing (HPC) of complex structure predictions with AlphaFold2-multimer on Discoba species, using LocalColabFold as AF2 implementation.

TaxIDs for Discoba species where added to Discoba database (DiscobaDB) to allow the generation of paired multiple sequence alignments (DiscobaMSA) and combining it with the corresponding ColabFoldMSA generated with ColabFold MSA server. It also contains an implementation to run monomer structure predictions and one to just retrive paired+unpaired DiscobaMSAs (without protein structure predictions).

  - For MSA generation only (no GPU required), take a look at this Colab notebook:
  *To implement...

  - For complex structure predictions (GPU required), take a look at this notebook:
  *To implement...

  - For monomeric  structure predictions (GPU required), take a look at this notebook:
  *To implement...

For non Discoba species, it works just as fine as ColabFold, but gives a systematic way to order AF2-multimer predictions.

To reproduce the data of the original article (see **"References"** section), take a look at the `README.md` file inside `acetylation_related_complexes` directory.

## Local installation (Ubuntu)
Before everything, update apt-get and make sure the dependencies are met:
  - Make sure curl, git and wget are installed:
    ```bash
    # Update apt-get
    sudo apt-get -y update

    # Install curl, git and wget
    sudo apt-get -y install curl git wget
    ```
  - Make sure you have LocalColabFold installed:
    ```bash
    # This must give no errors
    colabfold_batch -h
    ```
    If you get `colabfold_batch: command not found` error, install it following the instructions here: https://github.com/YoshitakaMo/localcolabfold.
  - Make sure you have Biopython package:
    ```bash
    pip install biopython    
    ```
  - Install emboss and hhsuite:
    ```bash
    sudo apt-get -y install emboss hhsuite
    ```

Once dependencies are met, clone the repo and perform the installation:

```bash
# Clone the repo and cd to DiscobaMultimer
git clone https://github.com/elviorodriguez/DiscobaMultimer.git
cd DiscobaMultimer

# Install DiscobaMultimer
bash install/install_MMseqs2_and_DiscobaDB.sh install
```

This will take a few minutes, as MMseqs2 is installed and DiscobaDB is built. In general, when you installed LocalColabFold, mmseqs is added to the path. In this case, MMseqs2 installation is skipped. After installation is completed, restart the shell (or `cd ~/` and then `source ~/.bashrc`).

Check the installation with:
```bash
# This will lead to an error and show the usage message for complex prediction
discoba_multimer_batch


# This will lead to an error and show the usage message for monomer prediction
discoba_monomer_batch
```

## Local installation (Other systems)
The program was successfully tested on Ubuntu 20.04, 22.04, 23.10 and some AWS Linux AMIs (_e.g._, Deep Learning AMI GPU TensorFlow 2.12.0 (Ubuntu 20.04) 20230324). I do not know if it is able to run on Windows subsystem for Linux or MAC. In theory, if you can satisfy all the dependencies, their executables are in the $PATH, and you have bash as interpreter, it must work. Let me know if you try.

## How does it work?
DiscobaMultimer requieres 2 positional arguments, a file with fasta sequences `database.fasta` (in general will be the target species proteome) and a file with combinations of **IDs separated with TABS** `IDs_table.txt` to search in the database:

```bash
discoba_multimer_batch [OPTIONS] <database.fasta> <IDs_table.txt>
discoba_monomer_batch [OPTIONS] <database.fasta> <IDs_table.txt>
```

The options will manage what you want to retrieve (MSA, MSA plot, AF2 predictions, etc) and the configuration of each run.

The fasta sequences on `database.fasta` must contain headers **ONLY with the ID of each sequence**. For example:

```bash
>Tb927.1.650
MDEIKWPLRLLQILRNHLEQVDHPHEPRTLPIAAKCGQGGVLAFLTDVDDDVRRRTGASE
CITACCRGLVLDMIEQCCTLLFLSSKERRVINMAAVQREKRVVARGVGKRRRGDQGETTE
ETAIGETDNSQRKCAAVLPLEYLLRLFVALPSLLAHYDKLGGSTMPAAFKQPLWNYVNAV
LDIMKDMTFLDPSEYIPLR
>Tb927.8.5320
MPLGTDGTPQRQQQQRKSTVDSAVALIESDSELMALFFRHISACPPHGPFKHYSALTFTE
RVRHWFPTHAPARAAEMGDAITPAVVLQLMEKYYDTKHLRSWPLLYVPAQIHLKDVDALI
EKQESLKPERATD
>Tb927.1.3400
MRHVDVLQSLMVRLVDPEEEEHDWLLISSLFYSLTKLMLPATRVKELYYHYLQVHGSRPA
MEAACKDFSNKLSEKLRQRIESTTTTTTATTLAGIHRSSVSGTSAAVGVLDVTKVSAGAA
GGEPGIPVGVLSGSPTVATSVWDESGFCGREAAFRKLTSEQKARWHNLRRELISPDILID
LVQQFTLVDGAAVFLAPVVPSTVMEFNGVRSGPYVTVIHQPLSLMCVKRRVLAARRDYEL
HKHQSGAYLPTGQNNVGVGSRKRERGNGGVTSRASAVSQPPHFSIATNSGEKGNVIRTLQ
ELEQAVWHITANCVMFNAPESYYPYVARKFALACVAIIDDYCTQRIAGV
```

For `IDs_table.txt`, each line will correspond to a prediction. For example:

```
Tb927.1.650  Tb927.8.5320
Tb927.1.650  Tb927.8.5320  Tb927.1.3400
Tb927.1.650  Tb927.1.650
```

This means that DiscobaMultimer will generate predictions for an heterodimer (first line), an heterotrimer (second line) and a homodimer (third line). Do not separate the the IDs with spaces or commas, only **TABS** will work.

NOTE: If you want to modify the headers of each sequence to other names, you can. Just do not add spaces or "-" to the new names. Also, make the corresponding change in the IDs file. If you work directly with IDs, you can change them later to the corresponding proteins symbols you want using the utility script `DiscobaMultimer/utils/annotate_AF2_models.sh`. Just call it with no arguments to see usage.


## Running batches of complex structure predictions
Let's use `database.fasta` and `IDs_table.txt` from before as example. First, create a directory for the project to store databases and IDs table files and `cd` to it. As we want to run a batch of complexes, we need to use `discoba_multimer_batch` alias. The most basic call would be as follows:

```
# Run batch of complex structures predictions with defauls
discoba_multimer_batch -ma database.fasta IDs_table.txt 2>&1 | tee report.log
```

The `-m` tag tells DiscobaMultimer to generate DiscobaMSAs for each ID line, and `-a` tells that AF2 predictions must be performed on those DiscobaMSAs. It is always recommended to redirect both stdout and stderr to a log file to keep track the progress and check for any errors (_e.g._: using `2>&1 | tee report.log`). DiscobaMultimer will first generate all the MSAs and, once completed, it will start to produce the AF2-multimer predictions.

The resulting filesystem will be as follows:
```bash
.
├── AF2                                           <------- AlphaFold2-multimer predictions directories
│   ├── Tb927.1.650__vs__Tb927.1.650
│   ├── Tb927.1.650__vs__Tb927.8.5320
│   └── Tb927.1.650__vs__Tb927.8.5320__vs__Tb927.1.3400
├── colabfold_MSA                                 <------- ColabFoldMSAs produced by ColabFold MSA server
│   ├── Tb927.1.650__vs__Tb927.1.650.a3m
│   ├── Tb927.1.650__vs__Tb927.8.5320.a3m
│   └── Tb927.1.650__vs__Tb927.8.5320__vs__Tb927.1.3400.a3m
├── database.fasta                                <------- Original FASTA database (proteome)
├── discoba_mmseqs_alignments                     <------- DiscobaDB monomer MSAs (MMseqs2 output)
│   ├── Tb927.1.3400
│   ├── Tb927.1.650
│   └── Tb927.8.5320
├── discoba_paired_unpaired                       <------- DiscobaMSA (A3M files)
│   ├── Tb927.1.650__vs__Tb927.1.650.a3m
│   ├── Tb927.1.650__vs__Tb927.8.5320.a3m
│   └── Tb927.1.650__vs__Tb927.8.5320__vs__Tb927.1.3400.a3m
├── IDs_table.txt                                 <------- Original IDs table
├── merged_MSA                                    <------- ColabFoldMSAs+DiscobaMSAs (A3M files)
│   ├── Tb927.1.650__vs__Tb927.1.650.a3m
│   ├── Tb927.1.650__vs__Tb927.8.5320.a3m
│   └── Tb927.1.650__vs__Tb927.8.5320__vs__Tb927.1.3400.a3m
└── report.log                                    <------- Log file
```

By default, discoba_multimer_batch runs `colabfold_batch` using the following options:

```
# discoba_multimer_batch default options
--num-models 5
--num-recycle 6
--rank iptm
--stop-at-score 80
--recycle-early-stop-tolerance 1.5
--num-relax 1
--use-gpu-relax
```
If you want to use a custom AF2 configuration, see **"Using custom AF2.config file"** section.

NOTE: All the filesystem is based on the `ID1 + "__vs__" + ID2 + "__vs__" + ...` notation. If you modify this filesystem nomenclature and then you run additional DiscobaMultimer pipelines with it, your pipeline will not work.

## Running batches of DiscobaMSA only predictions (without AF2 predictions)
Following with the same example as before, it is as easy as removing the `-a` flag.

```
# Run batch of complex structures predictions with defauls
discoba_multimer_batch -m database.fasta IDs_table.txt 2>&1 | tee report2.log
```

If you run the prediction in the same folder as before, DiscobaMultimer will know that the MSA predictions performed were already performed and they will be skipped. You can check it by looking at the log file:

```
cat report2.log
```

The resulting ColabFoldMSAs, DiscobaMSAs, and ColabFoldMSA+DiscobaMSA will be located in `colabfold_MSA`, `discoba_paired_unpaired` and `merged_MSA`, respectively. Additionally, Discoba MSAs for individual proteins (monomers) will be located in `discoba_mmseqs_alignments`.

NOTE: This is useful to save resources. For example, when you run DiscobaMultimer on AWS or other cloud computing service, you can first run the MSA section using only the `-m` flag on a low price instance (_e.g._,without GPU). After the MSAs are generated, you can switch to a higher capacity instance (with multiple GPUs) and run the pipeline again, but this time with both flags: `-ma`. As MSAs are already in the project directory, they will not be computed, and it will jump directly to AF2 section. This will save you a lot of money. For parallel processing using multiple GPUs, see  **"Using múltiple GPUs for parallel computing"** section.

## Running batches of monomeric structure predictions
If you are interested only in monomeric structures, you need to use the `discoba_monomer_batch` alias. 

```bash
# Run batch of monomeric structure predictions with defauls
discoba_monomer_batch -ma database.fasta IDs_table.txt 2>&1 | tee report3.log
```

The flags usage logic is the same as with complexes, but this time, the `IDs_table.txt` file must contain a single ID on each line. For example:

```
Tb927.11.7160
Tb927.10.13720
Tb927.4.1610
```

By default, discoba_monomer_batch runs `colabfold_batch` using the following options:

```
# discoba_monomer_batch default options
--num-models 5
--num-recycle 3
--rank plddt
--recycle-early-stop-tolerance 0.5
--num-relax 1
--use-gpu-relax
```

If you want to use a custom AF2 configuration, see **"Using custom AF2.config file"** section.

NOTE: `discoba_monomer_batch` does not support automatic parallel GPU processing yet (`-g` flag). For a manual parallel processing implementation, see **"Using múltiple GPUs for parallel computing"** section.

## Some additional options

### Using custom AF2.config file
Sometimes you want to run ColabFold's AF2 algorithm with custom options. You can create an `AF2.config` file and call it with the `-c` flag pointing to the `AF2.config` file path. For example, for a batch of complex structure predictions:

```
# Run batch of complex structures predictions with custom AF2 options
discoba_multimer_batch -ma -c AF2.config database.fasta IDs_table.txt 2>&1 | tee report4.log
```

You can find a sample of the AF2.config file in `DiscobaMultimer/scripts/AF2.config` path. They look like this:

```
# You can add comment lines to this file
--num-models 5
--num-recycle 20
--rank iptm
# Stops when iptm value has been reached
--stop-at-score 75
--recycle-early-stop-tolerance 8.0
# Only relax rank_001 prediction
--num-relax 1
--use-gpu-relax
--save-all
```

For a full list of available options, run `colabfold_batch -h`.

NOTE: When you run the pipeline in a project folder that already contains the AF2 predictions for the IDs file, the predictions will not be performed. You need to create a new project directory and run again the pipeline pointing to the MSA files with the `-i` flag (`-i ../old_project/merged_MSA`, see **"Performing AF2 predictions on already computed MSAs"**).

### Setting size restrictions to avoid memory errors
Some combinations of IDs can result in huge complexes that can not be processed due to GPU memory limitations. To skip the computation of these combinations, you can add the `-s MIN_MAX` flag with the minimum (MIN) and maximum (MAX) size range. For example, if you want to compute only combinations of IDs that results in a summ of residues less than 2000, you may use the following call:

```
# Run batch of complex structures predictions with custom AF2 options
discoba_multimer_batch -ma -s 0_2000 database.fasta IDs_table.txt 2>&1 | tee report5.log
```

This means that only IDs combinations that ends up with a summ of residues between 0 and 2000 will be computed. The rest will be ignored and an `ignored.txt` file will be created in the project directory. This file will contain a description of why the prediction was not performed and pointing to the MSA file that was skipped to keep track on them. An example `ignored.txt` file will look as follows:

```
ignored:AF2	reason:combined_size(2542)>max_size(2000)	msa_file:../small_networks/merged_MSA/Tb927.11.6350__vs__Tb927.11.6350.a3m
```

This means that AF2 prediction was ignored, because the combined size was of 2542 and the MAX size value was set to 2000.

### Performing AF2 predictions on already computed MSAs from another projects
If you already computed all the MSAs for an `IDs_table.txt` file using `-m` only flag, you can point to them from another project folder by removing `-m` and using the `-i` flag. Here you have 3 examples:

```
# Run using merged_MSAs from another project folder
discoba_multimer_batch -a -i ../other_project_folder_path/merged_MSA database.fasta IDs_table.txt 2>&1 | tee report6.log

# Run using ColabFoldMSAs from another project folder
discoba_multimer_batch -a -i ../other_project_folder_path/colabfold_MSA database.fasta IDs_table.txt 2>&1 | tee report7.log

# Run using DiscobaMSA (without ColabFoldMSA merging) from another project folder
discoba_multimer_batch -a -i ../other_project_folder_path/discoba_paired_unpaired database.fasta IDs_table.txt 2>&1 | tee report8.log
```

This will compute AF2 models using these A3M files instead of generating A3M files again.

NOTE: If you want to use your custom MSAs created with other pipelines, I recommended to just use `colabfold_batch` with the A3M file as input.

### Using múltiple GPUs for parallel computing
On systems with múltiple GPUs, you can use the `-g GPUs_to_use` flag to allow parallel processing of the `IDs_table.txt` file. For example:

```
# Run batch of complex structures predictions using automatic parallel processing
discoba_multimer_batch -ma -g 8 database.fasta IDs_table.txt 2>&1 | tee report9.log
```

This will use 8 GPUs in parallel. Internally, DiscobaMultimer will split `IDs_table.txt` into 8 and excecute each splitted file into a different GPU, while masking the other GPUs using CUDA_VISIBLE_DEVICES. This means that each GPU will not share VRAM, nor GPU usage with other GPUs.

Alternatively, you can do a more trivial parallel processing excecution by manually splitting `IDs_table.txt` and calling DiscobaMultimer on separate CUDA devices manually:

```
# Run batch of complex structures predictions wusing manual parallel processing (GPU ID 0)
CUDA_VISIBLE_DEVICES=0 discoba_multimer_batch -ma database.fasta IDs_table_1.txt 2>&1 | tee report_GPU0.log

# Run batch of complex structures predictions wusing manual parallel processing (GPUs IDs 1 and 2, they will share VRAM) 
CUDA_VISIBLE_DEVICES=1,2 discoba_multimer_batch -ma -g 8 database.fasta IDs_table_2.txt 2>&1 | tee report_GPU1and2.log

# Run batch of complex structures predictions wusing manual parallel processing (GPU ID 3)
CUDA_VISIBLE_DEVICES=3 discoba_multimer_batch -ma -g 8 database.fasta IDs_table_3.txt 2>&1 | tee report_GPU3.log
```

NOTE: Do not execute multiple `-g` flag calls on the same system, because DiscobaMultimer is not aware of which GPUs are already processing predictions and may lead to errors.

### Generate MSA plots for generated DiscobaMSAs
To get a visual representation of how many sequences are retrived for each target use the `-p` flag.

```
# Compute DiscobaMSAs and get their MSA plots
discoba_multimer_batch -mp database.fasta IDs_table.txt 2>&1 | tee report10.log
```

This will create an MSA plot for each ID line on `IDs_table.txt`. The results will be stored in `msa_plots` directory.

NOTE: MSA plot generation is not so time efficient. So, do not use it for large datasets.

### Generate RoseTTAFold 2-tack predictions (beta)
If you want to retrive contact information using the RoseTTAFold 2-track model, you can do it by adding the `-r` tag:

```
# Compute DiscobaMSAs and use them to retrive contact info with RF 2-track model
discoba_multimer_batch -mr database.fasta IDs_table.txt 2>&1 | tee report10.log
```

Some plots with contact information, a metrics file, 2 NPZ files, and some other stuff will be created in `RoseTTAFold_2track_results` directory for each IDs line on `IDs_table.txt`.

NOTE:
For this, you need to install RF 2-track:

```bash
# Go to home directory
cd ~/

# Install Miniconda, read agreement, type yes two times and source bashrc
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.3.1-0-Linux-x86_64.sh
bash Miniconda3-py38_23.3.1-0-Linux-x86_64.sh
source .bashrc

# Clone the repo
git clone https://github.com/RosettaCommons/RoseTTAFold.git
cd RoseTTAFold

# create conda environment for RoseTTAFold (NVIDIA driver compatible with cuda11)
conda env create -f RoseTTAFold-linux.yml

# Download model neural-net weigths (includes 2-track model)
wget https://files.ipd.uw.edu/pub/RoseTTAFold/weights.tar.gz
tar xfz weights.tar.gz

# Install other dependencies
./install_dependencies.sh
```

If you already have miniconda/anaconda installed, skip the first part. RF 2-track usage is still on beta, but you can give it a try if you are interested.
    

## References
- Original Discoba database paper:
  
- DiscobaMultimer application on Histone Acetylation Related Complexes preprint:
