# DiscobaMultimer
DiscobaMultimer is designed to perform highthroughput complex predictions with AlphaFold2-multimer on Discoba species, using LocalColabFold as AF2 implementation. TaxIDs for Discoba species where added to allow the generation of multiple sequence alignments using Discoba database (DiscobaMSA) for complexes, combining it or not with the corresponding ColabFoldMSA generated by ColabFold MSA server. There is also an implementation to run monomer structure predictions or just MSA retrivals.

For MSA generation only, no GPU is required. Take a look at this Colab notebook:
*To implement...

For complex structure predictions, GPU is required. Take a look at this notebook:
*To implement...

For monomeric structure predictions, GPU is required. Take a look at this notebook:
*To implement...

## Local installation
Before everything, make sure the dependencies are met:
  - Make sure curl, git and wget are installed:
    ```bash
    sudo apt-get -y update
    sudo apt-get -y install curl git wget
    ```
  - Make sure you have LocalColabFold installed:
    ```bash
    # This must give no errors
    colabfold_batch -h
    ```
    If you get `colabfold_batch: command not found` error, install it following the installation instructions here: https://github.com/YoshitakaMo/localcolabfold.
  - Make sure you have Biopython package
    ```bash
    pip install biopython    
    ```
  - Install emboss and hhsuite
    ```bash
    sudo apt-get -y install emboss hhsuite
    ```
  - miniconda/anaconda??? <------------ (CHECK THIS ONE!!!) ----------------

Once you are sure the dependencies are met, clone the repo and perform the installation:

```bash
# Clone the repo and cd to DiscobaMultimer
git clone https://github.com/elviorodriguez/DiscobaMultimer.git
cd DiscobaMultimer

# Install DiscobaMultimer
bash install/install_MMseqs2_and_DiscobaDB.sh
```

This will take a few minutes, as MMseqs2 is installed and DiscobaDB is built. In general, the istallation of LocalColabFold adds mmseqs to the path, and MMseqs2 installation is skipped. After installation is completed, restart the shell or `source ~/.bashrc`.

Check the installation with:
```bash
# This will show you the help message for complex prediction
discoba_multimer_batch


# This will show you the help message for monomer prediction
discoba_monomer_batch
```
