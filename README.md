# DiscobaMultimer
DiscobaMultimer is designed to perform highthroughput complex predictions with AlphaFold2-multimer on Discoba species, using LocalColabFold as AF2 implementation. TaxIDs for Discoba species where added to allow the generation of multiple sequence alignments using Discoba database (DiscobaMSA) for complexes, combining it or not with the corresponding ColabFoldMSA generated by ColabFold MSA server. There is also an implementation to run monomer structure predictions or just MSA retrivals.

For MSA generation only, no GPU is required. Take a look at this Colab notebook:
*To implement...

For complex structure predictions, GPU is required. Take a look at this notebook:
*To implement...

For monomeric structure predictions, GPU is required. Take a look at this notebook:
*To implement...

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
    If you get `colabfold_batch: command not found` error, install it following the installation instructions here: https://github.com/YoshitakaMo/localcolabfold.
  - Make sure you have Biopython package:
    ```bash
    pip install biopython    
    ```
  - Install emboss and hhsuite:
    ```bash
    sudo apt-get -y install emboss hhsuite
    ```
  - miniconda/anaconda???  <<<------------ (CHECK THIS ONE!!!) ----------------<<<

Once you are sure the dependencies are met, clone the repo and perform the installation:

```bash
# Clone the repo and cd to DiscobaMultimer
git clone https://github.com/elviorodriguez/DiscobaMultimer.git
cd DiscobaMultimer

# Install DiscobaMultimer
bash install/install_MMseqs2_and_DiscobaDB.sh install
```

This will take a few minutes, as MMseqs2 is installed and DiscobaDB is built. In general, the istallation of LocalColabFold adds mmseqs to the path, and MMseqs2 installation is skipped. After installation is completed, restart the shell (or `cd ~/` and then `source ~/.bashrc`).

Check the installation with:
```bash
# This will show you the help message for complex prediction
discoba_multimer_batch


# This will show you the help message for monomer prediction
discoba_monomer_batch
```

## How does it work?
DiscobaMultimer requieres 2 positional arguments, a file with fasta sequences `database.fasta` (in general will be the proteome) and a file with combinations of **IDs separated with TABS** `IDs_table.txt` to search in the database:

```bash
discoba_multimer_batch [OPTIONS] <database.fasta> <IDs_table.txt>
discoba_monomer_batch [OPTIONS] <database.fasta> <IDs_table.txt>
```

The options will manage what you want to retrieve (MSA, MSA plot, AF2 predictions, etc) and the configuration of each run.

The fasta sequences on `database.fasta` must contain headers **ONLY with the ID of each sequence**. For example:

```bash
>Tb927.11.7160
MEQQGDSEAKEHGELTRNGITDNAYDNSVGAFNRTILDRYAHWMTQLDQGACGLRIWKAE
ELWKRYERVIRVGRGSFGSVFIVYDTERKAYLTVKCMELLGKPGPALRSLSQPTLREVIL
LSQIDHPNVVRLIDYYLTSDGMLHMCMPIVSHDLVSLIRIWKMRGPRGESSLGRMPLPTV
KCVFRQLLRGLEYLHRRNIIHRDLKPSNVMLDDNGVVKIVDFGWARFVPRRWQGRLTGPP
CVVTYRPPEILLGGQCSFKYDSSIDIWSAGCILYEMLTGGKAFSNARNEQQALAAITDML
GSPSSRSEVYYGAAGGSRLRPSKRQPRNFEERCRMVNMSNESIDFLGEMLQLEPNSRKSA
SQLLGHSWFSTSPLPCEPEEVSLPGSNTYRLLERKRTR
>Tb927.10.13720
MAHVGQQRFGCYIGNIDRSVTLEVLRQVFSQCGTIVDCSLNGRDEDPYRYGFIDFATEDD
RARAMKYNGFTLAGRKIKVGISKGNVGRPEGYNNNPTPAPAASNTASSAGQHPQVPSQPQ
VPAVSVPASFLPGMVQQQQQQGATLLLQLLQQGAIDVNNLTAEQQQVLMASLLPQAPAAA
PGMHVMPPTAPMAYVPPPPQQPWGAPRGVYAGGPLGRPAPYTRPPANPQPPEETLKLREV
QRKQFLDVVRRDAEKYERKLAERNLKEGRTGSISGSEESSSDEEGEKGHRHSRRKIEGGE
DTEPEKSATLPMKTESDSVSCPMEVNCNNEGANISGEEAASNGNGGESNNNEDGDAIDCS
NVETEENNVDEEVENKC
>Tb927.4.1610
MIPPAKLNDFFNIVDDFLKKTFRDESFLFAAVKSRQYSESFPGEQLFFTSPPSATTDQPS
DESLQAGGGELRKTQNFMYFSPRLIFNRDGGYHGKVKLHSGVHVPQVCRFEQGIAVNSEG
LASGSVKVSDLLEGMEVKGRLAVNTIAPPSKDVWSVAMDYQRRDFYSTLNYQRNGLGSSD
LLVDCGTKFFNLLAGAGFERQKVSFLEQQDHTAQLDVLYAGVGFTGVNWSVGAKLVRAND
MWSAARIAFYQRVVPDTSVACAYNFDMEESRVHVSLGFSQGFRLRVPTILQQRACEQLDV
WTAILPFVGAFKAESGGLCAATIRGIFNGVVHWGLVAQKNVLVENSPIRFGLTLSVESG
```

For `IDs_table.txt`, each line will correspond to a prediction. For example:

```
Tb927.11.7160  Tb927.10.13720
Tb927.11.7160  Tb927.4.1610  Tb927.4.1610
Tb927.11.7160  Tb927.11.7160
```

This means that DiscobaMultimer will generate predictions for an heterodimer (first line), an heterotrimer (second line) and a homodimer (third line). Do not separate the the IDs with spaces or commas, only **TABS** will work.

NOTE: If you want to modify the headers of each sequence to other names, you can. Just do not add spaces or "-" to the new names. Also, make the corresponding change in the IDs file.


## Running batchs of complex structures predictions
Let's use `database.fasta` and `IDs_table.txt` from before as example. First, create a directory for the project to store databases and IDs table files and `cd` to it. As we want to run a batch of complexes, we need to use `discoba_multimer_batch` alias. The most basic call would be as follows:

```
# Run batch of complex structures predictions with defauls
discoba_multimer_batch -ma database.fasta IDs_table.txt 2>&1 | tee report.log
```

The `-m` tag tells DiscobaMultimer to generate DiscobaMSAs for each ID line, and `-a` tells that AF2 predictions must be performed on those DiscobaMSAs. It is always recommended to redirect both stdout and stderr to a log file to keep track the progress and check for any errors (_e.g._: using `2>&1 | tee report.log`). DiscobaMultimer will first generate all the MSAs and, once completed, it will start to produce the AF2-multimer predictions.

The resulting filesystem will be as follows:
```bash
.
├── AF2                                           <------- AlphaFold2-multimer predictions
│   ├── Tb927.11.7160__vs__Tb927.10.13720
│   ├── Tb927.11.7160__vs__Tb927.4.1610__vs__Tb927.4.1610
│   └── Tb927.11.7160__vs__Tb927.11.7160
├── colabfold_MSA                                 <------- ColabFoldMSAs produced by ColabFold MSA server
│   ├── Tb927.11.7160__vs__Tb927.10.13720.a3m
│   ├── Tb927.11.7160__vs__Tb927.4.1610__vs__Tb927.4.1610.a3m
│   └── Tb927.11.7160__vs__Tb927.11.7160.a3m
├── database.fasta                                <------- Original FASTA database (proteome)
├── discoba_mmseqs_alignments                     <------- DiscobaDB monomer MSAs (MMseqs2 output)
│   ├── Tb927.4.1610
│   ├── Tb927.10.13720
│   └── Tb927.11.7160
├── discoba_paired_unpaired                       <------- DiscobaMSA (A3M files)
│   ├── Tb927.11.7160__vs__Tb927.10.13720.a3m
│   ├── Tb927.11.7160__vs__Tb927.4.1610__vs__Tb927.4.1610.a3m
│   └── Tb927.11.7160__vs__Tb927.11.7160.a3m
├── IDs_table.txt                                 <------- Original IDs table
├── merged_MSA                                    <------- ColabFoldMSAs+DiscobaMSAs (A3M files)
│   ├── Tb927.11.7160__vs__Tb927.10.13720.a3m
│   ├── Tb927.11.7160__vs__Tb927.4.1610__vs__Tb927.4.1610.a3m
│   └── Tb927.11.7160__vs__Tb927.11.7160.a3m
└── report.log
```

By default, discoba_multimer_batch uses the following configuration to run `colabfold_batch`:

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
If you want to use a custom AF2 configuration, see **"Using custom AF2.options file"** section.

NOTE: All the filesystem is based on the `ID1 + "__vs__" + ID2 + "__vs__" + ...` notation. If you modify this filesystem nomenclature and then you run additional DiscobaMultimer pipelines with it, your pipeline will not work.

## Running batchs of DiscobaMSA only predictions (without AF2 predictions)
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

NOTE: This is useful to save resources. For example, when you run DiscobaMultimer on AWS or other cloud computing service, you can first run the MSA section using only the `-m` flag on a low price instance (_e.g._,without GPU). After the MSAs are generated, you can switch to a higher capacity instance (with multiple GPUs) and run the pipeline again, but this time with both flags: `-ma`. As the MSAs are already in the filesystem, they will not be computed, and it will jump directly to AF2 section. This will save you a lot of money. For parallel processing using multiple GPUs, see  **"Using múltiple GPUs for parallel computing"** section.

## Running batchs of monomer structures predictions


```
# discoba_monomer_batch default options
--num-models 5
--num-recycle 3
--rank plddt
--recycle-early-stop-tolerance 0.5
--num-relax 1
--use-gpu-relax
```




## Some additional options

### Using custom AF2.options file

### Setting size restrictions to avoid memory errors

### Performing AF2 predictions on already computed MSAs

### Using múltiple GPUs for parallel computing

### Generate MSA plots for generated DiscobaMSAs

### Generate RoseTTAFold 2-tack predictions (beta)

## References
- Original Discoba database:
  
- DiscobaMultimer application:
