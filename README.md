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
  - Make sure curl, git and wget are installed.
    ```bash
    sudo apt-get -y update
    sudo apt-get -y install curl git wget
    ```
  - a
  - a
  - a
  - a

```bash
# Clone the repo and cd to DiscobaMultimer
git clone https://github.com/elviorodriguez/DiscobaMultimer.git
cd DiscobaMultimer
```

```bash
# Make sure the dependencies are met
cat DiscobaMultimer
```
