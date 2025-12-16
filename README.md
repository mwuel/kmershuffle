## Kmershuffle
A small script counting k-mers and z-scoring randomly shuffled k-mers of same length on the same genes given a bed file of short binding sites of RNA-binding proteins.

The script takes as input a GTF(GFF)-file containing all annotated features of the analyzed genotype, a BED-file containing all binding sites of the RNA-binding proteins, and the full sequence of the annotated genotype in FASTA-format. All binding sites are assumed to be of the same length. The script also filters for valid genes via a regex that might need to be tuned for a given dataset. The tool counts k-mers in all binding sites and performs z-scoring on a given number of random shufflings of binding sites. For every shuffling, the same amount of binding sites are randomly chosen on every corresponding gene and k-mers are counted for the new binding sites. The shufflings are then used to compute the mean k-mer count, standard deviation and finally the z-score for every k-mer. The output is given as a csv-file.

## Software and Versions
- Python (3.11)
- numpy (2.3.5)
- pandas (2.3.3)
- pybedtools (0.12.0)
- pysam (0.23.3)


## Usage
The script can be run using the following command:
```
python --anno kmershuffle.py --anno gtf.gtf --binds binding_sites.bed --seq sequence.fa k=3 shufflings=100
```

Parameters *k* and *shufflings* are optional but can be changed as desired. K-mer length should, however, be chosen according to site lengths (e.g k <= site length). Binding sites can also be given in GFF-format. 

To test the script, a sample dataset has been uploaded to the repository. It takes the data from the following publication: https://doi.org/10.1111/tpj.16601
