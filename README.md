# GeneToCN
Gene copy number prediction from _k_-mer frequencies

## Setup:
Pull Github repository (or download only GeneToKmer.py and KmerToCN.py).
Also download glistmaker, glistquery, glistcompare, gmer_counter tools from GenomeTester4 toolkit (https://github.com/bioinfo-ut/GenomeTester4 - see the manual for GenomeTester4 and FastGT)

GeneToCN is implemented on Python 3 and requires following packages: biopython, pandas, matplotlib.

We suggest using anaconda and create a new virtual enivironment 
``` 
conda create --name <env_name> biopython pandas matplotlib
```

### Creating a _k_-mer list for reference genome (needed for GeneToKmer):
Outputs a binary file <output_prefix>_<k-mer length>.list, the suggested k-mer length for human genome is 25
```
glistmaker <reference_sequences> -w <k-mer length> -o <output_prefix>
```

### Before running the programs:
Activate the conda environment, if not already activated
```
conda activate <env_name>
```

## _K_-mer database creation
### _K_-mer databases for genes:
```
python GeneToKmer.py <region_file> <chr_reference_sequence>.fa <reference_kmer_list> -o <location_for_output_files> -i -gt <location_for_GenomeTester_source_code>
```
Needs to be run separately for each chromosome, requires the fasta sequence for only one chromosome. The region (gene) names have to be unique.

Example of the region file (fields separated by a space):
```
AMY1 G 1:103655290-103664554 1:103687415-103696680 1:103750406-103758690
AMY2A G 1:103616811-103625780
AMY2B G 1:103553815-103579534
```

### Databases of flanking referenence _k_-mers:

Option 1 (faster) - use premade _k_-mer databases for human genome (GRCh38):
```
python extract_flanking_kmers.py <region_file>
```
To change the default option for the number of flanking _k_-mers, use -_n_ <nr_of_flanking_kmers> (default: _n_=2000).

The script takes a subset of at least _n_/2 reference _k_-mers from both sides of each gene (a total of _n_ kmers) using the _k_-mers in "Kmer_db/Ref_kmers" for human genome. One _k_-mer database is created for each gene (may contain the same set of _k_-mers for closely located genes).

Option 2:

Add flanking regions to the region file to create _k_-mer databases when running GeneToKmer.py
```
Flanking_AMY F 1:103500000-103550000 1:103760000-103800000
```

## Copy number estimation

### Additional files needed
Cat together all the _k_-mer database files to be able to estimate all the copy numbers with the same run (the order is not important).

Create a file containing the original _k_-mer db locations for each gene. Genes using one flanking region should be separated with an empty row and flanking db always before the gene (or genes) for which it should be used. Fields separated by a space.

An example:
```
Flanking_AMY Flanking_AMY_NIPT_103M_kmers_cleaned.db
AMY1 AMY1_kmers.db
AMY2A AMY2A_kmers.db
AMY2B AMY/AMY2B_kmers.db

Flanking_NPY4R NPY4R/Flanking_NPY4R_NIPT_46M47M_kmers_cleaned.db
NPY4R NPY4R/NPY4R_kmers.db
```

### Running KmerToCN
The program will create new directories (counts/ counter_errors/ plots/) into the specified output directory containing the _k_-mer counts, possible errors from the gmer_counter and plots for visual inspection of the normalized _k_-mer frequencies in the gene regions. 

The copy number estimations are written in the outpur direction into file <output_prefix>_cn.txt

Use -i to print out some information while running
```
python KmerToCN.py -db <kmer_db> -kp <original_db_location_file> -s <fastq files> -d <output_directory> -o <output_prefix> -gm  <location_for_gmer_counter> -r <region_file>
```


## _K_-mer databases
Some _K_-mer databases and region files are included in the folder **Kmer_db**.

The _k_-mers in **Ref_kmers** folder can be used for choosing reference k-mers from flanking regions.

**Gene_kmers** folder contains _k_-mers used for the validation of the method.

All _k_-mer databases are text files formatted in a similar way with one _k_-mer per row: 
```
chr:position nr_of_kmers kmers
```

An example:
```
1:103656318     1       CTTCAAAGCAAAATGAAGCTCTTTT
1:103656319     1       TTCAAAGCAAAATGAAGCTCTTTTG
1:103656320     1       TCAAAGCAAAATGAAGCTCTTTTGG
1:103656321     1       CAAAGCAAAATGAAGCTCTTTTGGT
1:103656322     1       AAAGCAAAATGAAGCTCTTTTGGTT
```

## Cite
GeneToCN method as well as the results of the validation are described in the paper.

Pajuste, FD., Remm, M. GeneToCN: an alignment-free method for gene copy number estimation directly from next-generation sequencing reads. _Sci Rep_ **13**, 17765 (2023). https://doi.org/10.1038/s41598-023-44636-z

## Contact information
Fanny-Dhelia Pajuste, University of Tartu

fanny-dhelia.pajuste@ut.ee
