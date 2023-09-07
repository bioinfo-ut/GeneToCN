# GeneToCN
Gene copy number prediction from k-mer frequencies

## Setup:
Pull Github repository (or download only GeneToKmer.py and KmerToCN.py)
Also download glistmaker, glistquery, glistcompare, gmer_counter tools from GenomeTester4 toolkit (https://github.com/bioinfo-ut/GenomeTester4 - see the manual for GenomeTester4 and FastGT)

GeneToCN is implemented on Python 3 and requires following packages: biopython, pandas, matplotlib
We suggest using anaconda and create a new virtual enivironment 
```
module load <anaconda_location>  
conda create --name <env_name> biopython pandas matplotlib
```

### Creating a kmer list for reference genome (needed for GeneToKmer):
Outputs a binary file <output_prefix>_<k-mer length>.list
Suggested k-mer length is 25
```
glistmaker <reference_sequences> -w <k-mer length> -o <output_prefix>
```

## Run:
If conda is not already loaded, activate it
```
#module load <anaconda_location>
conda activate <env_name>
```

### Create k-mer lists 
```
python GeneToKmer.py <region_file> <chr_reference_sequence>.fa <reference_kmer_list> -o <location_for_output_files> -i -gt <location_for_GenomeTester_source_code>![image](https://github.com/bioinfo-ut/GeneToCN/assets/3532045/baa7b4ab-bee6-4390-81ab-db585def09de)
```

## K-mer databases
K-mer databases are included in the folder kmer db
The k-mers in ref kmers folder can be used for choosing reference k-mers from flanking regions. 
All k-mer databases are text files formatted in a similar way with on k-mer per row: 
```
chr:position nr_of_kmers kmers
```
