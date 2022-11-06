# ancient_dna

This repo hosts an early version of a tool that classifies reads as ancient or not, given a text file of known ancient kmers. 



### Usage

    python ancient_dna.py --ancient_kmers_set


### Input

The `data/` folder contains the following input files:

```
ancient_kmers : a text file where each row is a known anciant kmers
unknown_reads.fastq : a file with reads which we want to classify as ancient or not
```    

### Output 

The `output/` folder contains the output file with reads classified as ancient. 



### Dependencies
```
pip install biopython
pip install cython
pip install pybloomfiltermmap3
```

    
