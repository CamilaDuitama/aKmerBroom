# ancient_dna

This repo hosts an early version of a tool that classifies reads as ancient or not, given a text file of known ancient kmers. It does so in the following steps: 
1. Build a Bloom Filter from a text file of ancient kmers 
2. For a set of input reads, classify as ancient those reads which have greater than 5% of their kmers present in the Bloom Filter generated in step 1



### Usage

    python ancient_dna.py --ancient_kmers_set


### Input

The `data/` folder contains the following input files:

```
ancient_kmers : a text file where each row is a known ancient kmer
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

    
