# aKmerBroom: Ancient oral DNA decontamination using Bloom filters on k-mer sets

This tool identifies ancient reads, given a file of known ancient kmers. It does so in the following steps: 
1. Build an `ancient_kmers.bloom` filter from an ancient kmers text file (if such a Bloom filter does not yet exist).
2. For a set of input reads:
    1. Save those reads which have 2 consecutive kmer matches against `ancient_kmers.bloom`
    2. Kmerize the saved reads to generate a new set of ancient kmers, called "anchor kmers"
3. For the same set of input reads, identify matches against anchor kmers and classify each read with >50% matches as an ancient read.


### Usage
    # Use the ancient kmers bloom filter provided
    python akmerbroom.py --ancient_bloom

    or    

    # Use an ancient kmers text file 
    python akmerbroom.py --ancient_kmers_set

    


### Input

The `data/` folder should contain the following input files:

```
ancient_kmers.bloom : a bloom filter with ancient kmers
unknown_reads.fastq : a file with reads which we want to classify as ancient or not
[optional] ancient_kmers : a text file where each row is a known ancient kmer
```    

### Output 

The `output/` folder should contain the following output files:
```
annotated_reads.fastq                     # intermediate output
annotated_reads_with_anchor_kmers.fastq   # final output
```
The final output file has the following 4 fields in each record header: 
```
SeqId, ReadLen, isConsecutiveMatchFound, AnchorProportion
```   
By default, reads with `AnchorProportion` >= 0.5 (ie. 50%) are chosen as ancient reads. 


### Dependencies
```
pip install biopython
pip install cython
pip install pybloomfiltermmap3
```

    
