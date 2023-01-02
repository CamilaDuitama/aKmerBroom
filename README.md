# ancient_dna

This tool identifies ancient reads, given a file of known ancient kmers. It does so in the following steps: 
1. Build an `ancient_kmers.bloom` filter from an ancient_kmers text file (if such a Bloom Filter does not yet exist).
2. For a set of input reads:
    1. Save those reads which have 2 consecutive kmer matches against `ancient_kmers.bloom`
    2. Kmerize the saved reads to generate a new set of ancient kmers, called "seen kmers"
3. For the same set of input reads, for each read, identify matches against "seen kmers" and classify each read with >50% matches as an ancient read.


### Usage

    # Use the ancient kmers text file and run the method
    python ancient_dna.py --ancient_kmers_set

    # Or alternatively, use the ancient kmers bloom filter and run the method
    python ancient_dna.py --ancient_bloom


### Input

The `data/` folder contains the following input files:

```
ancient_kmers       : a text file where each row is a known ancient kmer
ancient_kmers.bloom : optionally, a bloom filter version of the ancient_kmers text file
unknown_reads.fastq : a file with reads which we want to classify as ancient or not
```    

### Output 

The `output/` folder contains the following outputs files:
```
annotated_reads.fastq                    # intermediate output
annotated_reads_with_seen_kmers.fastq    # final output
```
The final output file has the following 5 fields in each record header: 
```
SeqId, ReadLen, InitialAncientKmerProportion, isConsecutiveMatchFound, SeenProportion
```   
Reads with `SeenProportion` > 0.5 (ie. 50%) can be considered as ancient reads. 


### Dependencies
```
pip install biopython
pip install cython
pip install pybloomfiltermmap3
```

    
