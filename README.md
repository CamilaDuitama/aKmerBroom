# 🧹🦷rKmerBroom: Modified version of aKmerBroom, see original repo (temporary name to avoid shadowing)
## Read paper [here](https://www.cell.com/iscience/pdf/S2589-0042(23)02134-X.pdf)

`rKmerBroom` is a tool to decontaminate ancient DNA samples. It differes from aKmerBroom in the following aspects:
1. The input is a list of files (space separated). THere is no enforcing of a folder input name.
2. Output folder is created if needed, otherwise there is just a warning it already exists.
3. Each input file is split into two output files, both using the same name prefix. One is the decontaminated sequences, the other contains the contaminating sequences.
4. Each input file is processed in parallel.

![pipeline_svg.png](https://raw.githubusercontent.com/CamilaDuitama/aKmerBroom/main/pipeline_svg.png)**aKmerBroom pipeline:** First, an offline step is performed: a collection of samples representative from diverse sources is used to create a trusted set of oral kmers. The trusted collection indexes kmers that appear exclusively in modern and ancient oral samples, but not other samples from contaminant sources (see panel on the left called Collection of datasets). Then this set of oral kmers is used to decontaminate an input set of reads. The algorithm proceeds by looking up each read kmer inside the Bloom Filter of trusted oral kmers, and marking positions of matches. Reads having at least two consecutive matches to the Bloom Filter get passed to the construction of a set containing all kmers from such reads. Finally, the same input reads are scanned again using the aforementioned set, and reads having a proportion of kmer matches over a certain threshold are reported to be of ancient oral origin.

+ [Usage](#Usage)
+ [Input](#Input)
+ [Output](#Output)
+ [Dependencies](#Dependencies)
+ [Testing](#Testing)

## Usage
    # Use the ancient kmers bloom filter provided
    python akmerbroom.py --ancient_bloom

    or    

    # Use an ancient kmers text file 
    python akmerbroom.py --ancient_kmers_set



## Input

Other than fastq files (that you can provide using `find` command), either provide a kmer set from present DNA or a bloom filter built with kmers from present DNA.

## Output 

The output folder contains, for each file `{input}.fastq`, files `{input}_contamination.fastq` and `{input}_decontaminated.fastq.`.
The former contains sequences rejected (not ancient), the later the sequences accepected as ancient.

## Dependencies
```
pip install biopython
pip install cython
pip install pybloomfiltermmap3
pip install multiprocess
```

Or use the provided conda environment:
```bash
mamba create -n broom akmerbroom.yml
```

## Testing
The `tests/` folder contains a test dataset consisting of 3 fastq files a kmer set. Each fastq file contains only 2 sequences.
To run a test, use the following steps:

```bash
python3 ../rKmerBroom/akmerbroom.py --input $(find tests/*.fastq) --output output/ --present_kmers_set tests/kmer_set.txt
```
