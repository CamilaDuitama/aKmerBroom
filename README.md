# 🧹🦷rKmerBroom: Modified version of aKmerBroom, see original repo (temporary name to avoid shadowing)
## Read paper [here](https://www.cell.com/iscience/pdf/S2589-0042(23)02134-X.pdf)

`aKmerBroom` is a tool to decontaminate ancient DNA samples. This version differs from the original in the following aspects:
1. The input is a list of files (space separated). THere is no enforcing of a folder input name.
2. Output folder is created if needed, otherwise there is just a warning it already exists and asks user what to do, log the decision.
3. Each input file is split into two output files, both using the same name prefix. One is the decontaminated sequences, the other contains the contaminating sequences.
4. Each input file is processed in parallel.
5. Once the anchor kmers of each input file are computed, all the sets are combined into a union of the sets.

![pipeline_svg.png](https://raw.githubusercontent.com/CamilaDuitama/aKmerBroom/main/pipeline_svg.png)**aKmerBroom pipeline:** First, an offline step is performed: a collection of samples representative from diverse sources is used to create a trusted set of oral kmers. The trusted collection indexes kmers that appear exclusively in modern and ancient oral samples, but not other samples from contaminant sources (see panel on the left called Collection of datasets). Then this set of oral kmers is used to decontaminate an input set of reads. The algorithm proceeds by looking up each read kmer inside the Bloom Filter of trusted oral kmers, and marking positions of matches. Reads having at least two consecutive matches to the Bloom Filter get passed to the construction of a set containing all kmers from such reads. Finally, the same input reads are scanned again using the aforementioned set, and reads having a proportion of kmer matches over a certain threshold are reported to be of ancient oral origin.

+ [Usage](#Usage)
+ [Input](#Input)
+ [Output](#Output)
+ [Dependencies](#Dependencies)
+ [Testing](#Testing)

## Usage
```bash
usage: akmerbroom.py [options]

This program uses a reference (either a bloom filter or a kmer text file) to recognise targeted DNA (ancient or modern) reads to separate ancient DNA from modern DNA. aDNA will be stored in a "decontaminated" file and the modern DNA in the "contamination" for each sample.

options:
  -h, --help            show this help message and exit
  --bloom BLOOM         Used if a BloomFilter is provided (defaults to False)
  --bloom_capacity BLOOM_CAPACITY
                        If a BloomFilter is not provided, this sets the capacity of the bloom filter. This should be greater than the number of distinct kmers in the input file. Default to 2 billion.
  --kmers_set KMERS_SET
                        Used if a kmers set is provided (defaults to False).
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Set kmer size (defaults to 31)
  --n_consec_matches N_CONSEC_MATCHES
                        Set number of consec matches to classify read as anchor read, (defaults to 2).
  --anchor_proportion_cutoff ANCHOR_PROPORTION_CUTOFF
                        Set anchor kmer proportion, above which a read is classified as modern/ancient (defaults to 0.5)
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Path to input file(s), space-separated.
  -o OUTPUT, --output OUTPUT
                        Path to output folder, where you want aKmerBroom to write the results.
  -t THREADS, --threads THREADS
                        WARNING: right now, not used. Sorry, async is a pain. Number of threads to use, default to 1.
  -s, --single          Flag to decontaminate samples independently instead of pooling k-mers from multi-samples for decontamination.
  -m, --modern          Flag to indicate that reference is modern DNA (defaults to False).

EX: python3 akmerbroom.py -i $(find ./tests/*.fastq) -o ./output -t 2 --kmers_set kmers.txt

```

## Input

Other than fastq files (that you can provide using `find` command), either provide a kmer set from modern DNA or a bloom filter built with kmers from modern DNA.

## Output 

The output folder contains, for each file `{input}.fastq`, files `{input}_annotated.fastq`, `{input}_contamination.fastq` and `{input}_decontaminated.fastq.`.
The first contains annotated sequences, the second the rejected sequences (not ancient), the third contains the sequences accepted as ancient.

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
The `tests/` folder contains a test dataset consisting of 3 fastq files a kmer set. Each fastq file contains only 1 or 2 sequences.
To run a test, use the following steps:

```bash
python3 ../rKmerBroom/akmerbroom.py --input $(find tests/*.fastq) --output ./output/ --kmers_set tests/kmer_set.txt 
```
To verify the behavior of the `--modern` option:
```bash
python3 ../rKmerBroom/akmerbroom.py --input $(find tests/*.fastq) --output ./output_modern/ --kmers_set tests/kmer_set.txt --modern
```
Compare the output files : contamination and decontaminated content should be reversed.
To verify the behavior of the `--single` option:
```bash
python3 ../rKmerBroom/akmerbroom.py --input $(find tests/*.fastq) --output ./output_modern/ --kmers_set tests/kmer_set.txt --single
```
This time, the sequence number 5 (`seq5`) should not be recognised. This is because the reference is augmented separately for each sample, thus the kmers containing G, added to the reference pool because of `seq2`, will not be present in the reference at the time this sample is decontaminated.
Whether the sequence is in `contamination` or `decontaminated` depends on the use of the flag `--modern`. 
