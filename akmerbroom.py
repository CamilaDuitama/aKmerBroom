import argparse
from scripts import classify_reads, kmers


def main():
    """
    The following input files are required in ./data folder
    1. unknown.fastq : input file that needs to be classified
    2. ancient.bloom or ancient_kmers
    """
    parser = argparse.ArgumentParser(description='This program finds ancient reads')
    parser.add_argument('--ancient_bloom', help='Use if ancient BloomFilter provided (defaults to False)',
                        action='store_true', required=False)
    parser.add_argument("--ancient_bloom_capacity", type=int, help="If ancient BloomFilter is not provided, \
        This sets the capacity of the bloom filter. This should be greater than the number of distinct kmers \
        in ancient kmers input file")
    parser.add_argument('--ancient_kmers_set', help='Use if ancient kmers set provided (defaults to False)',
                        action='store_true', required=False)
    parser.add_argument('--kmer_size', help='Set kmer size (defaults to 31)', required=False)
    parser.add_argument('--anchor_proportion_cutoff', help="Set anchor kmer proportion, \
        above which a read is classified as ancient (defaults to 0.5)",
                        required=False)
    args = vars(parser.parse_args())

    bloom_filt = False
    bf_capacity = 0
    ancient_kmers = False

    # set kmer_size from argument or using default here
    if not args['kmer_size']:
        k_size = 31
    else:
        try:
            k_size = int(args['kmer_size'])
        except ValueError:
            print("Error : kmer_size provided is not an integer")
            kmers.exit_gracefully()

    # set bloom filter
    if args['ancient_bloom']:
        bloom_filt = True
        bf_capacity = None
    elif args['ancient_bloom_capacity']:
        try:
            bf_capacity = int(args['ancient_bloom_capacity'])
        except ValueError:
            print("Error : ancient_bloom_capacity provided is not an integer")
            kmers.exit_gracefully()
    else:
        bf_capacity = 1000 * 1000 * 1000 * 2

    # ancient kmers set
    if args['ancient_kmers_set']:
        ancient_kmers = True

    # set anchor proportion cutoff
    if not args['anchor_proportion_cutoff']:
        anchor_proportion_cutoff = 0.5
    else:
        anchor_proportion_cutoff = float(args['anchor_proportion_cutoff'])

    print("Started...")
    # declare defaults
    print("Using default of k=31 and input folder='data'")
    print("Shortlisting ancient reads")
    anchor_kmer_set = classify_reads.classify_reads(bloom_filt, bf_capacity, ancient_kmers, k_size)
    classify_reads.classify_reads_using_anchor_kmers(anchor_kmer_set, k_size, anchor_proportion_cutoff)
    print("Completed successfully")


if __name__ == "__main__":
    main()

