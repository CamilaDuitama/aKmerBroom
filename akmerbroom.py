import argparse
import os
import logging
from scripts import classify_reads, kmers


#TODO: faire un argument pour donner les inputs comme liste de fichiers ou un fof.txt
#TODO: parallelliser l'analyse des fichiers
def main():
    """
    The following input files are required in ./data folder
    1. unknown.fastq : input file that needs to be classified
    2. present.bloom or present_kmers
    """
    parser = argparse.ArgumentParser(description='This program finds removes present-day reads to only keep the ancient DNA.')
    parser.add_argument('--present_bloom', help='Use if present BloomFilter provided (defaults to False)',
                        action='store_true', required=False)
    parser.add_argument("--present_bloom_capacity", type=int, help="If present BloomFilter is not provided, \
        This sets the capacity of the bloom filter. This should be greater than the number of distinct kmers \
        in present kmers input file")
    parser.add_argument('--present_kmers_set', help='Use if present kmers set provided (defaults to False)',
                        action='store_true', required=False)
    parser.add_argument('--kmer_size', type=int, help='Set kmer size (defaults to 31)', required=False)
    parser.add_argument('--n_consec_matches', type=int, help="Set number of consec matches to classify read as anchor read, \
                                                   (defaults to 2)", required=False)
    parser.add_argument('--anchor_proportion_cutoff', help="Set anchor kmer proportion, \
        above which a read is classified as present (defaults to 0.5)",
                        required=False)
    parser.add_argument('--input', help = "Path to input folder or input file", required=False, default = "./data/")
    parser.add_argument("--output",
                        help="Path to output folder, where you want aKmerBroom to write the results.", required=True,
                        default="output")

    args = vars(parser.parse_args())

    bloom_filt = False
    bf_capacity = 0
    present_kmers = False

    # Create and configure logger
    logging.basicConfig(filename="aKmerBroom.log",
    format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG, filemode='w')

    # Creating an object
    logger = logging.getLogger()

    #check output folder does not exist else create folder
    if os.path.isdir(args["output"]):
        logger.warning("Output directory already exists")
    else:
        output=args["output"]
        os.mkdir(output)

    #check if input is folder or file else exit
    if os.path.isdir(args["input"]):
        input = [f for f in listdir(args["input"]) if isfile(join(args["input"],f))]
    elif os.path.isfile(args["input"]):
        input = args["input"]
    else :
        logger.error("Error: input does not exist.")
        kmers.exit_graccefully()

    # set kmer_size from argument or using default here
    if not args['kmer_size']:
        k_size = 31
    else:
        try:
            k_size = int(args['kmer_size'])
        except ValueError:
            logger.error("Error : kmer_size provided is not an integer")
            kmers.exit_gracefully()

    # set bloom filter
    if args['present_bloom']:
        bloom_filt = True
        bf_capacity = None
    elif args['present_bloom_capacity']:
        try:
            bf_capacity = int(args['present_bloom_capacity'])
        except ValueError:
            logger.error("Error : present_bloom_capacity provided is not an integer")
            kmers.exit_gracefully()
    else:
        bf_capacity = 1000 * 1000 * 1000 * 2

    # present kmers set
    if args['present_kmers_set']:
        present_kmers = True

    # set n_consec_matches cutoff
    if not args['n_consec_matches']:
        n_consec_matches = 2
    else:
        n_consec_matches = int(args['n_consec_matches'])

    # set anchor proportion cutoff
    if not args['anchor_proportion_cutoff']:
        anchor_proportion_cutoff = 0.5
    else:
        anchor_proportion_cutoff = float(args['anchor_proportion_cutoff'])

    logger.info("Started...")
    # declare defaults
    logger.info("Using default of k=31 and input folder='data'")
    logger.info("Shortlisting present reads")
    anchor_kmer_set = classify_reads.classify_reads(bloom_filt,
                                                    bf_capacity,
                                                    present_kmers,
                                                    k_size,
                                                    n_consec_matches,
                                                    output)
    classify_reads.classify_reads_using_anchor_kmers(anchor_kmer_set,
                                                     k_size,
                                                     anchor_proportion_cutoff,
                                                     output)
    logger.info("Completed successfully")


if __name__ == "__main__":
    main()

