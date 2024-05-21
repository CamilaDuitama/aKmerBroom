import argparse
import os
import logging
from scripts import classify_reads, kmers
from multiprocessing import Pool, Value, Process
from multiprocessing.managers import BaseManager
from pybloomfilter import BloomFilter
from ctypes import c_wchar_p



def main():
    """
    The following input files are required in ./data folder
    1. unknown.fastq : input file that needs to be classified
    2. present.bloom or present_kmers
    """
    parser = argparse.ArgumentParser(
        description='This program finds removes present-day reads to only keep the ancient DNA.',
        epilog="Ex: python3 akmerbroom.py -i test1.fastq test2.fastq -t 2 --present_kmers_set kmers.txt")
    parser.add_argument('--present_bloom', help='Use if present BloomFilter provided (defaults to False)',
                        type=str, action='store', required=False, default="")
    parser.add_argument("--present_bloom_capacity", type=int, help="If present BloomFilter is not provided, \
        This sets the capacity of the bloom filter. This should be greater than the number of distinct kmers \
        in present kmers input file")
    parser.add_argument('--present_kmers_set', help='Use if present kmers set provided (defaults to False)',
                        type=str, action='store', required=False, default="")
    parser.add_argument('--kmer_size', type=int, help='Set kmer size (defaults to 31)', required=False)
    parser.add_argument('--n_consec_matches', type=int, help="Set number of consec matches to classify read as anchor read, \
                                                   (defaults to 2)", required=False)
    parser.add_argument('--anchor_proportion_cutoff', help="Set anchor kmer proportion, \
        above which a read is classified as present (defaults to 0.5)",
                        required=False)
    parser.add_argument('--input', help="Path to input file(s), space-separated", required=True, nargs='+')
    parser.add_argument("--output",
                        help="Path to output folder, where you want aKmerBroom to write the results.", required=True,
                        default="output")
    parser.add_argument("-t", "--threads", help="Number of threads to use", default=1, type=int)

    args = vars(parser.parse_args())

    #bloom_filt = False
    bf_capacity = 0
    #present_kmers = False

    # Create and configure logger
    logging.basicConfig(filename="aKmerBroom.log",
                        format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG, filemode='w')

    # Creating an object
    logger = logging.getLogger()

    #check output folder does not exist else create folder
    if os.path.isdir(args["output"]):
        logger.warning("Output directory already exists")
    else:
        output = args["output"]
        os.mkdir(output)

    #check if input files exist
    for i in args["input"]:
        if not os.path.exists(i):
            logger.error(f"Error: input {i} does not exist.")
            kmers.exit_gracefully()

    if args["threads"] < 1:
        logger.error("The number of threads is not valid")
        kmers.exit_gracefully()

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
    if not os.path.isfile(args['present_bloom']) and not os.path.isfile(args['present_kmers_set']):
        logger.error("Error: please provide a kmer set or a bloom filter.")
        kmers.exit_gracefully()
        #bloom_filt = True

        bf_capacity = None
    elif args['present_bloom_capacity']:
        try:
            bf_capacity = int(args['present_bloom_capacity'])
        except ValueError:
            logger.error("Error : present_bloom_capacity provided is not an integer")
            kmers.exit_gracefully()
    else:
        bf_capacity = 1000 * 1000 * 1000 * 2

    #TODO:verifier que je peux virer ce bout la
    # present kmers set

    if os.path.isfile(args['present_kmers_set']):
        present_kmers = True

    # set n_consec_matches cutoffTrue
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
    logger.info("Shortlisting present reads")
    # builds the bloom filter once for all
    present_kmers_bf = classify_reads.getbloomFilter(args["present_bloom"], bf_capacity, args["present_kmers_set"],
                                                                  k_size, args["output"])

    # create Values for shared memory for the parameters
    shared_k_size = Value('i', k_size)
    shared_n_consec_matches = Value('i', n_consec_matches)
    shared_output = Value(c_wchar_p, args["output"])
    shared_anchor_proportion_cutoff = Value('f', anchor_proportion_cutoff)
    # creates a custom manager for a set in shared memory
    class CustomManager(BaseManager):
        # nothing
        pass

    CustomManager.register('set', set)

    manager = CustomManager()
    manager.start()
    shared_kmer_set = manager.set()

    # paralleled anchor sets construction
    #with Pool(processes=args["threads"]) as pool:
        #for i in args["input"]:
    #        print(i)
            #classify_reads.classify_reads(i, present_kmers_bf, shared_k_size.value, shared_n_consec_matches.value, shared_output.value)
        #    pool.apply_async(classify_reads.classify_reads,
        #                     args=(i, present_kmers_bf, shared_k_size.value, shared_n_consec_matches.value, shared_output.value),
        #                     callback=collect_sets)
    #    pool.close()
    #    pool.join()

    processes = [Process(target=classify_reads.classify_reads, args=(i, present_kmers_bf, shared_k_size.value, shared_n_consec_matches.value, shared_output.value, shared_kmer_set)) for i in args["input"]]
    for process in processes:
        process.start()
    for process in processes:
        process.join()
    for process in processes:
        process.terminate()

    #with Pool(processes=args["threads"]) as pool:
    #    for i in args["input"]:
            #classify_reads.classify_reads_using_anchor_kmers(i, present_kmers_bf, k_size, anchor_proportion_cutoff, args["output"])
    #        pool.apply_async(classify_reads.classify_reads_using_anchor_kmers, args=(i, shared_kmer_set,
    #                                                                                 shared_k_size.value,
    #                                                                                 shared_anchor_proportion_cutoff,
    #                                                                                 shared_output.value,))
    #    pool.close()
    #    pool.join()

    processes = [Process(target=classify_reads.classify_reads_using_anchor_kmers, args=(i, shared_kmer_set, shared_k_size.value, shared_anchor_proportion_cutoff.value, shared_output.value)) for i in args["input"]]
    for process in processes:
        process.start()
    for process in processes:
        process.join()
    for process in processes:
        process.terminate()

    logger.info("Completed successfully")


if __name__ == "__main__":
    main()
