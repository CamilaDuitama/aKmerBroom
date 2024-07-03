import argparse
import os
import logging
from scripts import classify_reads, kmers
from multiprocessing import Value, Process, Pool
from multiprocessing.managers import BaseManager
from ctypes import c_wchar_p


def main():
    """
    Main function.
    """
    parser = argparse.ArgumentParser(
        description='This program finds removes present-day reads to only keep the ancient DNA.',
        epilog="Ex: python3 akmerbroom.py -i test1.fastq test2.fastq -t 2 --present_kmers_set kmers.txt\n"
               "EX: python3 akmerbroom.py -i $(find ./*.fastq) -t 2 --present_kmers_set kmers.txt")
    parser.add_argument('--present_bloom', help='Use if present BloomFilter provided (defaults to False)',
                        type=str, action='store', required=False, default="")
    parser.add_argument("--present_bloom_capacity", type=int, help="If present BloomFilter is not provided, \
        This sets the capacity of the bloom filter. This should be greater than the number of distinct kmers \
        in present kmers input file")
    parser.add_argument('--present_kmers_set', help='Use if present kmers set provided (defaults to False)',
                        type=str, action='store', required=False, default="")
    parser.add_argument('-k', '--kmer_size', type=int, help='Set kmer size (defaults to 31)', required=False,
                        default=31)
    parser.add_argument('--n_consec_matches', type=int, help="Set number of consec matches to classify" \
                                                             " read as anchor read, (defaults to 2)", default=2,
                        required=False)
    parser.add_argument('--anchor_proportion_cutoff', help="Set anchor kmer proportion, \
        above which a read is classified as present (defaults to 0.5)", default=0.5, type=float,
                        required=False)
    parser.add_argument('-i', '--input', help="Path to input file(s), space-separated", required=True, nargs='+')
    parser.add_argument("-o", "--output",
                        help="Path to output folder, where you want aKmerBroom to write the results.", required=True,
                        default="output")
    parser.add_argument("-t", "--threads",
                        help="WARNING: right now, not used. Sorry, async is a pain. Number of threads to use",
                        default=1, type=int)

    args = vars(parser.parse_args())

    bf_capacity = 0

    # Create and configure logger
    logging.basicConfig(filename="aKmerBroom.log",
                        format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG, filemode='w')

    # Creating an object
    logger = logging.getLogger()

    #check output folder does not exist else create folder
    if os.path.isdir(args["output"]):
        print(f"Output directory {args['output']} already exists. Content might be overwritten.")
        logger.warning("Output directory already exists. Content might be overwritten.")
    else:
        output = args["output"]
        os.mkdir(output)

    #check if input files exist
    for i in args["input"]:
        if not os.path.exists(i):
            logger.error(f"Input {i} does not exist.")
            kmers.exit_gracefully()

    if args["threads"] < 1:
        logger.error("The number of threads is not valid.")
        kmers.exit_gracefully()

    # checks existence of present kmers in any form
    if not os.path.isfile(args['present_bloom']) and not os.path.isfile(args['present_kmers_set']):
        logger.error("Please provide a kmer set or a bloom filter.")
        kmers.exit_gracefully()
    elif args['present_bloom_capacity']:
        try:
            bf_capacity = int(args['present_bloom_capacity'])
        except ValueError:
            logger.error("Present_bloom_capacity provided is not an integer.")
            kmers.exit_gracefully()
    else:
        bf_capacity = 1000 * 1000 * 1000 * 2

    logger.info("Started...")
    # declare defaults
    logger.info("Shortlisting present reads")
    # builds the bloom filter once for all
    print("Getting the Bloom Filter ready ...")
    present_kmers_bf = classify_reads.getbloomFilter(args["present_bloom"], bf_capacity, args["present_kmers_set"],
                                                     args["kmer_size"], args["output"])

    # create Values for shared memory for the parameters
    shared_k_size = Value('i', args['kmer_size'])
    shared_n_consec_matches = Value('i', args['n_consec_matches'])
    shared_output = Value(c_wchar_p, args["output"])
    shared_anchor_proportion_cutoff = Value('f', args["anchor_proportion_cutoff"])

    # creates a custom manager for a set in shared memory
    class CustomManager(BaseManager):
        # nothing
        pass

    CustomManager.register('shared_set', set)
    manager = CustomManager()
    manager.start()
    shared_kmer_set = manager.shared_set()

    # with Pool(processes=args["threads"]) as pool:
    #     print(f"starting async")
    #     for i in args["input"]:
    #         print(i)
    #         test = pool.apply_async(classify_reads.classify_reads, args=(i, present_kmers_bf,
    #                                                               shared_k_size.value,
    #                                                               shared_n_consec_matches.value,
    #                                                               shared_output.value,
    #                                                               shared_kmer_set,))
    #     pool.close()
    #     pool.join()

    # with Pool(processes=args["threads"]) as pool:
    #     for i in args["input"]:
    #         pool.apply_async(classify_reads.classify_reads, args=(i, present_kmers_bf, shared_k_size.value, shared_n_consec_matches.value, shared_output.value, shared_kmer_set))
    #     pool.close()
    #     pool.join()

    processes = [Process(target=classify_reads.classify_reads, args=(
        i, present_kmers_bf, shared_k_size.value, shared_n_consec_matches.value, shared_output.value, shared_kmer_set))
                 for
                 i in args["input"]]
    for i in range(0, len(processes), args["threads"]):
        for process in processes[i:i + args["threads"]]:
            process.start()
    for process in processes:
        process.join()
    for process in processes:
        process.terminate()

    with Pool(processes=args["threads"]) as pool:
        for i in args["input"]:
            pool.apply_async(classify_reads.classify_reads_using_anchor_kmers, args=(i, shared_kmer_set,
                                                                                     shared_k_size.value,
                                                                                     shared_anchor_proportion_cutoff.value,
                                                                                     shared_output.value,))
        pool.close()
        pool.join()

    # processes = [Process(target=classify_reads.classify_reads_using_anchor_kmers, args=(i, shared_kmer_set, shared_k_size.value, shared_anchor_proportion_cutoff.value, shared_output.value)) for i in args["input"]]
    # for i in range(0, len(processes), args["threads"]):
    #     for process in processes[i:i + args["threads"]]:
    #         process.start()
    # for process in processes:
    #     process.join()
    # for process in processes:
    #     process.terminate()

    logger.info("Completed successfully")
    print("Done !")


if __name__ == "__main__":
    main()
