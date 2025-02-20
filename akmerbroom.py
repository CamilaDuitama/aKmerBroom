import argparse
import os
import logging
import sys

from scripts import classify_reads, kmers
from multiprocessing import Value, Process, Pool
from multiprocessing.managers import BaseManager
from ctypes import c_wchar_p


#TODO: make the single or multi-sample decontamination an option (through a (non-)shared object of anchor_kmers)
#TODO: arg pour dire si la ref est ancienne ou presente et en consequence retirer ce qui est reconnu ou le garder: DONE
#TODO: verifier que le multi-threading est tjrs pas possible sur la seconde passe

def main():
    """
    Main function.
    """
    parser = argparse.ArgumentParser(
        description='This program uses a reference (either a bloom filter or a kmer text file) to recognise targeted DNA \
        (ancient or modern) reads to separate ancient DNA from modern DNA. aDNA will be stored in a "decontaminated" file \
        and the modern DNA in the "contamination" for each sample.',
        usage='%(prog)s [options]',
        epilog="EX: python3 akmerbroom.py -i $(find ./tests/*.fastq) -o ./output -t 2 --kmers_set kmers.txt")
    parser.add_argument('--bloom', help='Used if a BloomFilter is provided (defaults to False)',
                        type=str, action='store', required=False, default="")
    parser.add_argument("--bloom_capacity", type=int, help="If a BloomFilter is not provided, \
                       this sets the capacity of the bloom filter. This should be greater than the number of distinct \
                       kmers in the input file. Default to 2 billion.")
    parser.add_argument('--kmers_set', help='Used if a kmers set is provided (defaults to False).',
                        type=str, action='store', required=False, default="")
    parser.add_argument('-k', '--kmer_size', type=int, help='Set kmer size (defaults to 31)', required=False,
                        default=31)
    parser.add_argument('--n_consec_matches', type=int, help="Set number of consec matches to classify \
                        read as anchor read, (defaults to 2).", default=2, required=False)
    parser.add_argument('--anchor_proportion_cutoff', help="Set anchor kmer proportion, \
                        above which a read is classified as modern/ancient (defaults to 0.5)", default=0.5, type=float,
                        required=False)
    parser.add_argument('-i', '--input', help="Path to input file(s), space-separated.", required=True, nargs='+')
    parser.add_argument("-o", "--output",
                        help="Path to output folder, where you want aKmerBroom to write the results.", required=True)
    parser.add_argument("-t", "--threads",
                        help="WARNING: right now, not used. Sorry, async is a pain. Number of threads to use, default to 1.",
                        default=1, type=int)
    parser.add_argument("-s", "--single", action="store_true", help="Flag to decontaminate samples independently \
                        instead of pooling k-mers from multi-samples for decontamination.")
    parser.add_argument("-m", "--modern", action = "store_true", help="Flag to indicate that reference \
                        is modern DNA (defaults to False).")

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
        user_response = ""
        while user_response.lower() != "y" or user_response.lower() != "n":
            user_response = input("Continue? (y/n) ")
            if user_response.lower() == "n":
                logger.warning("Output directory NOT overwritten.")
                kmers.exit_gracefully()
            else:
                logger.warning("Output directory overwritten.")
                break
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

    # checks existence of modern kmers in any form
    if not os.path.isfile(args['bloom']) and not os.path.isfile(args['kmers_set']):
        logger.error("Please provide a kmer set or a bloom filter.")
        kmers.exit_gracefully()
    elif args['bloom_capacity']:
        try:
            bf_capacity = int(args['bloom_capacity'])
        except ValueError:
            logger.error("bloom_capacity provided is not an integer.")
            kmers.exit_gracefully()
    else:
        bf_capacity = 1000 * 1000 * 1000 * 2

    logger.info("Started...")
    # declare defaults
    logger.info("Shortlisting modern reads")
    # builds the bloom filter once for all
    print("Getting the Bloom Filter ready ...")
    kmers_bf = classify_reads.getbloomFilter(args["bloom"], bf_capacity, args["kmers_set"],
                                                     args["kmer_size"], args["output"])

    # create Values for shared memory for the parameters
    shared_k_size = Value('i', args['kmer_size'])
    shared_n_consec_matches = Value('i', args['n_consec_matches'])
    shared_output = Value(c_wchar_p, args["output"])
    shared_anchor_proportion_cutoff = Value('f', args["anchor_proportion_cutoff"])
    shared_modern_DNA_flag = Value('i', args['modern'])

    # creates a custom manager for a set in shared memory
    class CustomManager(BaseManager):
        # nothing
        pass

    CustomManager.register('shared_set', set)
    manager = CustomManager()
    manager.start()
    shared_kmer_set = manager.shared_set()

    #TODO; modifier ce bout de code pour permettre de faire en single-sample
    #TODO: soit je fais un if ici pour passer en multi-threading soit je le mets dans la fonction : la modification necessaire serait de copier le bloom et de pas le partager a la fin, mais embarquer sur la suite avec la copie augmentee sans s'arreter pour join

    if args["single"]:
        # processes = [Process(target=classify_reads.classify_reads, args=(
        #     i, kmers_bf, shared_k_size.value, shared_n_consec_matches.value, shared_output.value, shared_kmer_set))
        #              for i in args["input"]]
        # for process in processes:
        #     process.start()
        # for process in processes:
        #     process.join()
        # for process in processes:
        #     process.terminate()
        #
        # with Pool(processes=args["threads"]) as pool:
        #     for i in args["input"]:
        #         pool.apply_async(classify_reads.classify_reads_using_anchor_kmers, args=(i, shared_kmer_set,
        #                                                                                  shared_k_size.value,
        #                                                                                  shared_anchor_proportion_cutoff.value,
        #                                                                                  shared_output.value,
        #                                                                                  shared_modern_DNA_flag.value))
        #     pool.close()
        #     pool.join()
        for i in args["input"]:
            anchor_kmers_set_single = set()
            classify_reads.classify_reads(i, kmers_bf, shared_k_size.value, shared_n_consec_matches.value, shared_output.value, anchor_kmers_set_single)
            classify_reads.classify_reads_using_anchor_kmers(i,anchor_kmers_set_single, shared_k_size.value,
                                                                                shared_anchor_proportion_cutoff.value,
                                                                                shared_output.value,
                                                                                shared_modern_DNA_flag.value)

    else :
        processes = [Process(target=classify_reads.classify_reads, args=(
            i, kmers_bf, shared_k_size.value, shared_n_consec_matches.value, shared_output.value, shared_kmer_set))
                     for i in args["input"]]
        for process in processes:
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
                                                                                         shared_output.value,
                                                                                         shared_modern_DNA_flag.value))
            pool.close()
            pool.join()






    # does not work and crashes everything ????????
    # with Pool(processes=args["threads"]) as pool:
    #     for i in args["input"]:
    #         pool.apply_async(classify_reads.classify_reads, args=(i, kmers_bf, shared_k_size.value,
    #                                                               shared_n_consec_matches.value,
    #                                                               shared_output.value,
    #                                                               shared_kmer_set))
    #     pool.close()
    #     pool.join()



    logger.info("Completed successfully")
    print("Done !")


if __name__ == "__main__":
    main()
