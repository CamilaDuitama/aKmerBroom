from __future__ import division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scripts import kmers
from pybloomfilter import BloomFilter
import logging

logger = logging.getLogger(__name__)

def getbloomFilter(bf, bf_capacity, present_kmers, kmer_size):
    """
    Function to build and return a Bloom filter of the kmers that must be removed. Loads the
    pre-built bloom filter if nothing is given by user.
    bf: bool, true if no argument was given.
    bf_capacity: size of the bloom filter to build.
    present_kmers: path to file of kmers in text format.
    kmer_size: int.
    """
    if bf:
        logger.info("Opening Bloom Filter of present kmers")
        present_kmers_bf = BloomFilter.open("data/present_kmers.bloom")
        logger.info("Done")
    else:
        logger.info("Need to make Bloom Filter of present k-mers")
        bf_filename = "data/present_kmers.bloom"
        present_kmers_bf = BloomFilter(bf_capacity, .001, bf_filename)

        if present_kmers: # if present kmers txt file exist
            present_kmers_file = "data/present_kmers"
            with open(present_kmers_file, 'r') as ac_kmers:
                first_line = ac_kmers.readline().rstrip()
                kmers.test_valid_kmer_format(first_line, kmer_size)
                ac_kmers.seek(0)
                for line in ac_kmers:
                    # present_kmers_bf.add(line[:kmer_size])
                    present_kmers_bf.add(line.rstrip())

        else:
            logger.error("Please provide an present kmer set")
            kmers.exit_gracefully()
        logger.info("Done creating bloom filter")
    return present_kmers_bf

#TODO: modifier la fonction ci-dessous pour virer les sequences 'annotees' comme contaminantes
def classify_reads(bf, bf_capacity, present_kmers, kmer_size, n_consecutive_matches, output):
    """
    Function that classifies seqeunces as ancient or present-day. Writes Ancient sequences in the output (fasta file format).
    bf: bool, true if no argument was given.
    bf_capacity: size of the bloom filter to build.
    present_kmers: path to file of kmers in text format.
    kmer_size: int .
    n_consecutive_matches: int, number of match between kmers of the current read and bloom filter to set an anchor.
    output: path to directory were output files will be written.
    """
    present_kmers_bf = getbloomFilter(bf, bf_capacity, present_kmers, kmer_size)
    unknown_reads_file = "data/unknown_reads.fastq"
    annotated_reads_file = open(output+"/annotated_reads.fastq", "w")
    read_count = 0
    anchor_kmer_set = set()

    for record in SeqIO.parse(unknown_reads_file, "fastq"):
        score = record.letter_annotations["phred_quality"]
        matches = []
        consecutive_matches_found = 0
        curr_kmer_set = set()
        read_count += 1
        if read_count % 100000 == 0:
            logger.info("No. of reads seen so far: "+str(read_count))

        to_kmerize_fwd = str(record.seq).upper()
        length = len(to_kmerize_fwd)
        reverse = kmers.reverse_complement(to_kmerize_fwd)
        count_of_all_kmers_in_this_read = 0
        count_of_present_kmers_in_this_read = 0
        for i in range(0, length - kmer_size + 1):
            kmer = to_kmerize_fwd[i:i + kmer_size]
            curr_kmer_set.add(kmer)
            count_of_all_kmers_in_this_read += 1
            if kmer in present_kmers_bf or reverse[length - kmer_size - i:length - i] in present_kmers_bf:
                count_of_present_kmers_in_this_read += 1
                matches.append(0)    # match to present kmer found, represented by '!'
            else:
                matches.append(12)   # match to present kmer not found, represented by '-'
        for i in range(kmer_size-1):
            matches.append(12)       # need to add an extra k-1 empty matches since kmers are now all out

        # test if consecutive_matches is greater than or equal to n_consecutive_matches parameter
        stringified_matches = ''.join(str(match) for match in matches)
        consecutive_matches_string = n_consecutive_matches*'0'
        if consecutive_matches_string in stringified_matches:
            consecutive_matches_found = 1

        if consecutive_matches_found:
            for kmer in curr_kmer_set:
                anchor_kmer_set.add(kmer)

        new_record = SeqIO.SeqRecord(seq=record.seq,
                                     id=record.id,
                                     description=record.id + " " + str(length) + " " + str(consecutive_matches_found),
                                     letter_annotations={'phred_quality': matches},
                                     )
        SeqIO.write(new_record, annotated_reads_file, "fastq")

    annotated_reads_file.close()
    return anchor_kmer_set


#TODO: changer les noms de merde des fichiers de sortie.
def classify_reads_using_anchor_kmers(anchor_kmer_set, kmer_size, anchor_proportion_cutoff,output):
    """
    anchor_kmer_set: set of kmers.
    kmers_size: int.
    kmer_size: int.
    anchor_proportion_cutoff: float, percentage of kmers from a read that are from contamination.
    output: path to directory were output files will be written.
    """
    ip_reads_file = output+"/annotated_reads.fastq"
    op_reads_file = open(output+"/annotated_reads_with_anchor_kmers.fastq", "w")
    read_count = 0
    for record in SeqIO.parse(ip_reads_file, "fastq"):
        score = record.letter_annotations["phred_quality"]
        read_count += 1
        if read_count % 100000 == 0:
            logger.info("No. of reads seen so far: " + str(read_count))

        to_kmerize_fwd = str(record.seq).upper()
        length = len(to_kmerize_fwd)
        reverse = kmers.reverse_complement(to_kmerize_fwd)
        count_of_all_kmers_in_this_read = 0
        count_of_anchor_kmers_in_this_read = 0
        for i in range(0, length - kmer_size + 1):
            kmer = to_kmerize_fwd[i:i + kmer_size]
            count_of_all_kmers_in_this_read += 1
            if kmer in anchor_kmer_set or reverse[length - kmer_size - i:length - i] in anchor_kmer_set:
                count_of_anchor_kmers_in_this_read += 1

        # compute anchor_proportion
        try:
            #print(count_of_anchor_kmers_in_this_read)
            anchor_proportion = round((count_of_anchor_kmers_in_this_read / count_of_all_kmers_in_this_read), 2)
        except ZeroDivisionError:
            anchor_proportion = 0

        # select present reads if they meet the proportion cutoff
        if anchor_proportion >= anchor_proportion_cutoff:
            new_record = SeqIO.SeqRecord(seq=record.seq,
                                         id=record.id,
                                         description=record.description + " " + str(anchor_proportion),
                                         letter_annotations={'phred_quality': score},
                                         )
            SeqIO.write(new_record, op_reads_file, "fastq")

    op_reads_file.close()
    return
