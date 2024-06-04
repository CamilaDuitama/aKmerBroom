from __future__ import division

import os.path
import time

from Bio import SeqIO
from scripts import kmers
from pybloomfilter import BloomFilter
import logging

logger = logging.getLogger(__name__)


def getbloomFilter(bf, bf_capacity, present_kmers_file, kmer_size, output):
    """
    Function to build and return a Bloom filter of the kmers that must be removed. Loads the
    pre-built bloom filter if nothing is given by user.
    bf: bool, true if no argument was given.
    bf_capacity: size of the bloom filter to build.
    present_kmers: path to file of kmers in text format.
    kmer_size: int.
    """
    if os.path.isfile(bf):
        logger.info("Opening Bloom Filter of present kmers")
        present_kmers_bf = BloomFilter.open(bf)
        logger.info("Done")
    else:
        logger.info("Need to make Bloom Filter of present k-mers")
        bf_filename = output + "/present_kmers.bloom"
        present_kmers_bf = BloomFilter(bf_capacity, .001, bf_filename)

        if os.path.isfile(present_kmers_file): # if present kmers txt file exist
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


def classify_reads(fastq_file, bloom_filter, kmer_size, n_consecutive_matches, output, shared_anchor_kmer_set):
    """
    Function that classifies sequences as ancient or present-day. Writes Ancient sequences in the output (fasta file format).
    :parameters:
    fastq_file : path to a fastq file
    bloom_filter: path to ptythonbloomfiltermmap3 filter
    kmer_size: int .
    n_consecutive_matches: int, number of match between kmers of the current read and bloom filter to set an anchor.
    shared_anchor_kmer_set: set of anchor kmers for reads, shared between processes
    :returns: path to directory were output files will be written.
    """
    logger.info(f"Classifying reads for input file {fastq_file}")
    print(f"Classifying reads for input file {fastq_file}", flush=True)
    annotated_reads_file = open(output+"/"+fastq_file.split("/")[-1].rstrip(".fastq")+"_annotated_reads.fastq", "w")
    read_count = 0
    anchor_kmer_set = set()

    for record in SeqIO.parse(fastq_file, "fastq"):
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
            if kmer in bloom_filter or reverse[length - kmer_size - i:length - i] in bloom_filter:
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
    shared_anchor_kmer_set.update(anchor_kmer_set)

def classify_reads_using_anchor_kmers(input_file, anchor_kmer_set_test, kmer_size, anchor_proportion_cutoff, output):
    """
    anchor_kmer_set: set of kmers.
    kmers_size: int.
    kmer_size: int.
    anchor_proportion_cutoff: float, percentage of kmers from a read that are from contamination.
    output: path to directory were output files will be written.
    """
    ip_reads_file = output +"/" + input_file.split("/")[-1].rstrip(".fastq") + "_annotated_reads.fastq"
    op_read_file_decontam = open(output +"/" + input_file.split("/")[-1].rstrip("annotated_reads.fastq") + "_decontaminated.fastq", "w")
    op_read_file_contam = open(output + "/" + input_file.split("/")[-1].rstrip("annotated_reads.fastq") + "_contamination.fastq", "w")
    print(f"Classifying reads for input file {ip_reads_file}, final round.", flush=True)
    anchor_kmer_set = anchor_kmer_set_test.copy()
    read_count = 0
    for record in SeqIO.parse(ip_reads_file, "fastq"):
        score = record.letter_annotations["phred_quality"]
        read_count += 1
        if read_count % 100000 == 0:
            print(f"No. of reads seen so far for file {ip_reads_file}: {read_count}")
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
            SeqIO.write(new_record, op_read_file_decontam, "fastq")
        else:
            SeqIO.write(record, op_read_file_contam, "fastq")

    op_read_file_decontam.close()
    op_read_file_contam.close()