from __future__ import division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scripts import kmers
from pybloomfilter import BloomFilter


def getbloomFilter(bf, bf_capacity, ancient_kmers, kmer_size):
    if bf:
        print("Opening Bloom Filter of ancient kmers")
        ancient_kmers_bf = BloomFilter.open("data/ancient_kmers.bloom")
        print("Done")
    else:
        print("Need to make Bloom Filter of ancient k-mers")
        bf_filename = "data/ancient_kmers.bloom"
        ancient_kmers_bf = BloomFilter(bf_capacity, .001, bf_filename)

        if ancient_kmers: # if ancient kmers txt file exist
            ancient_kmers_file = "data/ancient_kmers"
            with open(ancient_kmers_file, 'r') as ac_kmers:
                first_line = ac_kmers.readline().rstrip()
                kmers.test_valid_kmer_format(first_line, kmer_size)
                ac_kmers.seek(0)
                for line in ac_kmers:
                    # ancient_kmers_bf.add(line[:kmer_size])
                    ancient_kmers_bf.add(line.rstrip())

        else:
            print("Please provide an ancient kmer set")
            kmers.exit_gracefully()
        print("Done creating bloom filter")
    return ancient_kmers_bf


def classify_reads(bf, bf_capacity, ancient_kmers, kmer_size, ancient_proportion_cutoff):
    ancient_kmers_bf = getbloomFilter(bf, bf_capacity, ancient_kmers, kmer_size)
    unknown_reads_file = "data/unknown_reads.fastq"
    annotated_reads_file = open("output/annotated_reads.fasta", "w")
    read_count = 0

    for record in SeqIO.parse(unknown_reads_file, "fastq"):
        score = record.letter_annotations["phred_quality"]
        matches = []
        read_count += 1
        if read_count % 100000 == 0:
            print("No. of reads seen so far: ", read_count)
            #print("Current read ID is : ", record.id)

        to_kmerize_fwd = str(record.seq).upper()
        length = len(to_kmerize_fwd)
        reverse = kmers.reverse_complement(to_kmerize_fwd)
        count_of_all_kmers_in_this_read = 0
        count_of_ancient_kmers_in_this_read = 0
        curr_kmer_abundances = []
        for i in range(0, length - kmer_size + 1):
            kmer = to_kmerize_fwd[i:i + kmer_size]
            count_of_all_kmers_in_this_read += 1
            if kmer in ancient_kmers_bf or reverse[length - kmer_size - i:length - i] in ancient_kmers_bf:
                count_of_ancient_kmers_in_this_read += 1
                matches.append(0)    # match to ancient kmer found, represented by '!'
            else:
                matches.append(12)   # match to ancient kmer not found, represented by '-'
        for i in range(kmer_size-1):
            matches.append(12)       # need to add an extra k-1 empty matches since kmers are now all out

        try:
            proportion = round((count_of_ancient_kmers_in_this_read / count_of_all_kmers_in_this_read), 2)
        except ZeroDivisionError:
            proportion = 0

        new_record = SeqIO.SeqRecord(seq=record.seq,
                                     id=record.id,
                                     description=record.id + " " + str(length) + " " + str(proportion),
                                     letter_annotations={'phred_quality': matches},
                                     )
        SeqIO.write(new_record, annotated_reads_file, "fastq")

    annotated_reads_file.close()
    return

