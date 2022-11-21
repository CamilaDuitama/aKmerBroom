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

        #print("Total No. of k-mers from this read: ", count_of_all_kmers_in_this_read)
        #print("No. of ancient k-mers: ", count_of_ancient_kmers_in_this_read)

        try:
            proportion = round((count_of_ancient_kmers_in_this_read / count_of_all_kmers_in_this_read), 2)
        except ZeroDivisionError:
            proportion = 0
        #print("Ancient proportion is: ", proportion)

        # Don't apply cutoff at this stage
        #if proportion >= ancient_proportion_cutoff:
        SeqIO.write(SeqRecord(record.seq, id=record.id, description=str(length) + " " + str(proportion)),
                    annotated_reads_file, "fasta")
    annotated_reads_file.close()
    return

