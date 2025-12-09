from __future__ import division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scripts import kmers
from pybloomfilter import BloomFilter
import logging
import os
import time
import sys

logger = logging.getLogger(__name__)

def estimate_total_reads(filename):
    """Estimate total number of reads in FASTQ file by counting lines"""
    try:
        with open(filename, 'r') as f:
            lines = sum(1 for _ in f)
        return lines // 4  # FASTQ has 4 lines per read
    except IOError:
        return None

def format_time(seconds):
    """Format time in human readable format"""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{int(seconds//60)}m {int(seconds%60)}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}h {minutes}m"

def log_progress(read_count, total_reads, start_time, logger, interval=50000):
    """Log progress with percentage and ETA"""
    if read_count % interval == 0 or total_reads and read_count == total_reads:
        elapsed = time.time() - start_time
        
        if total_reads:
            percentage = (read_count / total_reads) * 100
            if read_count > 0 and elapsed > 0:
                reads_per_sec = read_count / elapsed
                remaining_reads = total_reads - read_count
                eta_seconds = remaining_reads / reads_per_sec if reads_per_sec > 0 else 0
                eta_str = f", ETA: {format_time(eta_seconds)}" if eta_seconds > 5 else ""
            else:
                eta_str = ""
            logger.info(f"Progress: {read_count:,}/{total_reads:,} reads ({percentage:.1f}%) - Elapsed: {format_time(elapsed)}{eta_str}")
        else:
            logger.info(f"Progress: {read_count:,} reads processed - Elapsed: {format_time(elapsed)}")

def getbloomFilter(bf, bf_capacity, ancient_kmers, kmer_size):
    if bf:
        logger.info("Opening Bloom Filter of ancient kmers")
        ancient_kmers_bf = BloomFilter.open("data/ancient_kmers.bloom")
        logger.info("Done")
    else:
        logger.info("Need to make Bloom Filter of ancient k-mers")
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
            logger.error("Please provide an ancient kmer set")
            kmers.exit_gracefully()
        logger.info("Done creating bloom filter")
    return ancient_kmers_bf


def classify_reads(bf, bf_capacity, ancient_kmers, kmer_size, n_consecutive_matches, output, input_file="data/unknown_reads.fastq", output_prefix=""):
    # Input validation
    if not os.path.exists(input_file):
        logger.error(f"Input file does not exist: {input_file}")
        kmers.exit_gracefully()
    
    if not os.path.isfile(input_file):
        logger.error(f"Input path is not a file: {input_file}")
        kmers.exit_gracefully()
    
    if kmer_size <= 0:
        logger.error(f"Invalid k-mer size: {kmer_size}. Must be positive integer.")
        kmers.exit_gracefully()
    
    if n_consecutive_matches <= 0:
        logger.error(f"Invalid consecutive matches: {n_consecutive_matches}. Must be positive integer.")
        kmers.exit_gracefully()
    
    # Validate output directory exists or create it
    if not os.path.exists(output):
        logger.error(f"Output directory does not exist: {output}")
        kmers.exit_gracefully()
    
    ancient_kmers_bf = getbloomFilter(bf, bf_capacity, ancient_kmers, kmer_size)
    unknown_reads_file = input_file
    
    # Create output filename with prefix
    if output_prefix:
        annotated_filename = f"{output}/{output_prefix}_annotated_reads.fastq"
    else:
        annotated_filename = f"{output}/annotated_reads.fastq"
    
    read_count = 0
    anchor_kmer_set = set()
    anchor_reads_count = 0
    
    # Estimate total reads for progress reporting
    total_reads = estimate_total_reads(unknown_reads_file)
    if total_reads:
        logger.info(f"Estimated total reads to process: {total_reads:,}")
    else:
        logger.info("Starting read processing (unable to estimate total)")
    
    start_time = time.time()
    
    try:
        # Use context manager for proper file handling
        with open(annotated_filename, "w") as annotated_reads_file:
            for record in SeqIO.parse(unknown_reads_file, "fastq"):
                score = record.letter_annotations["phred_quality"]
                matches = []
                consecutive_matches_found = 0
                curr_kmer_set = set()
                read_count += 1
                
                # Progress reporting
                log_progress(read_count, total_reads, start_time, logger)

                to_kmerize_fwd = str(record.seq).upper()
                length = len(to_kmerize_fwd)
                
                # Skip reads that are too short for k-mer extraction
                if length < kmer_size:
                    logger.warning(f"Read {record.id} is too short ({length} bp) for k-mer size {kmer_size}. Skipping.")
                    continue
                
                reverse = kmers.reverse_complement(to_kmerize_fwd)
                count_of_all_kmers_in_this_read = 0
                count_of_ancient_kmers_in_this_read = 0
                
                for i in range(0, length - kmer_size + 1):
                    kmer = to_kmerize_fwd[i:i + kmer_size]
                    curr_kmer_set.add(kmer)
                    count_of_all_kmers_in_this_read += 1
                    if kmer in ancient_kmers_bf or reverse[length - kmer_size - i:length - i] in ancient_kmers_bf:
                        count_of_ancient_kmers_in_this_read += 1
                        matches.append(0)    # match to ancient kmer found, represented by '0'
                    else:
                        matches.append(12)   # match to ancient kmer not found, represented by '12'
                        
                # Pad matches to account for k-mer positions
                for i in range(kmer_size-1):
                    matches.append(12)       # need to add an extra k-1 empty matches since kmers are now all out

                # test if consecutive_matches is greater than or equal to n_consecutive_matches parameter
                stringified_matches = ''.join(str(match) for match in matches)
                consecutive_matches_string = n_consecutive_matches*'0'
                if consecutive_matches_string in stringified_matches:
                    consecutive_matches_found = 1

                if consecutive_matches_found:
                    anchor_reads_count += 1
                    for kmer in curr_kmer_set:
                        anchor_kmer_set.add(kmer)

                new_record = SeqIO.SeqRecord(seq=record.seq,
                                             id=record.id,
                                             description=record.id + " " + str(length) + " " + str(consecutive_matches_found),
                                             letter_annotations={'phred_quality': matches},
                                             )
                SeqIO.write(new_record, annotated_reads_file, "fastq")
                
        # Final progress and summary
        elapsed_time = time.time() - start_time
        reads_per_sec = read_count / elapsed_time if elapsed_time > 0 else 0
        
        logger.info(f"\\n=== STEP 1 SUMMARY: Anchor Read Identification ===")
        logger.info(f"Total reads processed: {read_count:,}")
        logger.info(f"Anchor reads found: {anchor_reads_count:,} ({anchor_reads_count/read_count*100:.2f}%)")
        logger.info(f"Unique anchor k-mers: {len(anchor_kmer_set):,}")
        logger.info(f"Processing time: {format_time(elapsed_time)}")
        logger.info(f"Processing speed: {reads_per_sec:.1f} reads/second")
        
    except IOError as e:
        logger.error(f"Error writing to output file {annotated_filename}: {str(e)}")
        kmers.exit_gracefully()
    except Exception as e:
        logger.error(f"Unexpected error in classify_reads: {str(e)}")
        kmers.exit_gracefully()

    return anchor_kmer_set


def classify_reads_using_anchor_kmers(anchor_kmer_set, kmer_size, anchor_proportion_cutoff, output, output_prefix=""):
    # Input validation
    if not anchor_kmer_set:
        logger.warning("No anchor k-mers provided. Output file will be empty.")
    
    if kmer_size <= 0:
        logger.error(f"Invalid k-mer size: {kmer_size}. Must be positive integer.")
        kmers.exit_gracefully()
    
    if not (0 <= anchor_proportion_cutoff <= 1):
        logger.error(f"Invalid anchor proportion cutoff: {anchor_proportion_cutoff}. Must be between 0 and 1.")
        kmers.exit_gracefully()
    
    # Create input filename with prefix
    if output_prefix:
        ip_reads_file = f"{output}/{output_prefix}_annotated_reads.fastq"
        op_filename = f"{output}/{output_prefix}_annotated_reads_with_anchor_kmers.fastq"
    else:
        ip_reads_file = f"{output}/annotated_reads.fastq"
        op_filename = f"{output}/annotated_reads_with_anchor_kmers.fastq"
    
    # Validate input file exists
    if not os.path.exists(ip_reads_file):
        logger.error(f"Input file from previous step does not exist: {ip_reads_file}")
        kmers.exit_gracefully()
    
    read_count = 0
    ancient_reads_count = 0
    
    # Estimate total reads for progress reporting
    total_reads = estimate_total_reads(ip_reads_file)
    if total_reads:
        logger.info(f"Starting final classification on {total_reads:,} reads")
    else:
        logger.info("Starting final classification (unable to estimate total)")
    
    start_time = time.time()
    
    try:
        with open(op_filename, "w") as op_reads_file:
            for record in SeqIO.parse(ip_reads_file, "fastq"):
                score = record.letter_annotations["phred_quality"]
                read_count += 1
                
                # Progress reporting
                log_progress(read_count, total_reads, start_time, logger)

                to_kmerize_fwd = str(record.seq).upper()
                length = len(to_kmerize_fwd)
                
                # Skip reads that are too short for k-mer extraction
                if length < kmer_size:
                    logger.warning(f"Read {record.id} is too short ({length} bp) for k-mer size {kmer_size}. Skipping.")
                    continue
                    
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
                    logger.warning(f"Read {record.id} has no k-mers (length: {length}, kmer_size: {kmer_size})")

                # select ancient reads if they meet the proportion cutoff
                if anchor_proportion >= anchor_proportion_cutoff:
                    ancient_reads_count += 1
                    new_record = SeqIO.SeqRecord(seq=record.seq,
                                                 id=record.id,
                                                 description=record.description + " " + str(anchor_proportion),
                                                 letter_annotations={'phred_quality': score},
                                                 )
                    SeqIO.write(new_record, op_reads_file, "fastq")
                    
        # Final progress and summary
        elapsed_time = time.time() - start_time
        reads_per_sec = read_count / elapsed_time if elapsed_time > 0 else 0
        
        logger.info(f"\\n=== STEP 2 SUMMARY: Final Ancient Read Classification ===")
        logger.info(f"Total reads processed: {read_count:,}")
        logger.info(f"Ancient reads identified: {ancient_reads_count:,} ({ancient_reads_count/read_count*100:.2f}%)")
        logger.info(f"Classification cutoff: {anchor_proportion_cutoff}")
        logger.info(f"Processing time: {format_time(elapsed_time)}")
        logger.info(f"Processing speed: {reads_per_sec:.1f} reads/second")
        logger.info(f"Output written to: {op_filename}")
        
    except IOError as e:
        logger.error(f"Error writing to output file {op_filename}: {str(e)}")
        kmers.exit_gracefully()
    except Exception as e:
        logger.error(f"Unexpected error in classify_reads_using_anchor_kmers: {str(e)}")
        kmers.exit_gracefully()
        
    return