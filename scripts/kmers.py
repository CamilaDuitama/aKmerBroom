import sys


def exit_gracefully():
    sys.exit("^Program exited with error, please see the message above.^")


def reverse_complement(seq):
        complement_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return seq[::-1].translate(complement_map)


def test_valid_kmer_format(sample_line, kmer_size):
   if len(sample_line) != kmer_size:
       print("Format of kmer file does not match kmer size specified")
       exit_gracefully()

