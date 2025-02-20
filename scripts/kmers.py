import sys
import logging

logger = logging.getLogger(__name__)

def exit_gracefully():
    sys.exit("^Program exited with error, please see the message written in the file aKmerBroom.log located in the current working directory.^")


def reverse_complement(seq: str) -> str:
        complement_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return seq[::-1].translate(complement_map)


def test_valid_kmer_format(sample_line: str, kmer_size: int):
   if len(sample_line) != kmer_size:
       logger.error("Format of kmer file does not match kmer size specified")
       exit_gracefully()
