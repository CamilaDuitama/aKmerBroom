import sys
import logging
from collections import Counter
from typing import Dict

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


def count_bases(seq: str) -> Dict[str, int]:
    """
    Count the number of bases in a sequence.
    :param seq:
    :return: dictionary of bases counts
    """
    counts = Counter(seq)
    return counts


def test_C_proportion(seq: str):
    """
    Test the proportion of C in a sequence and computes the pvalue of the wilcoxon test.
    :param seq: str
    :return: pvalue of the test or boolean if the pvalue is significant enough.
    """
    counts = count_bases(seq)
    obs_prop_C = counts['C']/ len(seq)
    est_prop_C = round(0.4*len(seq))
    #TODO: verifier la proportion attendue de C dans la sequence
    #TODO: integrer un test stat entre les deux porportions de C
    #TODO: ajout parametre de la pvaleur filtre


def count_N(seq: str) -> bool:
    counts = count_bases(seq)
    prop_N = counts['N']/ len(seq)
    if prop_N > 0.5:
        cond = True
    else:
        cond = False
    return cond
