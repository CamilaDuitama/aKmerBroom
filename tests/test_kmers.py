#!/usr/bin/env python3
"""
Unit tests for aKmerBroom
"""

import unittest
import sys
import os

# Add the parent directory to the path so we can import the modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import kmers


class TestKmersModule(unittest.TestCase):
    """Test the kmers utility functions"""
    
    def test_reverse_complement_basic(self):
        """Test basic reverse complement functionality"""
        self.assertEqual(kmers.reverse_complement("ATCG"), "CGAT")
        self.assertEqual(kmers.reverse_complement("AAAA"), "TTTT")
        self.assertEqual(kmers.reverse_complement("GCGC"), "GCGC")
    
    def test_reverse_complement_lowercase(self):
        """Test reverse complement with lowercase input"""
        self.assertEqual(kmers.reverse_complement("atcg"), "cgat")
        self.assertEqual(kmers.reverse_complement("AtCg"), "cGaT")  # Mixed case gets transformed
    
    def test_reverse_complement_with_n(self):
        """Test reverse complement with N bases"""
        self.assertEqual(kmers.reverse_complement("ATCGN"), "NCGAT")
        self.assertEqual(kmers.reverse_complement("NNNNN"), "NNNNN")
    
    def test_reverse_complement_empty(self):
        """Test reverse complement with empty string"""
        self.assertEqual(kmers.reverse_complement(""), "")
    
    def test_reverse_complement_single_base(self):
        """Test reverse complement with single bases"""
        self.assertEqual(kmers.reverse_complement("A"), "T")
        self.assertEqual(kmers.reverse_complement("T"), "A")
        self.assertEqual(kmers.reverse_complement("G"), "C")
        self.assertEqual(kmers.reverse_complement("C"), "G")
    
    def test_valid_kmer_format_correct(self):
        """Test valid k-mer format validation with correct input"""
        # Should not raise any exception
        try:
            kmers.test_valid_kmer_format("ATCGATCGATCGATCGATCGATCGATCGATC", 31)
            kmers.test_valid_kmer_format("ATCG", 4)
        except SystemExit:
            self.fail("test_valid_kmer_format raised SystemExit unexpectedly")
    
    def test_valid_kmer_format_incorrect_length(self):
        """Test valid k-mer format validation with incorrect length"""
        # This should call sys.exit via exit_gracefully()
        with self.assertRaises(SystemExit):
            kmers.test_valid_kmer_format("ATCG", 5)  # k-mer is 4, but expecting 5
        
        with self.assertRaises(SystemExit):
            kmers.test_valid_kmer_format("ATCGATCG", 4)  # k-mer is 8, but expecting 4


class TestKmerExtraction(unittest.TestCase):
    """Test k-mer extraction logic (we'll need to extract this from classify_reads)"""
    
    def test_kmer_extraction_basic(self):
        """Test basic k-mer extraction"""
        sequence = "ATCGATCG"
        k = 3
        expected_kmers = ["ATC", "TCG", "CGA", "GAT", "ATC", "TCG"]
        
        # Extract k-mers manually (simulating the logic from classify_reads)
        extracted_kmers = []
        for i in range(len(sequence) - k + 1):
            extracted_kmers.append(sequence[i:i + k])
        
        self.assertEqual(extracted_kmers, expected_kmers)
    
    def test_kmer_extraction_edge_cases(self):
        """Test k-mer extraction edge cases"""
        # Sequence shorter than k
        sequence = "AT"
        k = 5
        extracted_kmers = []
        for i in range(len(sequence) - k + 1):
            extracted_kmers.append(sequence[i:i + k])
        self.assertEqual(extracted_kmers, [])  # Should be empty
        
        # Sequence equal to k
        sequence = "ATCG"
        k = 4
        extracted_kmers = []
        for i in range(len(sequence) - k + 1):
            extracted_kmers.append(sequence[i:i + k])
        self.assertEqual(extracted_kmers, ["ATCG"])


if __name__ == '__main__':
    unittest.main()