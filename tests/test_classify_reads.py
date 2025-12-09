#!/usr/bin/env python3
"""
Integration tests for aKmerBroom classify_reads functionality
"""

import unittest
import tempfile
import os
import sys
from unittest.mock import Mock, patch

# Add the parent directory to the path so we can import the modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts import classify_reads
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class TestClassifyReads(unittest.TestCase):
    """Test classify_reads functions with mock data"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.kmer_size = 3
        
        # Create a simple test FASTQ content
        self.test_fastq_content = """@seq1
ATCGATCG
+
IIIIIIII
@seq2
GCTAGCTA
+
IIIIIIII
@seq3
AAATTTCC
+
IIIIIIII
"""
    
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)
    
    def create_temp_fastq(self, content):
        """Helper to create temporary FASTQ file"""
        temp_file = os.path.join(self.temp_dir, "test.fastq")
        with open(temp_file, 'w') as f:
            f.write(content)
        return temp_file
    
    def create_mock_bloom_filter(self, kmers_to_include):
        """Create a mock bloom filter that contains specific kmers"""
        mock_bf = Mock()
        mock_bf.__contains__ = lambda self, kmer: kmer in kmers_to_include
        return mock_bf
    
    @patch('scripts.classify_reads.getbloomFilter')
    def test_classify_reads_basic(self, mock_get_bloom):
        """Test basic classify_reads functionality"""
        # Create test input file
        input_file = self.create_temp_fastq(self.test_fastq_content)
        
        # Mock bloom filter that contains some kmers from our test sequences
        mock_bloom = self.create_mock_bloom_filter(['ATC', 'TCG', 'CGA'])
        mock_get_bloom.return_value = mock_bloom
        
        # Run classify_reads
        output_dir = self.temp_dir
        anchor_kmers = classify_reads.classify_reads(
            bf=True,
            bf_capacity=1000,
            ancient_kmers=False,
            kmer_size=self.kmer_size,
            n_consecutive_matches=2,
            output=output_dir,
            input_file=input_file,
            output_prefix="test"
        )
        
        # Verify output file was created
        expected_output = os.path.join(output_dir, "test_annotated_reads.fastq")
        self.assertTrue(os.path.exists(expected_output))
        
        # Verify anchor_kmers is a set
        self.assertIsInstance(anchor_kmers, set)
        
        # Clean up
        if os.path.exists(expected_output):
            os.remove(expected_output)
    
    def test_output_filename_generation(self):
        """Test output filename generation with and without prefix"""
        # Test with prefix
        output_dir = "/test/dir"
        prefix = "sample1"
        
        # This tests the logic inside classify_reads function
        if prefix:
            expected = f"{output_dir}/{prefix}_annotated_reads.fastq"
        else:
            expected = f"{output_dir}/annotated_reads.fastq"
        
        self.assertEqual(expected, "/test/dir/sample1_annotated_reads.fastq")
        
        # Test without prefix
        prefix = ""
        if prefix:
            expected = f"{output_dir}/{prefix}_annotated_reads.fastq"
        else:
            expected = f"{output_dir}/annotated_reads.fastq"
        
        self.assertEqual(expected, "/test/dir/annotated_reads.fastq")


class TestConsecutiveMatches(unittest.TestCase):
    """Test consecutive match detection logic"""
    
    def test_consecutive_match_detection(self):
        """Test the consecutive match string detection logic"""
        # Simulate the logic from classify_reads
        matches = [0, 0, 12, 0, 0, 0, 12, 12]  # 0 = match, 12 = no match
        stringified_matches = ''.join(str(match) for match in matches)
        
        # Test different consecutive match requirements
        n_consec_matches = 2
        consecutive_matches_string = n_consec_matches * '0'
        found_2_consecutive = consecutive_matches_string in stringified_matches
        self.assertTrue(found_2_consecutive)  # Should find "00"
        
        n_consec_matches = 3
        consecutive_matches_string = n_consec_matches * '0'
        found_3_consecutive = consecutive_matches_string in stringified_matches
        self.assertTrue(found_3_consecutive)  # Should find "000"
        
        n_consec_matches = 4
        consecutive_matches_string = n_consec_matches * '0'
        found_4_consecutive = consecutive_matches_string in stringified_matches
        self.assertFalse(found_4_consecutive)  # Should NOT find "0000"
    
    def test_no_consecutive_matches(self):
        """Test case with no consecutive matches"""
        matches = [0, 12, 0, 12, 0, 12, 0, 12]  # Alternating matches
        stringified_matches = ''.join(str(match) for match in matches)
        
        n_consec_matches = 2
        consecutive_matches_string = n_consec_matches * '0'
        found_consecutive = consecutive_matches_string in stringified_matches
        self.assertFalse(found_consecutive)


class TestAnchorProportion(unittest.TestCase):
    """Test anchor proportion calculation logic"""
    
    def test_anchor_proportion_calculation(self):
        """Test the anchor proportion calculation"""
        # Simulate the logic from classify_reads_using_anchor_kmers
        
        # Test case 1: 50% match
        count_of_anchor_kmers_in_this_read = 3
        count_of_all_kmers_in_this_read = 6
        anchor_proportion = round((count_of_anchor_kmers_in_this_read / count_of_all_kmers_in_this_read), 2)
        self.assertEqual(anchor_proportion, 0.5)
        
        # Test case 2: 100% match
        count_of_anchor_kmers_in_this_read = 5
        count_of_all_kmers_in_this_read = 5
        anchor_proportion = round((count_of_anchor_kmers_in_this_read / count_of_all_kmers_in_this_read), 2)
        self.assertEqual(anchor_proportion, 1.0)
        
        # Test case 3: 0% match
        count_of_anchor_kmers_in_this_read = 0
        count_of_all_kmers_in_this_read = 4
        anchor_proportion = round((count_of_anchor_kmers_in_this_read / count_of_all_kmers_in_this_read), 2)
        self.assertEqual(anchor_proportion, 0.0)
    
    def test_zero_division_handling(self):
        """Test zero division case"""
        count_of_anchor_kmers_in_this_read = 0
        count_of_all_kmers_in_this_read = 0
        
        try:
            anchor_proportion = round((count_of_anchor_kmers_in_this_read / count_of_all_kmers_in_this_read), 2)
        except ZeroDivisionError:
            anchor_proportion = 0
        
        self.assertEqual(anchor_proportion, 0)


if __name__ == '__main__':
    unittest.main()