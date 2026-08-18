import unittest
import numpy as np

class TestUtils(unittest.TestCase):
        
        # Correctly identifies and replaces outliers in a 1D numpy array
    def test_hampel_filter(self):
        from dynamont.segmentation.utils import hampel
        signal = np.array([1, 1, 1, 10, 1, 1, 1])
        expectedOutput = np.array([1, 1, 1, 1, 1, 1, 1])
    
        hampel(signal)
    
        assert np.array_equal(signal, expectedOutput)
        
    # Count occurrences of A, C, G, T in a standard DNA sequence
    def test_count_nucleotides_in_standard_sequence(self):
        from dynamont.segmentation.utils import cnt_nts
        # Given
        sequence = "ACGTACGT"
    
        # When
        result = cnt_nts(sequence)
    
        # Then
        self.assertEqual(result["A"], 2)
        self.assertEqual(result["C"], 2)
        self.assertEqual(result["G"], 2)
        self.assertEqual(result["T"], 2)

    # Calculate correct ratios for a sequence with equal distribution of nucleotides
    def test_equal_distribution_of_nucleotides(self):
        from dynamont.segmentation.utils import cnt_nts_ratios
        # Given
        sequence = "ACGTACGTACGTACGT"
        expected_ratios = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    
        # When
        actual_ratios = cnt_nts_ratios(sequence)
    
        # Then
        self.assertEqual(expected_ratios, actual_ratios)

    # Writing a dictionary of kmer models with valid data to a file
    def test_write_kmer_models_with_valid_data(self):
        from dynamont.segmentation.utils import write_kmer_model
        # Given
        import os
        import tempfile
    
        kmer_models = {
            'ACGT': (1.5, 0.2),
            'TGCA': (2.1, 0.3),
            'GATC': (1.8, 0.25)
        }

        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_filepath = temp_file.name

        # When
        write_kmer_model(temp_filepath, kmer_models)

        # Then
        with open(temp_filepath, 'r') as f:
            content = f.readlines()

        self.assertEqual(content[0], 'kmer\tlevel_mean\tlevel_stdv\n')
        self.assertEqual(len(content), 4)  # Header + 3 kmers

        # Check each kmer line is correctly formatted
        for kmer in kmer_models:
            mean, stdev = kmer_models[kmer]
            expected_line = f'{kmer}\t{mean}\t{stdev}\n'
            self.assertIn(expected_line, content)

        # Clean up
        os.unlink(temp_filepath)

    # Read a valid TSV file with kmer, level_mean, and level_stdv columns
    def test_read_valid_kmer_model_file(self):
        from dynamont.segmentation.utils import read_kmer_model
        # Given
        import tempfile
        import os
        import pandas as pd

        test_data = pd.DataFrame({
            'kmer': ['AAAAA', 'AAAAC', 'AAAAG'],
            'level_mean': [100.0, 110.0, 120.0],
            'level_stdv': [10.0, 11.0, 12.0]
        })

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as temp_file:
            test_data.to_csv(temp_file.name, sep='\t', index=False)
            temp_filepath = temp_file.name

        # When
        result = read_kmer_model(temp_filepath)

        # Then
        expected = {
            'AAAAA': (100.0, 10.0),
            'AAAAC': (110.0, 11.0),
            'AAAAG': (120.0, 12.0)
        }

        self.assertEqual(len(result), 3)
        for kmer, (mean, stdv) in expected.items():
            self.assertIn(kmer, result)
            self.assertEqual(result[kmer][0], mean)
            self.assertEqual(result[kmer][1], stdv)

        # Cleanup
        os.unlink(temp_filepath)

    # Creating a SegmentationError with a read parameter
    def test_segmentation_error_with_read_parameter(self):
        from dynamont.segmentation.utils import SegmentationError
        # Given
        read = "read123"
    
        # When
        error = SegmentationError(read)
    
        # Then
        self.assertEqual(error.read, "read123")
        self.assertEqual(error.message, "No segmentation calculated for read123")
        self.assertEqual(str(error), "No segmentation calculated for read123")

if __name__ == '__main__':
    unittest.main()