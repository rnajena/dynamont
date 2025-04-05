from subprocess import PIPE
import unittest
import numpy as np

class TestFileIO(unittest.TestCase):

    def test_feed_pipe(self):
        from src.python.segmentation.FileIO import feedPipe
        from unittest.mock import Mock
        mockPipe = Mock()
        mockPipe.communicate.return_value = ("result", "")
        signal = np.array([1.23456, 2.34567, 3.45678])
        read = "test_read"
        result = feedPipe(signal, read, mockPipe)
        expectedCookie = "1.23456,2.34567,3.45678\ntest_read\n"
        # Assert the communicate method was called with the expected input
        mockPipe.communicate.assert_called_once_with(input=expectedCookie)
        # Use self.assertEqual for consistency with unittest
        self.assertEqual(result, ("result", "", mockPipe.returncode))
        
    def test_format_segmentation_output(self):
        from src.python.segmentation.FileIO import formatSegmentationOutput
        output = "M2,0,1.000000;M3,2,1.000000;M4,871,0.895648;"
        sigOffset = 0
        lastIndex = 1000
        read = "ACGTACGTACGT"
        expected = np.array([
            [0, 2, 9, 'C', 'TACGT', 'M', 1.0, 'NA'],
            [2, 871, 8, 'A', 'GTACG', 'M', 1.0, 'NA'],
            [871, 1000, 7, 'T', 'CGTAC', 'M', 0.895648, 'NA']
        ], dtype=object)
        result = formatSegmentationOutput(output, sigOffset, lastIndex, read, 5)
        assert np.array_equal(result, expected)
        
    def test_format_segmentation(self):
        from src.python.segmentation.FileIO import formatSegmentation
        readid = "read123"
        signalid = "signal456"
        segmentation = np.array([[1, 2, 3], [4, 5, 6]])
        expected_output = "read123,signal456,1,2,3\nread123,signal456,4,5,6\n"
        result = formatSegmentation(readid, signalid, segmentation)
        assert result == expected_output
        
        # Correctly identifies and replaces outliers in a 1D numpy array
    def test_hampel_filter(self):
        from src.python.segmentation.FileIO import hampelFilter
        signal = np.array([1, 1, 1, 10, 1, 1, 1])
        expectedOutput = np.array([1, 1, 1, 1, 1, 1, 1])
    
        hampelFilter(signal)
    
        assert np.array_equal(signal, expectedOutput)
        
    # Count occurrences of A, C, G, T in a standard DNA sequence
    def test_count_nucleotides_in_standard_sequence(self):
        from src.python.segmentation.FileIO import countNucleotides
        # Given
        sequence = "ACGTACGT"
    
        # When
        result = countNucleotides(sequence)
    
        # Then
        self.assertEqual(result["A"], 2)
        self.assertEqual(result["C"], 2)
        self.assertEqual(result["G"], 2)
        self.assertEqual(result["T"], 2)

    # Calculate correct ratios for a sequence with equal distribution of nucleotides
    def test_equal_distribution_of_nucleotides(self):
        from src.python.segmentation.FileIO import countNucleotideRatios
        # Given
        sequence = "ACGTACGTACGTACGT"
        expected_ratios = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    
        # When
        actual_ratios = countNucleotideRatios(sequence)
    
        # Then
        self.assertEqual(expected_ratios, actual_ratios)

    # Writing a dictionary of kmer models with valid data to a file
    def test_write_kmer_models_with_valid_data(self):
        from src.python.segmentation.FileIO import writeKmerModels
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
        writeKmerModels(temp_filepath, kmer_models)

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
        from src.python.segmentation.FileIO import readKmerModels
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
        result = readKmerModels(temp_filepath)

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

    # Successfully opens a subprocess with a valid script path using unittest.mock
    def test_open_cpp_script_with_valid_path(self):
        from unittest.mock import patch, MagicMock
        # Given
        with patch('src.python.segmentation.FileIO.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_popen.return_value = mock_process
            script_path = ['path/to/valid/script.sh']
    
            # When
            from src.python.segmentation.FileIO import openCPPScript
            result = openCPPScript(script_path)
    
            # Then
            mock_popen.assert_called_once_with(
                script_path, 
                stdout=PIPE, 
                stdin=PIPE, 
                stderr=PIPE, 
                text=True
            )
            self.assertEqual(result, mock_process)

    # Creating a SegmentationError with a read parameter
    def test_segmentation_error_with_read_parameter(self):
        from src.python.segmentation.FileIO import SegmentationError
        # Given
        read = "read123"
    
        # When
        error = SegmentationError(read)
    
        # Then
        self.assertEqual(error.read, "read123")
        self.assertEqual(error.message, "No segmentation calculated for read123")
        self.assertEqual(str(error), "No segmentation calculated for read123")

    # Successfully calculates Z value when valid signal, read, params, and script are provided
    def test_calcz_returns_z_value_when_valid_inputs(self):
        from src.python.segmentation.FileIO import calcZ
        # Given
        from unittest.mock import Mock, patch
    
        signal = np.array([1.0, 2.0, 3.0])
        read = "ACGT"
        params = {"param1": 0.1, "param2": 0.2}
        script = "/path/to/script"
        model = "test_model"

        mock_pipe = Mock()
        with patch('src.python.segmentation.FileIO.openCPPScriptCalcZ', return_value=mock_pipe) as mock_openCPPScriptCalcZ, \
             patch('src.python.segmentation.FileIO.feedPipe', return_value=("3.14", "", 0)) as mock_feedPipe:
        
            # When
            result = calcZ(signal, read, params, script, model)
    
            # Then
            mock_openCPPScriptCalcZ.assert_called_once_with(script, params, model)
            mock_feedPipe.assert_called_once_with(signal, read, mock_pipe)
            self.assertEqual(result, 3.14)

    # Successfully feed signal and read to pipe and return expected results
    def test_feed_pipe_returns_expected_results(self):
        from src.python.segmentation.FileIO import feedPipe
        # Given
        import numpy as np
        from unittest.mock import MagicMock
        from subprocess import Popen

        signal = np.array([1.0, 2.0, 3.0])
        read = "test_read"
        mock_pipe = MagicMock(spec=Popen)
        mock_pipe.communicate.return_value = ("expected_result", "")
        mock_pipe.returncode = 0

        # When
        result, errors, returncode = feedPipe(signal, read, mock_pipe)

        # Then
        expected_input = "1.0,2.0,3.0\ntest_read\n"
        mock_pipe.communicate.assert_called_once()
        self.assertEqual(mock_pipe.communicate.call_args[1]['input'], expected_input)
        self.assertEqual(result, "expected_result")
        self.assertEqual(errors, "")
        self.assertEqual(returncode, 0)
        mock_pipe.kill.assert_called_once()

if __name__ == '__main__':
    unittest.main()