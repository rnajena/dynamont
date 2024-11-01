import unittest

class TestSum(unittest.TestCase):

    def test_valid_ndarray_and_string_input(self):
        import numpy as np
        from src.dynamont.FileIO import feedPipe
        from unittest.mock import Mock
        mock_pipe = Mock()
        mock_pipe.communicate.return_value = ("result", None)
        signal = np.array([1.23456, 2.34567, 3.45678])
        read = "test_read"
        result = feedPipe(signal, read, mock_pipe)
        expected_cookie = "1.23456,2.34567,3.45678\ntest_read\n"
        # Assert the communicate method was called with the expected input
        mock_pipe.communicate.assert_called_once_with(input=expected_cookie)
        # Use self.assertEqual for consistency with unittest
        self.assertEqual(result, "result")
        
    def test_parses_well_formed_output(self):
        import numpy as np
        from src.dynamont.FileIO import formatSegmentationOutput
        output = "M2,0,1.000000;M3,2,1.000000;M4,871,0.895648;"
        signal_offset = 0
        lastIndex = 1000
        read = "ACGTACGTACGT"
        expected = np.array([
            [0, 2, 9, 'C', 'TACGT', 'M', 1.0, 'NA'],
            [2, 871, 8, 'A', 'GTACG', 'M', 1.0, 'NA'],
            [871, 1000, 7, 'T', 'CGTAC', 'M', 0.895648, 'NA']
        ], dtype=object)
        result = formatSegmentationOutput(output, signal_offset, lastIndex, read)
        assert np.array_equal(result, expected)
        
    def test_format_segmentation_correctly_formats_array(self):
        import numpy as np
        from src.dynamont.FileIO import formatSegmentation
        readid = "read123"
        signalid = "signal456"
        segmentation = np.array([[1, 2, 3], [4, 5, 6]])
        expected_output = "read123,signal456,1,2,3\nread123,signal456,4,5,6\n"
        result = formatSegmentation(readid, signalid, segmentation)
        assert result == expected_output
        
        # Correctly identifies and replaces outliers in a 1D numpy array
    def test_identifies_and_replaces_outliers(self):
        import numpy as np
        from src.dynamont.FileIO import hampel_filter
    
        signal = np.array([1, 1, 1, 10, 1, 1, 1])
        expected_output = np.array([1, 1, 1, 1, 1, 1, 1])
    
        filtered_signal = hampel_filter(signal)
    
        assert np.array_equal(filtered_signal, expected_output)
        
if __name__ == '__main__':
    unittest.main()