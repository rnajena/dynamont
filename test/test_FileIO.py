import unittest

class TestSum(unittest.TestCase):

    def testFeedPipe(self):
        import numpy as np
        from src.dynamont.FileIO import feedPipe
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
        
    def testFormatSegmentationOutput(self):
        import numpy as np
        from src.dynamont.FileIO import formatSegmentationOutput
        output = "M2,0,1.000000;M3,2,1.000000;M4,871,0.895648;"
        sigOffset = 0
        lastIndex = 1000
        read = "ACGTACGTACGT"
        expected = np.array([
            [0, 2, 9, 'C', 'TACGT', 'M', 1.0, 'NA'],
            [2, 871, 8, 'A', 'GTACG', 'M', 1.0, 'NA'],
            [871, 1000, 7, 'T', 'CGTAC', 'M', 0.895648, 'NA']
        ], dtype=object)
        result = formatSegmentationOutput(output, sigOffset, lastIndex, read)
        assert np.array_equal(result, expected)
        
    def testFormtSegmentation(self):
        import numpy as np
        from src.dynamont.FileIO import formatSegmentation
        readid = "read123"
        signalid = "signal456"
        segmentation = np.array([[1, 2, 3], [4, 5, 6]])
        expected_output = "read123,signal456,1,2,3\nread123,signal456,4,5,6\n"
        result = formatSegmentation(readid, signalid, segmentation)
        assert result == expected_output
        
        # Correctly identifies and replaces outliers in a 1D numpy array
    def testHampelFilter(self):
        import numpy as np
        from src.dynamont.FileIO import hampelFilter
    
        signal = np.array([1, 1, 1, 10, 1, 1, 1])
        expectedOutput = np.array([1, 1, 1, 1, 1, 1, 1])
    
        filteredSignal = hampelFilter(signal)
    
        assert np.array_equal(filteredSignal, expectedOutput)
        
if __name__ == '__main__':
    unittest.main()