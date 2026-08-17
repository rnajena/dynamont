from collections import OrderedDict

class TestSegment():

    def test_hampel_three_point_fast_path_matches_reference(self):
        import numpy as np
        from python.segmentation.utils import hampel

        def reference(signal):
            original = signal.copy()
            window = original[:3].copy()
            for i in range(1, len(signal) - 2):
                median = np.median(window)
                mad = 1.4826 * np.median(np.abs(window - median))
                if np.abs(original[i] - median) > 3.0 * mad:
                    signal[i] = median
                window[:-1] = window[1:]
                window[-1] = original[i + 2]

        values = np.array([1.0, 1.0, 50.0, 1.0, 1.0, 1.0, 75.0, 1.0])
        expected = values.copy()
        actual = values.copy()
        reference(expected)
        hampel(actual)

        np.testing.assert_array_equal(actual, expected)

    # Parse command line arguments with all required parameters
    def test_parse_with_all_required_parameters(self):
        # Given
        import sys
        from argparse import Namespace
        from python.segmentation.segment import parse
    
        # Save original sys.argv
        original_argv = sys.argv.copy()
    
        # When
        sys.argv = [
            'dynamont-resquiggle',
            '-r', '/path/to/raw.pod5',
            '-b', '/path/to/basecalls.bam',
            '-o', '/path/to/output',
            '--mode', 'basic',
            '--model_path', '/path/to/model',
            '-p', 'rna002'
        ]
    
        try:
            # Then
            args = parse()
            assert isinstance(args, Namespace)
            assert args.raw == '/path/to/raw.pod5'
            assert args.basecalls == '/path/to/basecalls.bam'
            assert args.outfile == '/path/to/output'
            assert args.mode == 'basic'
            assert args.model_path == '/path/to/model'
            assert args.pore == 'rna002'
            assert args.processes > 0  # Default value should be set
            assert args.qscore == 0.0  # Default value
        finally:
            # Restore original sys.argv
            sys.argv = original_argv

    def test_close_raw_cache_closes_all_cached_readers(self):
        from python.segmentation import segment

        class FakeReader:
            def __init__(self, name):
                self.name = name
                self.closed = False

            def close(self):
                self.closed = True

        reader1 = FakeReader('one')
        reader2 = FakeReader('two')
        segment.RAW_CACHE = OrderedDict([
            ('a', reader1),
            ('b', reader2),
        ])

        segment.close_raw_cache()

        assert segment.RAW_CACHE == OrderedDict()
        assert reader1.closed is True
        assert reader2.closed is True