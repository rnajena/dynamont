from unittest.mock import MagicMock, patch, mock_open
from dynamont.segmentation.segment import get_raw

from collections import OrderedDict
from unittest.mock import MagicMock, patch


class TestSegment:

    def setup_method(self):
        import dynamont.segmentation.segment as segment

        segment.RAW_CACHE = OrderedDict()
        segment.RAW_CACHE_SIZE = 10

    def test_get_raw_opens_and_caches_path(self):
        import dynamont.segmentation.segment as segment

        raw = MagicMock()

        with patch(
            "python.segmentation.segment.open_pod5",
            return_value=raw,
        ) as mock_open:
            result = segment.get_raw("/tmp/read.pod5")

        assert result is raw
        mock_open.assert_called_once_with("/tmp/read.pod5")
        assert segment.RAW_CACHE["/tmp/read.pod5"] is raw

    def test_get_raw_returns_cached_value(self):
        import dynamont.segmentation.segment as segment

        raw = MagicMock()
        segment.RAW_CACHE["/tmp/read.pod5"] = raw

        with patch(
            "python.segmentation.segment.open_pod5"
        ) as mock_open:
            result = segment.get_raw("/tmp/read.pod5")

        assert result is raw
        mock_open.assert_not_called()

    def test_get_raw_moves_cached_path_to_end(self):
        import dynamont.segmentation.segment as segment

        raw1 = MagicMock()
        raw2 = MagicMock()

        segment.RAW_CACHE["/tmp/read1.pod5"] = raw1
        segment.RAW_CACHE["/tmp/read2.pod5"] = raw2

        result = segment.get_raw("/tmp/read1.pod5")

        assert result is raw1
        assert list(segment.RAW_CACHE) == [
            "/tmp/read2.pod5",
            "/tmp/read1.pod5",
        ]

    def test_get_raw_evicts_oldest_entry_when_cache_is_full(self, monkeypatch):
        import dynamont.segmentation.segment as segment

        monkeypatch.setattr(segment, "RAW_CACHE_SIZE", 2)

        old_raw = MagicMock()
        existing_raw = MagicMock()
        new_raw = MagicMock()

        segment.RAW_CACHE["/tmp/read1.pod5"] = old_raw
        segment.RAW_CACHE["/tmp/read2.pod5"] = existing_raw

        with patch(
            "python.segmentation.segment.open_pod5",
            return_value=new_raw,
        ) as mock_open:
            result = segment.get_raw("/tmp/read3.pod5")

        assert result is new_raw
        mock_open.assert_called_once_with("/tmp/read3.pod5")

        assert list(segment.RAW_CACHE) == [
            "/tmp/read2.pod5",
            "/tmp/read3.pod5",
        ]

        old_raw.close.assert_called_once()

    def test_get_raw_ignores_error_when_closing_evicted_entry(
        self, monkeypatch
    ):
        import dynamont.segmentation.segment as segment

        monkeypatch.setattr(segment, "RAW_CACHE_SIZE", 1)

        old_raw = MagicMock()
        old_raw.close.side_effect = RuntimeError("close failed")

        new_raw = MagicMock()

        segment.RAW_CACHE["/tmp/read1.pod5"] = old_raw

        with patch(
            "python.segmentation.segment.open_pod5",
            return_value=new_raw,
        ):
            result = segment.get_raw("/tmp/read2.pod5")

        assert result is new_raw
        assert list(segment.RAW_CACHE) == ["/tmp/read2.pod5"]
        old_raw.close.assert_called_once()

    def test_listener_writes_header_results_and_errors(self):
        from dynamont.segmentation.segment import listener
        outfile = "/tmp/results.csv.zst"

        queue = MagicMock()
        queue.get.side_effect = [
            b"read1,signal1,0,10,0,A,motif,state,0.99,polished\n",
            "error: failed to segment read2",
            b"read3,signal3,20,30,20,G,motif,state,0.95,polished\n",
            "kill",
        ]

        raw_file = MagicMock()
        output = MagicMock()
        stream_writer = MagicMock()
        stream_writer.__enter__.return_value = output

        compressor = MagicMock()
        compressor.stream_writer.return_value = stream_writer

        pbar = MagicMock()
        tqdm_cm = MagicMock()
        tqdm_cm.__enter__.return_value = pbar

        with (
            patch("python.segmentation.segment.zstd.ZstdCompressor", return_value=compressor) as mock_compressor,
            patch("python.segmentation.segment.open", mock_open()) as mock_open_file,
            patch("python.segmentation.segment.tqdm", return_value=tqdm_cm),
        ):
            # Make the outfile open return our mocked raw file.
            mock_open_file.return_value.__enter__.return_value = raw_file

            listener(queue, outfile)

        mock_compressor.assert_called_once_with(level=3)
        compressor.stream_writer.assert_called_once_with(raw_file)

        # Header + two successful results should be written.
        assert output.write.call_count == 3

        output.write.assert_any_call(
            b"readid,signalid,start,end,basepos,base,motif,state,"
            b"posterior_probability,polish\n"
        )
        output.write.assert_any_call(
            b"read1,signal1,0,10,0,A,motif,state,0.99,polished\n"
        )
        output.write.assert_any_call(
            b"read3,signal3,20,30,20,G,motif,state,0.95,polished\n"
        )

        # The error should go to the .errors file.
        mock_open_file.assert_any_call("/tmp/results.errors", "a")

        error_handle = mock_open_file.return_value.__enter__.return_value
        error_handle.write.assert_called_once_with(
            "error: failed to segment read2\n"
        )

        # Three actual results were processed.
        assert pbar.update.call_count == 3
        pbar.update.assert_any_call(1)

        pbar.set_postfix.assert_called_once_with(errors=1)

    def test_hampel_three_point_fast_path_matches_reference(self):
        import numpy as np
        from dynamont.segmentation.utils import hampel

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
        from dynamont.segmentation.segment import parse
    
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
        from collections import OrderedDict
        from dynamont.segmentation import segment

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