class TestSegment():

    # Returns memory usage in MB for the current process
    def test_returns_memory_usage_in_mb(self, mocker):
        from src.python.segmentation.segment import get_memory_usage
        import psutil
        # Given
        mock_process = mocker.Mock()
        mock_memory_info = mocker.Mock()
        mock_memory_info.rss = 104857600  # 100 MB in bytes
        mock_process.memory_info.return_value = mock_memory_info
        mocker.patch('psutil.Process', return_value=mock_process)
        mocker.patch('src.python.segmentation.segment.getpid', return_value=12345)
    
        # When
        result = get_memory_usage()
    
        # Then
        assert result == 100.0  # 100 MB
        psutil.Process.assert_called_once_with(12345)
        mock_process.memory_info.assert_called_once()

    # Parse command line arguments with all required parameters
    def test_parse_with_all_required_parameters(self):
        # Given
        import sys
        from argparse import Namespace
        from src.python.segmentation.segment import parse
    
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
            '-p', 'rna_r9'
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
            assert args.pore == 'rna_r9'
            assert args.processes > 0  # Default value should be set
            assert args.qscore == 0.0  # Default value
        finally:
            # Restore original sys.argv
            sys.argv = original_argv