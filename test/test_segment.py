class TestSegment():

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