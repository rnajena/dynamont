class TestTrain:
    # Parse command line arguments with all required arguments provided
    def test_parse_with_all_required_arguments(self):

        # Given
        import sys
        from argparse import Namespace
        from unittest.mock import patch
        from src.python.segmentation.train import parse
    
        test_args = [
            'dynamont-train',
            '-r', '/path/to/raw',
            '-b', '/path/to/basecalls.bam',
            '-o', '/path/to/outdir',
            '-p', 'rna_r9',
            '--mode', 'basic',
            '--model_path', '/path/to/model'
        ]
    
        # When
        with patch.object(sys, 'argv', test_args):
            args = parse()
    
        # Then
        assert isinstance(args, Namespace)
        assert args.raw == '/path/to/raw'
        assert args.basecalls == '/path/to/basecalls.bam'
        assert args.outdir == '/path/to/outdir'
        assert args.pore == 'rna_r9'
        assert args.mode == 'basic'
        assert args.model_path == '/path/to/model'
        assert args.batch_size == 24  # Default value
        assert args.epochs == 1  # Default value
        assert args.qscore == 10  # Default value

  