source ~/.bashrc
conda activate science

python ~/projects/dynamont/src/f5c/estimateModelsFromF5C.py --raw /data/fass5/reads/rna_modifications/raw/A-RNA_20201014_FAO07649/A-RNA_20201014_FAO07649.pod5 --fastx /data/fass5/reads/rna_modifications/basecalls/A-RNA_20201014_FAO07649_dorado_server/A-RNA_20201014_FAO07649.fastq --out /home/yi98suv/projects/dynamont/data/bastiIVTA/ --polya /data/fass5/jannes/rna_modifications/IVT/polyA/A-RNA_20201014_FAO07649.csv --segmentationPickle /data/fass5/reads/rna_modifications/segmentation/A-RNA_20201014_FAO07649/resquiggle/segmentation.pickle