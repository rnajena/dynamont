source ~/.bashrc
conda activate science

python ~/projects/dynamont/src/f5c/estimateModelsFromF5C.py --raw /data/fass5/reads/rna_modifications/raw/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/20231019_1641_MN21435_FAX28269_36c48ee6/pod5/RNAmod_20231019_EcoliWT_R1_FAX28269_WG.pod5 --fastx /data/fass5/reads/rna_modifications/basecalls/RNAmod_20231019_EcoliWT_R1_FAX28269_WG_dorado_server/qualityfiltered/RNAmod_20231019_EcoliWT_R1_FAX28269_WG_q25.fastq --out /home/yi98suv/projects/dynamont/data/ecoliWTModels/ --polya /data/fass5/jannes/rna_modifications/kai_ecoli/segmentation/k12_wt/polya.csv --segmentationPickle /data/fass5/reads/rna_modifications/segmentation/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/resquiggle/segmentation.pickle

python ../extend_model.py -a ACGT -k 5 pAmodels.model.filled pAmodels.model.extended

# added a quality filter
python ~/projects/dynamont/src/f5c/estimateModelsFromF5C.py --raw /data/fass5/reads/rna_modifications/raw/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/20231019_1641_MN21435_FAX28269_36c48ee6/pod5/RNAmod_20231019_EcoliWT_R1_FAX28269_WG.pod5 --fastx /data/fass5/reads/rna_modifications/basecalls/RNAmod_20231019_EcoliWT_R1_FAX28269_WG_dorado_server/qualityfiltered/RNAmod_20231019_EcoliWT_R1_FAX28269_WG_q25.fastq --out /home/yi98suv/projects/dynamont/data/ecoliWTModels_q15/ --polya /data/fass5/jannes/rna_modifications/kai_ecoli/segmentation/k12_wt/polya.csv --segmentationPickle /data/fass5/reads/rna_modifications/segmentation/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/resquiggle/segmentation.pickle -q 15
