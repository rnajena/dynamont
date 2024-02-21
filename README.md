# ont_segmentation

## training

<!-- ```bash
python src/train.py --raw /data/fass5/reads/rna_modifications/raw/A-RNA_20201014_FAO07649/20201014_FAO07649_A_RNA/20201014_1657_MN34689_FAO07649_ca95a109/fast5/ --fastq /data/fass5/reads/rna_modifications/basecalls/A-RNA_20201014_FAO07649_dorado_server/A-RNA_20201014_FAO07649.fastq --out /data/fass5/jannes/ont_segmentation/training_A-RNA_20201014/
``` -->

- first training
```bash
python train.py --raw /data/fass5/reads/rna_modifications/raw/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/20231019_1641_MN21435_FAX28269_36c48ee6/pod5/RNAmod_20231019_EcoliWT_R1_FAX28269_WG.pod5 --fastq /data/fass5/reads/rna_modifications/basecalls/RNAmod_20231019_EcoliWT_R1_FAX28269_WG_dorado_server/RNAmod_20231019_EcoliWT_R1_FAX28269_WG.fastq --out /data/fass5/jannes/ont_segmentation/training_ecoli_wt/ --batch_size 5 --polya /data/fass5/jannes/rna_modifications/kai_ecoli/segmentation/k12_wt/polya.csv
```

- second training
- removed prefix and suffix
- switched to multiprocessing in training (#threads = batch_size)
```bash
python train_mp.py --raw /data/fass5/reads/rna_modifications/raw/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/RNAmod_20231019_EcoliWT_R1_FAX28269_WG/20231019_1641_MN21435_FAX28269_36c48ee6/pod5/RNAmod_20231019_EcoliWT_R1_FAX28269_WG.pod5 --fastq /data/fass5/reads/rna_modifications/basecalls/RNAmod_20231019_EcoliWT_R1_FAX28269_WG_dorado_server/RNAmod_20231019_EcoliWT_R1_FAX28269_WG.fastq --out /data/fass5/jannes/ont_segmentation/training_ecoli_wt_basic/ --batch_size 16 --polya /data/fass5/jannes/rna_modifications/kai_ecoli/segmentation/k12_wt/polya.csv
```
