# Upstream processing of full-length 16S rRNA and PhIP-seq data

This repository contains upstream sequencing-data processing scripts for the study *Epitope-resolved precision serology for pediatric Helicobacter pylori infection*.

## Full-length 16S rRNA amplicon processing

`run_qiime2_full_length_16s.sh` processes demultiplexed PacBio HiFi/CCS full-length 16S rRNA amplicon reads using QIIME 2. It performs primer removal, quality filtering, denoising, chimera removal, and taxonomic assignment. Outputs include ASV counts, taxonomic annotations, relative-abundance tables at the ASV and phylum-to-genus levels, and quality-control summaries for downstream microbiome analysis.

## PhIP-seq processing

`PhIP-seq_15wHP.py` processes paired-end sequencing data from the HpScan PhIP-seq assay. It performs read quality control and merging, flanking-sequence removal, translation into peptide sequences, barcode-based sample assignment, and peptide counting. Sample counts are combined with input-library counts to produce a peptide-by-sample count matrix (`Count.txt.gz`) for downstream antibody-reactivity analysis.

Supporting files:

- `barcode_new.txt`: barcode-to-sample mapping used to assign sequencing reads to individual samples.
- `HP_15W.txt`: reference peptide table containing peptide sequences and input-library counts.

