#!/usr/bin/env bash
# PacBio HiFi/CCS full-length 16S: demultiplexed FASTQ -> relative abundance.
# Prerequisites: activated QIIME 2 amplicon environment; one FASTQ per sample;
# both 16S-specific primers are still present. NOT for ONT reads or subreads.
set -Eeuo pipefail
trap 'echo "ERROR: pipeline stopped at line $LINENO." >&2' ERR

# Paths are relative to the working directory unless overridden.
RAW_DIR="${RAW_DIR:-$PWD/raw_fastq}"
DB_DIR="${DB_DIR:-$PWD/db}"
REF_SEQS="${REF_SEQS:-$DB_DIR/silva-138-99-seqs.qza}"
REF_TAX="${REF_TAX:-$DB_DIR/silva-138-99-tax.qza}"
OUT="${OUT:-$PWD/qiime2_16s_results}"
THREADS="${THREADS:-16}"

# Initial QC settings; inspect retention statistics before accepting results.
MIN_LEN="${MIN_LEN:-1000}"
MAX_LEN="${MAX_LEN:-1600}"
MAX_EE="${MAX_EE:-2}"

# Biological primer sequences in their original 5' -> 3' orientation.
# Exclude the outer adapters/barcodes. q2-dada2 reverse-complements REV internally.
FWD='AGRGTTYGATYMTGGCTCAG'
REV='RGYTACCTTGTTACGACTT'

# Conservative example taxonomy settings, NOT species-identification cutoffs.
TAX_ID="${TAX_ID:-0.97}"
TAX_COV="${TAX_COV:-0.90}"
TAX_CONSENSUS="${TAX_CONSENSUS:-0.70}"

for cmd in qiime python; do
    command -v "$cmd" >/dev/null || { echo "Missing command: $cmd" >&2; exit 1; }
done
[[ -d "$RAW_DIR" ]] || { echo "FASTQ directory not found: $RAW_DIR" >&2; exit 1; }
for f in "$REF_SEQS" "$REF_TAX"; do
    [[ -s "$f" ]] || { echo "Reference artifact not found: $f" >&2; exit 1; }
done
[[ ! -e "$OUT" ]] || { echo "Output already exists; select a NEW OUT directory: $OUT" >&2; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer." >&2; exit 1; }
python -c 'import qiime2, biom, numpy, pandas'
qiime dada2 denoise-ccs --help >/dev/null

mkdir -p "$OUT"/{qc,abundance,tables,logs}
exec > >(tee "$OUT/logs/pipeline.log") 2>&1
qiime info > "$OUT/logs/qiime-info.txt"
printf 'FWD=%s\nREV=%s\nMIN_LEN=%s\nMAX_LEN=%s\nMAX_EE=%s\nTAX_ID=%s\nTAX_COV=%s\nTAX_CONSENSUS=%s\n' \
    "$FWD" "$REV" "$MIN_LEN" "$MAX_LEN" "$MAX_EE" "$TAX_ID" "$TAX_COV" "$TAX_CONSENSUS" \
    > "$OUT/logs/parameters.txt"
qiime tools peek "$REF_SEQS" > "$OUT/logs/reference-sequences.txt"
qiime tools peek "$REF_TAX" > "$OUT/logs/reference-taxonomy.txt"

# 1. Build a TAB-delimited V2 single-end manifest.
# Filename without .fastq[.gz] or .fq[.gz] becomes the sample ID.
python - "$RAW_DIR" "$OUT/manifest.tsv" <<'PY_MANIFEST'
import csv
import re
import sys
from pathlib import Path

raw, destination = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
suffixes = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
rows, seen = [], set()
for path in sorted(raw.iterdir()):
    if not path.is_file():
        continue
    suffix = next((s for s in suffixes if path.name.endswith(s)), None)
    if suffix is None:
        continue
    sample = path.name[:-len(suffix)]
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
        raise SystemExit(f'Use simple, nonempty sample filenames: {path.name}')
    if sample.lower() in {'id', 'sampleid', 'sample-id', 'sample id', 'featureid', 'feature-id', 'feature id'}:
        raise SystemExit(f'Reserved QIIME 2 sample ID: {sample}')
    if sample in seen:
        raise SystemExit(f'Duplicate sample ID: {sample}; combine files appropriately first.')
    if path.stat().st_size == 0:
        raise SystemExit(f'Empty FASTQ file: {path}')
    seen.add(sample)
    rows.append((sample, str(path.resolve())))
if not rows:
    raise SystemExit(f'No FASTQ files found directly inside {raw}')
with destination.open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['sample-id', 'absolute-filepath'])
    writer.writerows(rows)
print(f'Manifest: {len(rows)} samples; one file per sample.')
PY_MANIFEST

# 2. Import CCS reads as single-end sequences and inspect raw quality.
qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$OUT/manifest.tsv" \
    --input-format SingleEndFastqManifestPhred33V2 \
    --output-path "$OUT/demux.qza"

qiime demux summarize \
    --i-data "$OUT/demux.qza" \
    --output-dir "$OUT/qc/demux"

# 3. Remove primers + exterior tails, orient reads, filter, denoise, remove chimeras.
# Do NOT pre-trim these primers with cutadapt and then run this same command.
# No fixed-length truncation; the length filter acts on primer-removed reads.
# --output-dir captures all outputs, including extra diagnostics in newer QIIME 2.
qiime dada2 denoise-ccs \
    --i-demultiplexed-seqs "$OUT/demux.qza" \
    --p-front "$FWD" \
    --p-adapter "$REV" \
    --p-max-mismatch 2 \
    --p-no-indels \
    --p-trim-left 0 \
    --p-trunc-len 0 \
    --p-trunc-q 0 \
    --p-min-len "$MIN_LEN" \
    --p-max-len "$MAX_LEN" \
    --p-max-ee "$MAX_EE" \
    --p-pooling-method independent \
    --p-chimera-method consensus \
    --p-min-fold-parent-over-abundance 3.5 \
    --p-n-reads-learn 1000000 \
    --p-n-threads "$THREADS" \
    --output-dir "$OUT/dada2" \
    --verbose

qiime metadata tabulate \
    --m-input-file "$OUT/dada2/denoising_stats.qza" \
    --o-visualization "$OUT/qc/denoising-stats.qzv"

qiime feature-table summarize \
    --i-table "$OUT/dada2/table.qza" \
    --output-dir "$OUT/qc/table-before-taxonomy-filter"

# 4. Assign taxonomy against MATCHED full-length reference sequences/taxonomy.
qiime feature-classifier classify-consensus-vsearch \
    --i-query "$OUT/dada2/representative_sequences.qza" \
    --i-reference-reads "$REF_SEQS" \
    --i-reference-taxonomy "$REF_TAX" \
    --p-perc-identity "$TAX_ID" \
    --p-query-cov "$TAX_COV" \
    --p-min-consensus "$TAX_CONSENSUS" \
    --p-maxaccepts 10 \
    --p-strand both \
    --p-threads "$THREADS" \
    --output-dir "$OUT/taxonomy"

# 5. Exclude explicit non-target annotations, but RETAIN Unassigned and Archaea.
# Relative-abundance denominators include all retained features, not just named genera.
qiime taxa filter-table \
    --i-table "$OUT/dada2/table.qza" \
    --i-taxonomy "$OUT/taxonomy/classification.qza" \
    --p-exclude 'mitochondria,chloroplast,Eukaryota' \
    --p-mode contains \
    --o-filtered-table "$OUT/non-target-filtered.qza"

# Only remove zero-count samples. This is NOT a sequencing-depth QC cutoff.
qiime feature-table filter-samples \
    --i-table "$OUT/non-target-filtered.qza" \
    --p-min-frequency 1 \
    --o-filtered-table "$OUT/abundance/ASV-counts.qza"

qiime feature-table filter-seqs \
    --i-data "$OUT/dada2/representative_sequences.qza" \
    --i-table "$OUT/abundance/ASV-counts.qza" \
    --o-filtered-data "$OUT/representative-sequences-filtered.qza"

qiime feature-table summarize \
    --i-table "$OUT/abundance/ASV-counts.qza" \
    --output-dir "$OUT/qc/table-after-taxonomy-filter"

# 6. ASV relative abundance; no rarefaction is performed.
qiime feature-table relative-frequency \
    --i-table "$OUT/abundance/ASV-counts.qza" \
    --o-relative-frequency-table "$OUT/abundance/ASV-relative.qza"

# 7. Collapse COUNTS first, then calculate within-sample relative abundances.
# Assumes a standard seven-rank SILVA taxonomy. Default export stops at genus.
for entry in 2:phylum 3:class 4:order 5:family 6:genus; do
    level="${entry%%:*}"
    rank="${entry#*:}"
    qiime taxa collapse \
        --i-table "$OUT/abundance/ASV-counts.qza" \
        --i-taxonomy "$OUT/taxonomy/classification.qza" \
        --p-level "$level" \
        --o-collapsed-table "$OUT/abundance/${rank}-counts.qza"
    qiime feature-table relative-frequency \
        --i-table "$OUT/abundance/${rank}-counts.qza" \
        --o-relative-frequency-table "$OUT/abundance/${rank}-relative.qza"
done

# 8. Export plain TSV matrices and audit sample retention.
# QIIME 2 calculates all abundances; Python only exports and validates the results.
python - "$OUT" <<'PY_EXPORT'
import sys
from pathlib import Path
import biom
import numpy as np
import pandas as pd
from qiime2 import Artifact, Metadata

out = Path(sys.argv[1])
def read_table(path):
    return Artifact.load(str(path)).view(biom.Table).to_dataframe(dense=True)

all_counts = read_table(out / 'dada2/table.qza')
all_counts.to_csv(out / 'tables/ASV-counts-before-taxonomy-filter.tsv',
                  sep='\t', index_label='FeatureID')
for name in ('ASV', 'phylum', 'class', 'order', 'family', 'genus'):
    for kind in ('counts', 'relative'):
        frame = read_table(out / f'abundance/{name}-{kind}.qza')
        values = frame.to_numpy()
        if not np.isfinite(values).all() or (values < 0).any():
            raise SystemExit(f'Invalid abundance values: {name}-{kind}')
        if kind == 'relative' and not np.allclose(frame.sum(axis=0), 1.0,
                                                  rtol=0, atol=1e-8):
            raise SystemExit(f'Sample relative abundances do not sum to 1: {name}')
        frame.to_csv(out / f'tables/{name}-{kind}.tsv', sep='\t',
                     index_label='FeatureID' if name == 'ASV' else 'Taxon')

taxonomy = Artifact.load(str(out / 'taxonomy/classification.qza')).view(pd.DataFrame)
taxonomy.to_csv(out / 'tables/taxonomy.tsv', sep='\t', index_label='FeatureID')
stats = Artifact.load(str(out / 'dada2/denoising_stats.qza')).view(Metadata).to_dataframe()
stats.to_csv(out / 'qc/denoising-stats.tsv', sep='\t', index_label='sample-id')

ids = pd.read_csv(out / 'manifest.tsv', sep='\t', dtype=str)['sample-id']
kept = read_table(out / 'abundance/ASV-counts.qza')
audit = pd.DataFrame(index=pd.Index(ids, name='sample-id'))
audit['nonchimeric_reads'] = all_counts.sum(axis=0).reindex(audit.index, fill_value=0)
audit['retained_reads'] = kept.sum(axis=0).reindex(audit.index, fill_value=0)
audit['removed_by_taxonomy_filter'] = audit['nonchimeric_reads'] - audit['retained_reads']
audit['included_in_relative_tables'] = audit.index.isin(kept.columns)
audit.to_csv(out / 'qc/sample-retention.tsv', sep='\t')
print(f'Done. TSV rows = ASVs/taxa; columns = samples; relative values = 0 to 1.\n{out / "tables"}')
PY_EXPORT
