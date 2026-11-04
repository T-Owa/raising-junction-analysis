# RAISING downstream analysis — public source release

Version **1.0.0**. Paired-read classification and donor–host junction analysis for the bundled **800-bp donor and its fixed 3′/5′ final primers**. This is a study-specific downstream workflow, not a general structural-variant caller. It is separate from the sequencing provider's mapping/annotation software.

[日本語ガイド](README_ja.md) · [Methods and definitions](docs/METHODS.md) · [Validation](docs/VALIDATION.md) · [Reference provenance](docs/REFERENCE_PROVENANCE.md)

## Included analyses

1. Classify **every input read pair** into four mutually exclusive categories.
2. Discover conservative donor–host junction templates using exact independent anchors.
3. Reassess reads against existing nearby templates, allowing bounded substitutions and small indels.
4. Export library-level microhomology (MH), donor-terminal loss and intervening-sequence distributions, with denominators and support-threshold sensitivity analyses.

No original sequencing reads, subject/sample records, vendor result tables, personal paths or original execution logs are included. The donor/primer references and generated demonstration reads are not human-subject sequencing data. See [sharing and provenance](docs/SHARING.md) before publishing newly generated results.

## Install

Use Linux or Ubuntu under WSL. The scripts use POSIX file locking and do not run directly in native Windows. macOS may work when the dependencies are available, but was not validated for this release.

```bash
micromamba create -y -f environment.yml --strict-channel-priority
micromamba activate raising-public
python tools/check_environment.py
python tests/run_all.py
```

`environment.yml` records known study Python/Biopython/BWA/BLAST version requirements. It is an installation specification, **not a complete platform-specific lockfile**. Record the solved environment for each analysis. The runtime requires Python, Biopython, BWA, BLAST+ (`blastn` and `makeblastdb`) and samtools. The workflow does not call fastp, pandas or pysam.

## Run the synthetic demonstration first

The host sequence in this demonstration is generated from a fixed random seed; it is not GRCh38. The small example checks execution and intended categories, not accuracy on biological material.

```bash
python examples/make_demo.py --outdir ./demo-work
python pipelines/fourway/raising_batch.py \
  --root ./demo-work --host-fasta ./demo-work/ref/synthetic_host.fa \
  --threads 2 --allow-index-build
python pipelines/junction/raising_junction_batch.py \
  --root ./demo-work --expected-count 2
```

The junction runner prints the completed batch directory, under `demo-work/results/junction_batch_v1/`. Use that **specific completed run** for the next command:

```bash
python tools/prepare_figure_data.py \
  --input-dir ./demo-work/results/junction_batch_v1/REPLACE_WITH_COMPLETED_RUN \
  --output-dir ./demo-work/figure-data \
  --metadata ./demo-work/sample_metadata.csv
```

The demonstration's expected values are described in `demo-work/expected_demo.json`. The generator refuses a nonempty destination.

## Analyze study libraries

Use a **new analysis root**, separate from previous private runs. Sanitizing source text changes source hashes; the original provenance checks intentionally prevent reuse of incompatible completed chunks. Do not edit old manifests to bypass these checks.

```bash
mkdir -p ./work/vendor_results/Raw_fastq ./work/ref
```

Place or link paired raw files into `work/vendor_results/Raw_fastq/`. The side is inferred from `_3prm_` or `_5prm_` in each sample name. For example:

```text
sample001_3prm_S1_L001_R1_001.fastq.gz
sample001_3prm_S1_L001_R2_001.fastq.gz
sample001_5prm_S2_L001_R1_001.fastq.gz
sample001_5prm_S2_L001_R2_001.fastq.gz
```

R1/R2 record identifiers and record counts must match. Multiple lanes/parts are supported. Duplicate input identifiers and missing mates fail validation. The batch runner supports `.fastq` and `.fastq.gz`. Retain the raw files for junction reassessment.

Place the chosen genomic FASTA in `work/ref/host.fa.gz` (uncompressed FASTA is also supported). Use the **same exact reference file and assembly across compared libraries**; see the reference provenance notes before reproducing a published run.

```bash
# Inspect all detected library/mate combinations.
python pipelines/fourway/raising_batch.py \
  --root ./work --host-fasta ./work/ref/host.fa.gz --plan-only

# Process all pairs. Index creation is explicitly authorized on the first run.
python pipelines/fourway/raising_batch.py \
  --root ./work --host-fasta ./work/ref/host.fa.gz \
  --threads 4 --allow-index-build

# Inspect the libraries available for junction analysis.
python pipelines/junction/raising_junction_batch.py --root ./work --dry-run

# Replace 2 with the actual intended number of libraries (not biological samples).
python pipelines/junction/raising_junction_batch.py --root ./work --expected-count 2
```

The four-way batch automatically resumes compatible completed work when the same command is rerun. Junction analysis also reuses only compatible completed results. The 10,000-pair chunk size limits working memory; it is **not random subsampling**. Initial human-reference indexing is considerably larger than the synthetic demonstration. Indexing can be largely single-threaded; `--threads` is not a guarantee of parallel indexing.

Use `--samples sample001_3prm_S1 sample001_5prm_S2` to select libraries. When changing the selected four-way batch scope, use a separate `--outdir` for its aggregate tables, as documented by the CLI. Do not combine 3′ and 5′ libraries as independent biological replicates.

## Figure-ready tables

The complete four-way aggregate is:

```text
work/results/all_libraries_fourway_batch_v1/all_figure_counts.csv
```

Its `count`, `total_pairs`, `fraction` and `percent` describe **input read-pair composition**. `partial_*` tables indicate an incomplete batch and must not be treated as the final complete result.

For junction figures, use the specific completed batch directory printed by the junction runner:

```bash
python tools/prepare_figure_data.py \
  --input-dir ./work/results/junction_batch_v1/REPLACE_WITH_COMPLETED_RUN \
  --output-dir ./work/figure-data \
  --metadata ./examples/metadata.csv
```

Copy `examples/metadata.csv` and replace its sample identifiers with your actual identifiers before using it. The metadata argument is optional. No condition or biological replicate is inferred from a sample name. Conflicting metadata fail validation.

The exporter accepts native `.tsv.gz` output or decompressed `.tsv`; provide only one format for each table. It requires the complete batch summary and all three matching merged tables. Key outputs:

| File | Use |
|---|---|
| `figure_values.csv` | Library-level fractions and donor-loss medians, with explicit denominators |
| `metrics.csv` | Integer counts underlying those fractions |
| `mh_distribution.csv` | Direct-junction MH-length distribution |
| `deletion_distribution.csv` | Minimum/maximum donor-terminal-loss endpoint distributions |
| `insertion_distribution.csv` | Lengths of intervening-sequence candidates |
| `qc.csv` | Recovery, unresolved support and singleton checks |
| `alternative_review.csv` | Alternative substitution/small-indel representations requiring interpretation |
| `candidate_annotations.csv` | Original candidate annotations and reassessed support |
| `analysis_summary.json` | Source hashes, settings and definitions |

For the descriptive primary comparison, select `subset=nuclear_all`, `min_support=1`, `support_basis=all_supported_pairs`, retaining 3′ and 5′ libraries separately. Use `min_support=2` as a sensitivity analysis; two read pairs are not necessarily two independent molecules. The `nuclear_all` subset includes nuclear/unplaced reference contigs and excludes mitochondrial calls; it is not restricted to canonical chromosomes. The other subsets are defined in `docs/METHODS.md`.

**Interpretation:** recovered sequence patterns are not validated independent integration events; MH does not identify a repair pathway. Insert-containing junctions have MH **unassessed/NA**, not zero. Donor-terminal loss is not host-genome deletion. Reassessment supports existing sequence templates; it does not automatically resolve every physical indel or discover all complex junctions. The fixed profile's nominal direct-junction donor-loss limits are 39 bp (3′) and 45 bp (5′), and its intervening-sequence search limit is 20 bp.

## Release and citation

The public package version is separate from internal analysis-engine versions: four-way 3′ `1.0.0`, four-way 5′ `1.1.0`, exact-anchor pilot `0.1.0`, junction refinement `0.1.1`. See `docs/CHANGELOG.md` for packaging changes and validation scope.

The package includes a proposed MIT license. Public author/repository/DOI metadata are deliberately left to the authors in `docs/CITATION.cff.example` and `docs/CODE_AVAILABILITY.md`. No repository or DOI has been registered by generating this archive. Archive the actual paper release and cite its version-specific persistent identifier.
