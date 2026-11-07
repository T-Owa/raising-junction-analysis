# Reviewer guide: RAISING downstream analysis

Software version: **1.0.0**. Repository: [T-Owa/raising-junction-analysis](https://github.com/T-Owa/raising-junction-analysis).

This guide maps the analyses to their executable entry points and outputs. Commands are run from the repository root, whether obtained from GitHub or extracted from the review archive.

## Scope of the supplied material

The archive contains the analysis source, the fixed **800-bp donor** and final primers, synthetic-data generation, and tests. It supports inspection of the classification and junction rules and execution of a complete synthetic FASTQ-to-summary example.

Study FASTQs, completed study junction tables and the study host-genome FASTA are **not included**. The synthetic example tests specified behaviors; it is not a replay of the biological experiment or a sensitivity/specificity benchmark. The exact study host-FASTA release and checksum are not recorded here; see [Reference provenance](REFERENCE_PROVENANCE.md).

## Analysis-to-output map

| Analysis | Entry point | Principal outputs and unit of analysis |
|---|---|---|
| Four-way classification | `pipelines/fourway/raising_batch.py` | `all_figure_counts.csv`; counts and fractions of all input read pairs, separately by library |
| Donor–host junction discovery and reassessment | `pipelines/junction/raising_junction_batch.py` | `batch_summary.json` and merged candidate, reassessment and template-comparison tables |
| Junction summaries | `tools/prepare_figure_data.py` | `figure_values.csv`, `metrics.csv`, MH/donor-loss/intervening-sequence distributions; supported candidate sequence groups within each library |
| Quality control and alternative representations | Same exporter | `qc.csv`, `alternative_review.csv`, `candidate_annotations.csv` and `analysis_summary.json` |

The junction batch runner coordinates exact-anchor discovery and reassessment against existing templates. It is the entry point for routine analysis; the versioned engine directories retain the underlying implementations.

## Recommended verification

Use **Linux or Ubuntu under WSL**. Native Windows lacks the required POSIX file locking. Install the environment as described in the [README](../README.md). Then run:

```bash
python tools/check_environment.py
python tests/run_all.py --include-external
```

The second command runs the six unit/synthetic suites and the actual-tool round trip using BWA, BLAST+, samtools and Biopython. Successful completion exits with status 0 and prints an external-round-trip success message. Some tests deliberately produce validation errors as part of their expected failure cases.

The round trip uses temporary output directories. To retain CSV files for inspection, follow the persistent [synthetic demonstration](../README.md#run-the-synthetic-demonstration-first) in the README.

Expected results for each synthetic library (one 3′ library and one 5′ library):

| Quantity | Expected value per library |
|---|---:|
| Input read pairs | 14 |
| Donor rearrangement supported | 2 |
| Other host-mappable | 8 |
| Short donor/primer candidate | 2 |
| Other/unresolved | 2 |
| Supported junction sequence patterns | 4 |
| Direct-junction patterns | 3 |
| Intervening-sequence patterns | 1 |

The four junction patterns cover direct joining, 5-bp donor loss, a 3-bp MH-compatible junction and a 3-bp intervening sequence. Each has two supporting pairs and is retained at support thresholds of one and two. Intervening-sequence MH is unassessed (NA).

## Statistical units and interpretation

- Four-way proportions use **all input read pairs** as their denominator.
- Junction summaries count supported **candidate sequence groups**, not validated independent integration events. Supporting reads do not establish independent molecules.
- The descriptive primary setting is `nuclear_all`, `min_support=1`, `support_basis=all_supported_pairs`; support of at least two pairs is a sensitivity analysis. Keep 3′ and 5′ libraries separate.
- MH fractions use direct junctions; intervening-sequence candidates are excluded from that denominator. MH does not identify a repair pathway.
- Donor-terminal loss is distinct from host-genome deletion. Equivalent boundaries are retained as ranges.

Full thresholds, coordinate conventions and definitions are in [Methods](METHODS.md).

## Reproducing study results

A study rerun requires the paired study FASTQs, the exact host reference, the fixed donor/primer profile, and sample-to-condition/biological-replicate metadata. A summary-only rerun requires the completed `batch_summary.json` and all three matching merged junction tables, as described in the README. Sample labels are not used to infer biological conditions.

[Validation](VALIDATION.md) distinguishes recorded Linux engine validation from the additional submission-package checks. The reported 48-library aggregation regression used study tables that are not bundled and cannot be independently repeated from this archive alone.

File identities are recorded in `SHA256SUMS.txt` and `release_manifest.json`. On Linux/WSL, verify the extracted source before generating outputs:

```bash
sha256sum -c SHA256SUMS.txt
```
