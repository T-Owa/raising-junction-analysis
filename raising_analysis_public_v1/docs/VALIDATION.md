# Validation of public release 1.0.0

Validated on **2026-10-24**, Linux x86_64, Python **3.12.14**, Biopython **1.88**, BWA **0.7.19-r1273**, BLAST+ **2.17.0+** and samtools **1.24**. The installed compiled tools were used for the synthetic round trip; no external binaries are distributed in this source package.

## Synthetic/unit suites

All **107 tests passed**, with no skips:

| Suite | Tests |
|---|---:|
| Four-way batch, pairing, 5′ orientation, resumability and failure handling | 22 |
| Exact-anchor junction discovery | 26 |
| Existing-template reassessment, equivalent boundaries and alignment alternatives | 36 |
| Junction batch completeness and provenance checks | 7 |
| Figure aggregation, metadata conflicts, compressed input, empty libraries and invalid counts | 12 |
| Deterministic synthetic data generation and overwrite protection | 4 |

The portable runner uses isolated unittest discovery for every suite. This includes the reassessment test module, which has no standalone `__main__` entry point. Some tests deliberately trigger failure cases; their expected error messages are not failed test cases. Four-way unit tests mock external alignment calls where appropriate; these are separate from the actual-tool test below.

```bash
python tests/run_all.py
```

## Actual-tool synthetic FASTQ-to-Figure round trip

Passed with the actual BWA, BLAST+ and samtools executables and Biopython 1.88:

```bash
python tests/run_all.py --include-external
```

The demonstration generates a 100-kb random toy-host reference and two synthetic libraries (3′ and 5′), each containing 14 paired reads. Each library recovered the intended four-class counts:

| Category | Expected and observed pairs per library |
|---|---:|
| Donor rearrangement supported | 2 |
| Other host-mappable | 8 |
| Short donor/primer candidate | 2 |
| Other/unresolved | 2 |

Each library recovered four supported host-junction patterns with two read pairs per pattern: direct joining with no donor loss/MH, 5-bp donor loss, a 3-bp MH-compatible junction with donor-loss bounds 0–3, and a 3-bp intervening sequence whose MH remained unassessed. The primary all-supported nuclear subset recovered these at both support thresholds. Count denominators, required donor loss and insertion-interface MH=NA were explicitly checked after the complete pipeline.

An independent synthetic round trip with Biopython 1.85 also passed, including all 12 Figure-summary combinations across two libraries, three subsets and two thresholds. This does not establish equivalence for every possible dataset or every software version.

## Regression against existing study tables

The public aggregator was privately rerun on the existing 48-library completed junction tables. All eight CSV exports matched previous cell values exactly after excluding the five metadata-label columns that intentionally no longer infer labels from sample names:

| Export | Compared rows |
|---|---:|
| Figure values | 288 |
| Metrics | 288 |
| Library QC | 48 |
| MH distribution | 4,608 |
| Donor-loss distribution | 12,384 |
| Intervening-sequence length distribution | 5,760 |
| Alternative representation review | 944 |
| Candidate annotations | 101,337 |

This regression validates aggregation consistency, not the biological truth of each upstream candidate. Original data and private validation logs are deliberately excluded from the source distribution.

## Scope of the claims

- The junction engine source files were retained byte-for-byte from the previously used version.
- Four-way classification functions, thresholds, priority and reference sequences were preserved; changes are generic diagnostic text/defaults and an explicit batch option to request existing index-building functionality.
- All 31 included Python source files parse, and the donor reference and primer coordinates were checked.
- A source-only privacy audit found no original usernames, workstation paths, original sample identifiers or instrument read identifiers in the distributed files.
- The original full study FASTQs were **not** remapped in this build environment. Do not describe these tests as full biological validation, a measured false-positive rate or proof of complete sensitivity.
- The installation recipe was not solved as a complete micromamba environment during this build; the specified major analysis-tool versions were independently exercised. Retain a solved environment export for the actual manuscript run.

See `release_manifest.json` and `SHA256SUMS.txt` for the packaged file identities. If code, settings or references change after release, rerun the relevant checks and assign a new release version.
