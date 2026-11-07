# Analysis definitions

This document describes the bundled downstream analysis. The provider's wet-laboratory RAISING protocol and proprietary annotation pipeline are outside its scope. Reference identity and the validated software environment are described in [Reference provenance](REFERENCE_PROVENANCE.md) and [Validation](VALIDATION.md).

## Input and technical processing

R1/R2 FASTQ pairs are processed in their supplied order. The workflow validates paired identifiers, counts and sequence/quality lengths, rejects repeated normalized read identifiers, and hashes the input files. All input pairs are processed in resumable chunks; splitting into chunks is not sampling. Each 3′ or 5′ library is analyzed separately. The sample identifier determines the donor side using the `_3prm` or `_5prm` marker; experimental condition and biological replicate are not inferred by the public Figure exporter.

The final donor-specific primer (P3) is recognized near the R1 start with up to two mismatches. The appropriate P3 sequence and its reverse complement are considered by the screening helper; primer-aligned donor-reference bases do not contribute to independent donor-junction anchors. A poly(A) tract is technically trimmed only when its adjacent complete library adapter is recognized within the permitted mismatch tolerance. R2 start adapter/poly(T) features and low-quality 3′ edges are handled separately. Runs of at least ten identical bases are masked for host evidence. Original sequence and trimming coordinates are retained in local evidence outputs.

## Four-way classification

Cleaned sequences are aligned to the study donor with BLAST+ (`blastn-short`, word size 7, dust/soft masking disabled), and each mate is independently aligned to a combined host-plus-donor reference with BWA-MEM (`-a -Y -T 20`). R1/R2 evidence is combined at the **read-pair** level; the BWA calls are not conventional jointly aligned paired-end mappings.

Donor rearrangement support requires qualifying donor-alignment blocks with at least 30 independent matching bases on each side and at least 95% identity. Final-primer matches and shared query overlap are excluded from the independent-anchor counts. The screen permits up to 50 bp between donor blocks and up to 20 bp of overlap and requires the specified reference-order/orientation discontinuity (including the 20-bp rewind criterion for applicable same-orientation joins). These are candidates for donor–donor joining or donor rearrangement, not proof of a chromosome-integrated concatemer.

Accepted host evidence requires mapping quality of at least 20 (excluding unavailable MAPQ=255), at least 95% estimated identity and at least 30 independent high-quality matching host bases after excluding donor-covered and low-complexity positions. MD/CIGAR information is used to count matching bases; where an MD tag is absent, a conservative NM-based count is used. Secondary alignments are excluded; supplementary evidence can contribute. An accepted host alignment is not itself a fully resolved donor–host junction.

Each pair is assigned once, in this priority order:

| Machine-readable category | Short Figure label | Definition |
|---|---|---|
| `donor_rearrangement_supported` | Concatemer candidates | At least one mate supports a donor junction/rearrangement, including pairs that also support host mapping |
| `host_mappable_other` | Other host-mapped | Host-mapping support without a higher-priority donor-rearrangement call |
| `short_donor_primer_candidate` | Short donor/primer candidates | Neither preceding category; R1 matches the defined donor/primer-prefix followed by poly(A) pattern |
| `other_unresolved` | Unclassified | None of the preceding criteria are met |

The short-product pattern uses an exact donor/primer prefix of at least 15 nt, beginning at P3 or up to two bases into P3, followed by the first qualifying poly(A) run of at least 15 nt. Adapter confirmation is recorded separately and is not required for every short-product candidate. This class does not establish free donor, a donor monomer, absence of concatemerization or a biological donor end. The unresolved class is not a confirmed non-concatemer group. The supplementary both-supported/host-supported flags overlap categories and must not be summed with the four exclusive categories.

The four-way denominator is all input read pairs in each library, including unresolved pairs. PCR-amplified read abundance does not measure cells, independent molecules or absolute integration frequency.

## Donor–host junction discovery

The pilot uses previously accepted host evidence to obtain local genomic comparison windows, then searches for a directly sequenced donor-to-host transition within an individual cleaned read. It requires independent exact donor and host anchors of at least 20 and 30 bp, respectively, with quality and complexity checks and overlap with the original accepted host alignment. Shared MH and primer bases do not supply independent anchors. Up to 20 bp of intervening sequence are permitted. Unresolved or ambiguous evidence is retained rather than assigned zero MH or zero deletion.

Junction coordinates are zero-based boundaries **between bases**. Intervals are zero-based, half-open unless an output explicitly says otherwise. Donor loss is measured relative to the relevant outward-facing end of the supplied donor. Sequence-compatible alternative boundaries arising from MH are retained as minimum/maximum donor-loss and host-boundary values. These bounds are not statistical confidence intervals. The current fixed primer/anchor geometry gives nominal direct-junction donor-loss limits of 39 bp for 3′ libraries and 45 bp for 5′ libraries. These are analysis ascertainment limits, not biological maximum deletions.

An intervening sequence is a sequence-description candidate. It is not automatically nontemplated insertion or a proven physical insertion event. Host-genome deletion length cannot be determined from a single donor–host boundary; the opposite host boundary and a compatible structural reconstruction would be required.

## Existing-template reassessment

The reassessment stage organizes nearby host-boundary candidates (default cluster span 10 bp) and compares supporting reads to existing junction templates using the bundled Biopython alignment procedure. Bounded substitutions and small indels are permitted under the fixed settings (default identity at least 95%, at most three edited bases and three gap bases, quality threshold 20, boundary guard 3). Related templates are compared over common endpoints so that unequal comparison spans do not create spurious sequence ambiguity.

This step reassesses support for **existing** templates; it does not perform unrestricted discovery of novel junctions, reconstruct whole insertions, merge both donor ends into one insertion or automatically relabel every alternative deletion/insertion representation as a substitution. Alternative representations are exported for review. Template clustering is not proof of an independent insertion locus or clone.

The original `pair_support` field is retained as provenance. Final support is taken from `supported_pairs` in the sequence-group table or the corresponding `supported_pair_count_once_per_group` in the reassessment table, not by adding old and new counts. Each library's identical candidate sequence group is counted once for Figure summaries. The public exporter explicitly rejects multi-template equivalent groups rather than arbitrarily choosing among conflicting annotations; the validated study tables contain single-member groups.

## Figure denominators and subsets

The primary descriptive setting is `nuclear_all`, support at least one pair, with a support-at-least-two sensitivity analysis. These are **supported candidate sequence-pattern counts**, not independently validated insertion-event counts. Nuclear/unplaced reference contigs are included; mitochondrial calls are tracked separately. Three subsets are exported:

- `nuclear_all`: threshold on total `supported_pairs`.
- `nuclear_with_simple`: threshold on `without_prior_dd_or_complex_host_pairs`; a group may also contain donor-rearrangement/complex support.
- `nuclear_only_simple`: no supported prior-donor-rearrangement or complex-host pairs; threshold on simple-support pairs.

“Simple” refers to the available support classification; it does not prove that the entire integrated molecule is a monomer.

For MH proportions, the denominator is supported **direct-junction** patterns, excluding intervening-sequence candidates. MH=0, MH=1 and MH≥2 categories partition that denominator. The interfaces of insert-containing junctions have MH **not assessed (NA)** and are never included as MH=0. MH is a sequence feature and cannot alone establish a DNA-repair pathway.

For donor-terminal-loss proportions, the same direct-junction denominator is used. Minimum loss >0 denotes a required donor loss under all equivalent boundaries; minimum=maximum=0 denotes exact zero; minimum=0 and maximum>0 denotes boundary ambiguity. Both endpoint distributions and endpoint medians are exported, including zero-valued direct junctions. Intervening-sequence-candidate prevalence instead uses **all supported patterns** as denominator. Its length distribution uses the intervening-sequence subset as denominator.

Library fractions are not pooled read-level biological replicates. Compare independently obtained biological samples using the experimental design; 3′ and 5′ libraries from one sample remain paired technical views of that sample. A failure to resolve a junction is missing structural information, not absence of integration.
