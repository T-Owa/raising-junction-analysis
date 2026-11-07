# Reference provenance

## Donor profile

The bundled `RI_donor_800bp.fasta` is the 800-bp study donor. Both four-way engines carry the same sequence. The final primers are listed in `reference/primers.csv`; the implementation validates their occurrence and orientation. This profile is fixed in code. Replacing only the FASTA is not a supported way to analyze another donor or primer set.

## Host genome

Study analyses used a human GRCh38 primary-assembly FASTA. An Ensembl release-99 download URL was supplied during setup:

https://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

The exact release and checksum of the host FASTA used in the study are not recorded in this repository. The URL and filename alone do not establish the identity of that file. Exact reproduction therefore requires the retained study FASTA or its verified checksum and source metadata.

The provider's reported mapping reference (GRCh38.p13) and annotation database (GRCh38.99) describe its own workflow. They do not independently establish the exact FASTA used by the custom BWA workflow. Vendor annotations are not inputs to the custom junction caller.

The custom workflow appends the synthetic donor as `__RAISING_DONOR_800__` and builds a combined reference. Record the exact host FASTA bytes, source URL/release, donor sequence and checksums. Compare analyses only when these choices are compatible. Mitochondrial names `MT`/`chrM` are distinguished from the nuclear/unplaced subset; inspect other contig naming conventions before reuse.

An optional local record can be produced with standard tools, for example:

```bash
sha256sum ./work/ref/host.fa.gz
micromamba list -n raising-public --explicit > ./work/environment-explicit.txt
python tools/check_environment.py --output ./work/software-versions.json
```

Review environment exports before public sharing: channel URLs can contain local/private configuration. The bundled environment recipe pins the known study Python/Biopython/BWA/BLAST versions; a complete export of the original resolved environment is not included.
