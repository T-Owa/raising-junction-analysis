# Code availability and citation

This repository provides the source code for paired-read classification, donor–host junction discovery and reassessment, and figure-data aggregation. Installation instructions and a synthetic demonstration are in the [README](../README.md); analysis definitions are in [Methods](METHODS.md).

Repository: [T-Owa/raising-junction-analysis](https://github.com/T-Owa/raising-junction-analysis). The source archive can be inspected and run without a GitHub account.

## Citation

Identify the software as **RAISING downstream analysis**, and cite the repository URL together with the version and commit used for the analysis. When reviewing an archive without Git history, `SHA256SUMS.txt` and `release_manifest.json` identify its exact file contents. The package version is recorded in `VERSION`; engine versions are listed in the README and `release_manifest.json`.

If the version used has an archived software DOI, include that version-specific DOI. No DOI is recorded in this repository's metadata.

## License and data

The source code is distributed under the [MIT license](../LICENSE). Study sequencing reads and sample-specific results are not included. The bundled donor/primer references and synthetic demonstration support execution of the workflow; exact reproduction of a study run also requires its inputs, reference FASTA, settings and software environment. See [reference provenance](REFERENCE_PROVENANCE.md) and [validation](VALIDATION.md).
