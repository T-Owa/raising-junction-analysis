# Analysis outputs and sharing

This repository contains source code, the study donor/primer profile, documentation and synthetic tests. Study sequencing reads and sample-specific results are not included. Demonstration sample names are synthetic.

The workflow records paths, input identities, source hashes, tool versions, read identifiers and sequences in local analysis outputs. These records support reproducibility and validation of resumed runs. Keep the original manifests with the corresponding analysis results.

Before sharing outputs, review logs, manifests and read-level evidence tables for local paths and sample or sequence information that is not intended for release. Prepare a separate sharing copy while preserving the original provenance records. The supplied `.gitignore` excludes common analysis products but does not remove files already tracked in Git.
