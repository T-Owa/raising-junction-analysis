# Source anonymization and analysis provenance

This distribution contains source code, the study's synthetic donor/primer profile, documentation and synthetic tests. It omits original sample-specific results, raw reads, vendor reports, local usernames/hostnames, original absolute paths and old execution logs. Example sample names do not represent study participants or experimental records.

The workflow intentionally records paths, input identities, source hashes, tool versions and read identifiers in **new local analysis outputs**. These records protect reproducibility and safe resumption. They have not been disabled to make the source distribution anonymous.

Consequently, the source distribution is ready for review as code, but newly generated results are **not automatically anonymized**. Keep raw FASTQs, SAM/BAM files, evidence tables, manifests and logs out of a public code repository. The supplied ignore file provides common exclusions; it does not retroactively remove files already tracked in Git.

For a source-data deposit, prepare a separate curated copy. Remove irrelevant local paths, use approved sample aliases and evaluate whether sequence/read-level data may be shared. Preserve a private mapping and the original provenance files. Do not modify the working manifests to evade their input/source checks.

Names of scientific methods, references, primers and external software are retained because they are necessary to understand and reproduce the analysis. Author attribution and identifiers for a public release can be filled in separately from workstation identifiers.
