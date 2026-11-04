# Code availability — manuscript template

Use this text only after replacing the bracketed fields with the actual public resources:

> Custom scripts for paired-read classification, donor–host junction reassessment and figure-data aggregation are available at [repository URL]. The version used for this study (public release [version], commit [commit identifier]) is archived at [version-specific DOI]. Installation instructions, analysis definitions and a synthetic demonstration are included in the repository.

This archive has not been uploaded to a public repository and has no assigned DOI. Public author attribution can be added without including private workstation paths or instrument identifiers. Complete `CITATION.cff.example` with the intended public authors and repository URL, then rename it `CITATION.cff`.

Recommended release procedure:

1. Review the source-only archive, README, license and intended public attribution.
2. Upload the extracted source directory to the chosen repository, keeping work/results/raw-data folders outside the repository.
3. Tag the release used for the manuscript, archive it in an appropriate persistent software repository, and record the version-specific DOI.
4. Cite that exact release in the manuscript and describe the actual reference FASTA and software environment separately.

Official guidance:

- Nature Portfolio reporting/code-availability policy: https://www.nature.com/nature/editorial-policies/reporting-standards
- GitHub guidance on citation and Zenodo archiving: https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content
- Zenodo software metadata: https://help.zenodo.org/docs/github/describe-software/

Journal requirements differ; publishing a source archive does not replace raw-data accession, reference provenance or the manuscript's Methods.
