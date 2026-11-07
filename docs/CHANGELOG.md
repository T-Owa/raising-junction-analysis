# Public release 1.0.0

The public release number describes packaging, documentation and aggregation portability; it does not replace the historical per-engine analysis versions recorded in outputs.

- Included original four-way 3′ engine 1.0.0 and 5′ engine 1.1.0, pilot 0.1.0 and refinement 0.1.1.
- Added explicit four-way batch `--allow-index-build`, forwarded to the existing engine option for fresh analysis roots. Classification thresholds and category priority are unchanged.
- Made Figure aggregation accept arbitrary library counts/sample labels and either compressed or decompressed merged tables. Explicit metadata replace name-derived experimental labels; conflicting values fail.
- Retained original numerical Figure definitions. Aggregation rejects unsupported multi-template equivalence groups instead of silently counting them as independent events.
- Added a synthetic demonstration, portable test runner, dependency check and source-sharing guidance.

This release has different four-way source hashes from earlier study runs. Use a new analysis root when rerunning with this version. Completed merged junction tables may instead be summarized directly with the standalone exporter without resuming or modifying the original engine runs.
