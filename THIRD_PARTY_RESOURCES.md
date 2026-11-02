# Scope of MIT and third-party resources

The MIT License applies to the authors' original implementation, tests and
original documentation. It does not grant new rights to third-party source
datasets, third-party software, or published signature resources.

The package imports NumPy and pandas as dependencies and does not bundle their
source code. They retain their respective upstream licences.

`bulk_score_tools/data/quah_reference.npz` derives from the Quah GSE188737
reference and preserves 33 selected genes with complete-human RNA denominators.
`bulk_score_tools/data/paper_presets.json` combines source-derived gene
memberships with the study's fixed coefficients. The exact source names and
membership are recorded inside that file. `resource_provenance.json` records
the frozen input hashes and normalization.

Compact archived fixtures contain source-derived numerical data and provenance,
including Quah and Simkin mixture inputs. `examples/` contains deterministic
synthetic mixtures derived from the bundled reference. These materials retain
their source attribution and applicable data-use conditions; the software MIT
licence is not a blanket licence for them. Their redistribution conditions must
be verified before public publication of a data-inclusive release.

No third-party licence has been substituted, waived, or inferred from public
accessibility. Data accessions, normalization metadata and input hashes accompany the bundled resources.
