# Bulk score tools

Version 0.1.0.

This package implements two fixed-score log-RNA corrections and an exact-enumeration, reference-conditioned gene leakage diagnostic. It packages the existing study calculations; it does not train a model or claim a new validated estimator.

## Quick start

Use Python 3.10 or newer with NumPy and pandas. From this extracted directory:

```bash
python -m pip install -r requirements.txt
python quickstart.py
python -m unittest discover -s tests -v
python validation/validate_archived_fixtures.py
```

The runtime and tests require no internet access once dependencies are installed. Optional installation is `python -m pip install .`; otherwise imports work from this directory. The declared dependency ranges describe the API requirements; the validation records identify the environment actually tested, not every version or platform.

`quickstart.py` reads `examples/bulk_cpm.csv` and `examples/q_rna.csv`. These are 12 deterministic synthetic mixtures assembled from the bundled Quah reference; they are not patients or new clinical validation. The example uses preset P4, known RNA fractions, the known-loading correction, and an explicitly chosen +0.01 CPM pseudocount. The **public function default remains alpha=0**. Outputs are written to `example_output/` and include raw/corrected scores, both correlations, complete exclusions, per-gene leakage, weighted signature indices and provenance. Override the output directory with `--output PATH`.

## Three main inputs

```python
import pandas as pd
from bulk_score_tools import analyze_bulk_scores

bulk = pd.read_csv('examples/bulk_cpm.csv', index_col='sample_id')
q = pd.read_csv('examples/q_rna.csv', index_col='sample_id')
result = analyze_bulk_scores(
    bulk, q, 'P4', alpha=0.01, bulk_scale='cpm', q_provenance='known'
)
print(result['correlations'])
print(result['gene_leakage'].head())
```

1. `bulk`: rows are sample IDs, columns are exact gene symbols; values are on a **linear** CPM, TPM, probability or counts scale. Supply the scale explicitly when it differs from CPM. Log expression is not an accepted input.
2. `q`: rows match the bulk IDs; columns contain the complete RNA-fraction composition, sum to one, and retain their original values. The two target fractions must be positive. Cell fractions are not RNA fractions and require an upstream RNA-content conversion. Label q as known, estimated or unspecified; a label itself does not establish estimation accuracy.
3. `signatures`: a paper preset P1-P6 or exactly two named fixed definitions. Custom definitions provide `target`, `genes`, optional fixed positive `weights`, and optional `offset`. Omitted weights mean equal-weight mean-log coefficients. They do not fit a new z score.

Presets use the existing pair-specific shared-member removal and frozen Bill coefficients and mean offsets. P1-P4 are the original four pairs; P5-P6 use the compatible Danaher macrophage list. These six combinations share two fibroblast and three myeloid lists and are not six independent biological validations or native implementations of the signature source methods. See `bulk_score_tools/data/paper_presets.json` within the package for every gene and coefficient. For custom signatures shared genes are reported and retained, never automatically removed.

Presets explicitly load the bundled Quah reference unless `compute_leakage=False`. It has 7 donors, compartments `Fibroblast`, `Macrophage`, `Epithelial`, and 33 selected genes with full-human denominators. Custom signatures do not receive an implicit reference: without an explicit reference their leakage status is `not_estimated`. There is no automatic synonym or compartment relabelling. For example, a `Malignant` background cannot be silently relabelled `Epithelial` to use this reference. Independent-source scoring can disable leakage and use its actual q labels.

## Correction definitions and limits

For fixed coefficients w, the raw score is `offset + sum(w_g * log2(CPM_g + alpha))`. The default correction is

`corrected = raw - sum(w_g) * log2(q_target)`.

With positive target-specific expression, no non-target expression, true positive RNA fractions, no pseudocount and matching fixed scoring, this known-loading subtraction recovers the corresponding pure-compartment score pointwise. It does not require independence of state and composition. Real data do not establish these assumptions merely by supplying q; leakage, estimated q and positive pseudocounts affect recovery. A positive pseudocount defines its own regularized score target.

Optional `method='joint_ols'` regresses each raw score on an intercept and both log2 target RNA fractions and returns the residuals. It reports the design rank, residual degrees of freedom and residual variance. It may remove state variation associated with composition and **does not share the pointwise guarantee** above. Constant fractions yield an intercept-only effective design; correlations remain unchanged if otherwise defined. Undefined correlations are returned as `None` with reasons, never as zero.

## Reference-conditioned leakage

For each input q row, target and gene, the diagnostic enumerates paired donor states or the Cartesian product of compartment donors. It calculates the non-target/total ratio **within each state before averaging** over positive-total states. A zero/zero state is undefined. Each result reports both the valid-state count and the full state count; a gene with no valid state stays undefined. This differs from a ratio of average reference expression. If any gene has undefined leakage, the complete-signature index remains undefined; weights are not renormalized over the remaining genes.

The score-level outputs are the actual weighted sum of gene means, its weight-normalized counterpart, and the mean of the two signature indices. They are descriptive indices, **not** error bounds, calibrated risk scores or patient-specific measured leakage. The weighted complete-signature extension did not establish a consistent cross-donor-rule diagnostic in the study. No high/low risk threshold is applied.

The default exact enumeration cap is 100,000 reference states. Exceeding it raises an error before enumeration; there is no stochastic approximation. The bundled product reference has 343 states. Valid-state counts can be gene-specific; undefined expression is not filled with zero. All selected probability profiles retain the complete human denominator.

## Counts, zeros and missing values

For counts, `full_human_totals` must be a Series indexed by the same sample IDs. It must contain full-human library totals, even if bulk includes only selected genes. The function never uses the selected-gene column sum as a full library size and does not renormalize a reduced matrix.

Zero selected expression at alpha=0 and nonpositive target q fail by default. Explicit `invalid_policy='complete_case'` removes their common row union, returns all excluded IDs/reasons, and uses the same retained rows for raw and corrected correlations. Negative/nonfinite expression, missing genes, duplicate IDs, incorrect q sums or units always fail. No epsilon, clipping of q, gene imputation or gene selection is performed. A +0.01 pseudocount means +0.01 CPM (not +0.01 to probabilities). A missing gene is not a zero count.

## Outputs and manuscript correspondence

`correlations` reports raw and corrected correlations, sample counts, and reasons for undefined results. `scores`, `valid_mask` and `exclusions` provide sample-level scores and the common retained sample set. `gene_leakage`, `score_leakage` and `pair_leakage` provide the reference-conditioned diagnostics. `metadata` records scales, weights, fraction provenance, reference identifiers and hashes, reference-state rules, and OLS diagnostics.

The function does not return confidence intervals. Uncertainty in estimated RNA fractions and transfer between reference sources is not automatically propagated into the correlations.

Manuscript Proposition 2 describes the known-loading correction; Supplementary Methods SM8 describes this interface, input scales and output definitions. `joint_ols` is the fitted adjustment described separately in the manuscript.

## Reproducibility and provenance

See `API_SCHEMA.md` for the exact input and return schema. Tests use deterministic gold calculations, equivalent input units, permutations, weight rescaling, masking, undefined states, paired/product references and cost guardrails. Run the tests and fixture replay commands above to check the implementation. VALIDATION.json summarizes the checks run on this package.

Bundled `bulk_score_tools/data/quah_reference.npz` includes selected probabilities, selected observed counts, full-human totals, donor and compartment labels, and cell counts. `resource_provenance.json` records original input SHA-256 hashes and normalization. These resources come from the existing frozen study inputs, not new annotation or processing. `paper_presets.json` records gene membership, coefficients, original shared-member removal and source list names.

The authors' original software and original documentation are licensed under the MIT License (see LICENSE). Bundled third-party reference data and source-derived signature resources retain their original conditions; MIT does not relicense those resources (see THIRD_PARTY_RESOURCES.md). Source repository: [https://github.com/tjw0613-cell/bulk-score-tools](https://github.com/tjw0613-cell/bulk-score-tools). Study archive DOI: [10.5281/zenodo.23886566](https://doi.org/10.5281/zenodo.23886566). No remote service, user credentials, machine-specific source-data paths or network requests are used at runtime.

`MANIFEST_SHA256.json` lists files and SHA-256 hashes for this software directory. The archive-root `FILE_MANIFEST_SHA256.csv` covers the full public archive; run `python software/verify_archive.py --strict` from the archive root to verify an unmodified extraction.
