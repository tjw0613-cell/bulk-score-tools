# Public API and stable return schema

`from bulk_score_tools import analyze_bulk_scores, load_paper_reference, load_paper_signatures`

`analyze_bulk_scores(bulk, q, signatures, *, reference=None, method='known_loading', bulk_scale='cpm', alpha=0.0, full_human_totals=None, q_provenance='unspecified', q_scale='rna', reference_law='product', invalid_policy='raise', compute_leakage=True, max_reference_states=100000)`

`bulk` and `q` are pandas DataFrames with samples as rows and unique, nonmissing IDs. `q` is reordered by IDs to match bulk; neither rows nor values are silently removed or renormalized. Gene symbols and compartment names are exact. A correlation requires at least three common samples and nonzero numerical variance. `alpha` is in per-million units (0.01 means log2(CPM + 0.01)).

Custom `signatures` is an insertion-ordered mapping with exactly two entries, e.g. `{'F': {'target':'Fibroblast', 'genes':['G1'], 'weights':[0.5], 'offset':-1.0}, 'M': {...}}`. Weights are strictly positive fixed coefficients in the given gene order (or a mapping of exact gene names to coefficients). Missing weights use 1 / gene count. The score is `offset + sum(weight * log2(expression_per_million + alpha))`. No data-dependent weights are fitted. Shared genes remain and are reported.

Preset strings P1-P6 load the frozen Bill coefficients and offset and explicitly attach the bundled Quah reference unless leakage is disabled. Bundled compartments are Fibroblast, Macrophage and Epithelial. Presets are fixed-membership rescorings, not native MCP-counter/Danaher algorithms.

Custom reference is a dictionary: `profiles` (donor x compartment x selected gene float probability array), `donors`, `compartments`, `genes` (unique nonmissing labels), `normalization='full_human_probability'`, nonempty `version`, and a nonempty `provenance` dictionary. It may additionally carry `full_human_totals`, `counts`, `n_cells`. The selected probabilities are not renormalized. All q compartments must match the reference exactly. A reference with no full-human denominator must not be labelled as this normalization.

Return dictionary (`schema_version='1.0'`):

- `correlations`: keys `raw` and `corrected`; each has `r` (float or None), `status` (`ok` or `undefined`), `reason` (None or a diagnostic string), `n` (common support).
- `scores`: DataFrame in original bulk row order; columns `<signature>__raw`, `<signature>__corrected`. Excluded rows have NaN for both methods.
- `valid_mask`: Boolean Series on original sample index.
- `exclusions`: DataFrame with `sample_id`, `reasons` (semicolon-separated); no exclusions yields an empty table with these columns.
- `coverage`: dictionary with `n_input`, `n_used`, `n_excluded`, `selected_gene_count`, `signature_gene_counts`.
- `gene_leakage`: DataFrame with `sample_id`, `signature`, `target`, `gene`, `lambda`, `valid_states`, `total_states`. Gene lambda is the mean of state-specific non-target / total ratios on positive-total states only. Rows with zero valid states have NaN, not zero. Excluded input samples do not have leakage rows.
- `score_leakage`: DataFrame with `sample_id`, `signature`, `target`, `weighted_sum`, `normalized_mean`, `beta`, `defined_gene_count`, `total_gene_count`, `status`. Any undefined gene makes the complete-score index undefined; no renormalization over available genes.
- `pair_leakage`: DataFrame with `sample_id`, `weighted_pair_mean`, `normalized_pair_mean`, `status`; arithmetic mean of the two complete-score indices.
- `metadata`: `method`, `alpha`, `bulk_scale`, `normalization_source`, `q_scale`, `q_provenance`, `signature_order`, resolved `signatures` (genes, weights, offset, beta), `shared_genes`, `preset`, `leakage_status`, `reference` (version, hash, law, donor/compartment/state counts, provenance), `ols` (or None), and interpretation notes. Custom signatures without a reference yield `leakage_status='not_estimated'` and empty leakage tables; disabled leakage yields `disabled`.
- `metadata['ols']` for joint OLS: `design_columns`, `rank`, `residual_df`, `rank_deficient`, and `residual_variance` keyed by signature. Fitted OLS residuals do not inherit the pointwise known-loading guarantee. Fixed q may produce rank 1 and still return defined correlations if residual variance and degrees of freedom remain positive.

The function raises `ValueError` for invalid units, negative/nonfinite inputs or derived scores, duplicate or unmatched labels, missing selected genes, malformed signatures/reference, and too many exact reference states. With `invalid_policy='complete_case'`, only zero selected expression at alpha=0 and nonpositive target q may be excluded, using one shared mask. It never rescales q, imputes genes, adds an epsilon, or infers full-human totals from selected columns. Correlations first rescale each score by the larger of 1 and its maximum absolute value to avoid overflow. Corrected scores also use the original raw-score scale to detect cancellation remnants. The centered norm must exceed 64 times float64 machine epsilon times the larger of 1 and the rescaled uncentered norm. Thus an exact design-span score does not acquire a spurious correlation from floating-point residuals.

Leakage indices are reference-conditioned descriptive outputs. They are not high/low risk categories, calibrated probabilities, patient-specific measured leakage, or proven predictors of correlation error. No guarantee is inferred merely from supplying numeric q.
