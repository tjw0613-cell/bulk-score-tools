These are compact fixtures extracted from eight frozen natural-profile simulation cohorts.
They are not new biological data or a new benchmark run.

Coverage: Quah and Simkin; paired and product donors; fixed_equal and variable_variable scenarios; cohort index0;118 mixture rows per case.

bulk_probability.csv contains33 selected genes already divided by the original complete-human compartment RNA totals before the original q-weighted mixing. Row sums are below1 because other human genes are not exported. Do not renormalize these selected columns.

rna_q.csv preserves original compartment labels. Quah uses Fibroblast/Macrophage/Epithelial; Simkin uses Fibroblast/Macrophage/Malignant. Correction parity runs use compute_leakage=False and do not silently relabel the Simkin background as the bundled Quah Epithelial reference. The separate six-preset leakage regression test uses Quah's own compartment axis and q=(.2,.2,.6), product donor law.

expected_correlations.csv records raw, known_loading and joint_ols separately. The source method label logq_offset is the archived fixed-loading oracle. The saved small-alpha field known_q_OLS means fitted joint OLS, never fixed-loading correction. Atalpha.01 no saved offset was assumed; an independent direct subtraction gives the expectation.

Atalpha0, selected bulk zeros produce strict API errors. The fixture additionally records finite-source pure-target availability, which is a different property. Explicit complete-case checks use the same mask for both score sides and both raw/corrected outputs; they are not compared with the original full-support target.

The six fixed signature definitions include actual Bill coefficients and fixed centering. Saved correlation arrays omitted fixed centering constants; Pearson correlations and OLS with an intercept are unchanged by that convention.

To validate from the software directory:
python validation/validate_archived_fixtures.py

Quickstart-compatible safe example: Quah_paired_variable_variable, preset P4, bulk_scale='probability', alpha=0. It retains the bundled Quah compartment labels and complete scoring support. For reference-diagnostic regression across all six presets, usealpha.01 so reference test inputs do not fail due to selected bulk zeros.
