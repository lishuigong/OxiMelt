# OxiMelt
OxiMelt is a lightweight, formula-only calculator for oxide melting temperature based on the frozen oxygen-referenced three-parameter analytical relation developed in this work. Version 1.1 also implements the frozen family-specific residual-transfer refinement for the predeclared HfO2-rich oxygen-vacancy three-cation family when binary out-of-fold support is sufficient.

The user enters an oxide chemical formula; OxiMelt automatically:
1. parses the composition,
2. assigns formal cation oxidation states using the conventional O(-II) oxidation-state convention,
3. selects the corresponding frozen simple-oxide melting-point references,
4. calculates the simple-oxide reference temperature `T0`,
5. calculates all pairwise `E_strain`, `E_chem`, and `Omega_ij` terms,
6. evaluates the frozen global 3P prediction,
7. checks applicability and the predeclared Hf defect-family rule,
8. when eligible and sufficiently supported, applies the frozen binary-OOF residual transfer,
9. returns both the global 3P value and the recommended final value.

No model retraining is performed and no three-cation target label is used to construct the local correction.

Global frozen equation

```text
Tm = T0 + sum_(i<j) xi*xj*[b0 + b1*E_strain,ij + b2*E_chem,ij]
```

with

```text
E_strain,ij = Kbar_ij*rbar_ij*(|ri-rj|/rbar_ij)^2
E_chem,ij   = (Ibar_ij - 3.5*Ecohbar_ij)*(|chi_i-chi_j|/(3.44-chibar_ij))^2
```

Frozen coefficients:

- `b0 = -792.5130673853114 K`
- `b1 = 18.945566079895755`
- `b2 = 233.31975001858584`
- `lambda = 3.5`
- `chi_O = 3.44`

## Hf-family local residual transfer

For the predeclared three-cation HfO2-rich oxygen-vacancy family, OxiMelt evaluates the same frozen rule used in the paper. Binary support residuals are experimental melting temperature minus **chemical-system out-of-fold** global-3P prediction. For each non-Hf dopant `D_k`, the binary residual is matched at the target Hf cation fraction by exact matching, conservative interpolation, or a very-near neighbor.

```text
alpha_k = x_Dk / (1 - x_Hf)
Delta_local = sum_k alpha_k * Delta_Hf-Dk^OOF
Tm_local = Tm_3P + Delta_local
```

The correction is applied only if:

- the target is a three-cation composition,
- `x_Hf >= 0.55`,
- Hf has formal oxidation state +4,
- the oxygen-deficiency proxy is at least `0.005`,
- at least one dopant has formal oxidation state below +4,
- frozen binary-residual support covers at least `80%` of the dopant share.

The published support-matching limits are frozen as:

- exact Hf-fraction tolerance: `0.01`
- maximum nearest Hf-fraction distance: `0.03`
- maximum interpolation span: `0.15`
- supported shares are **not** renormalized.

This is a family-specific transfer rule, not a general oxygen-vacancy thermodynamic correction.

Desktop GUI

Requires Python 3.10+.

Clone the repository using the HTTPS URL shown under GitHub's **Code** button, then run:

```bash
cd OxiMelt
python app.py
```

Windows users can also double-click:

```text
run_oximelt.bat
```

The interface reports:

- recommended predicted melting temperature in K and °C,
- global 3P prediction,
- family-specific local correction when applicable,
- prediction mode and support coverage,
- formal oxidation states and cation fractions,
- `T0` and pair-by-pair contributions,
- applicability-domain information,
- batch CSV prediction.

Command-line usage

After installation:

```bash
pip install -e .
oximelt predict Lu3Al5O12
```

Hf-family examples:

```bash
# two-cation OOF-calibrated family estimate
oximelt predict Hf0.6Lu0.4O1.8

# three-cation family-specific residual transfer
oximelt predict Hf0.6Lu0.2Sc0.2O1.8
```

Batch prediction:

```bash
oximelt batch examples/batch_example.csv predictions.csv
```

The input CSV must contain a `formula` column.

## Output semantics

- `global_Tm_K`: frozen global 3P prediction.
- `local_correction_K`: family-specific binary-OOF residual transfer; zero when not applied.
- `corrected_Tm_K`: `global_Tm_K + local_correction_K`.
- `predicted_Tm_K`: recommended output; identical to `corrected_Tm_K`.
- `prediction_mode`: indicates whether the result is global-only or family-corrected.
- `hf_local_support_coverage`: fraction of dopant share with valid frozen binary support.

For ordinary compositions, the local correction is zero and the recommended result equals the global 3P result.

For a supported two-cation Hf-defect composition, the local correction is taken from the frozen chemical-system OOF residual profile at the target Hf fraction (exact match, interpolation, or a predeclared nearest-support rule). This output is explicitly labeled **OOF calibrated**. Three-cation Hf-defect compositions retain the composition-weighted residual-transfer semantics.

## Applicability labels

- REFERENCE ONLY — single-cation oxide; output is the frozen simple-oxide reference.
- TWO-CATION DOMAIN — the cation pair occurs in the two-cation training chemistry.
- PAIR EXTRAPOLATION — a two-cation pair absent from the two-cation training chemistry.
- THREE-CATION ZERO-SHOT DOMAIN — global frozen zero-shot application.
- TWO-CATION Hf DEFECT FAMILY — OOF CALIBRATED — the Hf-defect rule is satisfied and a frozen binary OOF residual profile is available at the target composition. This is a family-calibrated estimate, not an independent external-validation prediction.
- THREE-CATION Hf DEFECT FAMILY — LOCAL TRANSFER — the predeclared defect-family rule is satisfied and binary residual support is sufficient.
- HIGHER-ORDER EXTRAPOLATION — mathematically evaluable, but >3 cations were not directly validated in the paper.
- UNAVAILABLE / outside current domain — no unique formal-valence assignment, missing frozen reference/property data, peroxide/superoxide-like chemistry, or unsupported composition.

### Important oxidation-state note

O(-II) is used only as a conventional formal oxidation-state convention for cation-valence assignment, reference-state matching, and electroneutrality. It is **not** interpreted as a fixed actual local charge of -2e in every oxide.

## Reproducibility resources

```text
data/model.json
data/element_properties.csv
data/simple_oxide_reference_states.csv
data/training_pair_landscape.csv
data/hf_local_residual_support.csv
```

`hf_local_residual_support.csv` contains only frozen binary chemical-system OOF residual support. For two-cation Hf-defect inputs it can be used as an OOF-calibrated family estimate; for three-cation Hf-defect inputs it is used for composition-weighted residual transfer. It is not generated from three-cation target labels.

The calculator performs no fitting and no network access.

## Tests

```bash
python -m unittest discover -s tests -v
```

Regression tests cover the global literature-consistency examples, formula parsing, rejection of unsupported peroxide chemistry, the Hf0.6Lu0.4O1.8 binary OOF-calibration case, all six published three-cation Hf-family local-transfer cases, and insufficient-support fallback behavior.

## Scientific status

Global model-level validation reported in the associated work:

- two-cation chemical-system Group-CV RMSE: 224.13 K
- three-cation frozen zero-shot RMSE: 164.81 K
- three-cation MAE: 133.37 K
- three-cation R²: 0.842

Family-specific refinement reported in the associated work:

- Hf defect-family subset, global 3P RMSE: 436.65 K (`n = 6`)
- same subset after frozen local residual transfer: 9.23 K
- all 176 three-cation records with the local rule applied only where eligible/supported: 143.75 K RMSE

The 143.75 K hybrid value is a layered analysis metric and does not replace the paper's primary global zero-shot metric of 164.81 K.

## Citation

Please cite the associated paper when using OxiMelt. Update `CITATION.cff` with the final bibliographic metadata after publication.

## License

Software source code is released under the MIT License.

Data notice: before public release, verify redistribution permissions for third-party elemental-property data from the original source. If redistribution is restricted, replace `data/element_properties.csv` with an installer/download step or distribute only data for which redistribution is authorized.
