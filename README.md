# OxiMelt

OxiMelt estimates oxide melting temperatures from chemical formulas using a frozen three-coefficient analytical relation. It combines a composition-weighted simple-oxide reference with effective cation-pair interactions. A nominal-vacancy-scaled residual correction is available for supported HfO2-rich compositions.

Associated manuscript: **From Binary to Multication Oxides: Composable Cation-Pair Interactions for Melting-Point Prediction**.

## Run

Python 3.10 or later is required. The calculation uses the Python standard library; the graphical interface also requires Tkinter.

```bash
python app.py
python -m oximelt.cli predict Lu3Al5O12
python -m oximelt.cli batch examples/batch_example.csv predictions.csv
```

On Windows, `run_oximelt.bat` opens the graphical interface. On Linux or macOS, run `bash run_oximelt.sh`. The command line works without Tkinter. For an editable installation, use `pip install -e .`; the entry points are `oximelt` and `oximelt-gui`.

Batch input must contain a `formula` column. Decimal coefficients, Unicode subscripts and parentheses are supported. Use a single stoichiometric formula, such as `Hf0.6Sc0.4O1.8`. Mixture notation such as `0.25Sc2O3-0.75HfO2` identifies a support record but is not an accepted input formula.

## Global equation

```text
T_global = T0 + sum_(i<j) x_i*x_j*Omega_ij
Omega_ij = b0 + b1*E_strain,ij + b2*E_chem,ij
E_strain,ij = Kbar_ij*rbar_ij*(abs(r_i-r_j)/rbar_ij)^2
E_chem,ij = (Ibar_ij - lambda*Ecohbar_ij)
            * (abs(chi_i-chi_j)/(chi_O-chibar_ij))^2
```

The fixed coefficients are `b0 = -792.5130673853114`, `b1 = 18.945566079895755`, and `b2 = 233.31975001858584`, with `lambda = 3.5` and `chi_O = 3.44`. The coefficient units follow the descriptor definitions in the manuscript. Input properties and oxide references are supplied in `data/`. No fitting or network access is performed during prediction.

## Hf-family correction

```text
fV = (2 - n_O / sum_cations(n_i)) / 2
   = sum_cations[x_i * (4 - z_i)] / 4
alpha_k = x_Dk / (1 - x_Hf)
Delta_local = sum_supported_k alpha_k * (fV_target / fV_reference_k) * residual_OOF_k
T_family = T_global + Delta_local
```

The signed OOF residual is the experimental support temperature minus its chemical-system out-of-fold prediction. The target must contain two or three cation elements, have `x_Hf >= 0.55`, Hf at formal valence +4, at least one lower-valence dopant, and `fV > 0`.

For each dopant, the eligible same-element, same-valence support with the closest nominal vacancy fraction is selected. Duplicate fractions use the median signed residual; an equal-distance tie selects the lower fraction. The residual is scaled by the vacancy-fraction ratio and weighted by the dopant fraction on the non-Hf sublattice.

Support coverage must be at least 0.80. Unsupported terms contribute zero, without renormalizing the supported weights. Compositions failing the family or coverage criteria use the global estimate. Targets outside a dopant's support range are labeled concentration extrapolation. Nearest-support selection is not interpolation and can introduce discontinuities when the selected support changes. Compositions with four or more cation species use only the global equation and are labeled higher-order extrapolation.

Nominal vacancy fractions are stoichiometric proxies, not measured defect concentrations. The residual does not isolate a vacancy entropy or formation-energy contribution. Independent evidence across vacancy concentrations is needed to establish the accuracy of proportional scaling.

## Outputs

The interface displays a selected estimate and retains the global result, pair contributions and support information for inspection. CSV and JSON outputs contain:

- `global_Tm_K/C`: global three-parameter prediction.
- `local_correction_K`: applied family correction; zero if the rule is inactive.
- `corrected_Tm_K/C`: global prediction plus the applied correction.
- `predicted_Tm_K/C`: selected estimate, numerically equal to the corrected fields.
- `nominal_oxygen_vacancy_site_fraction`: target stoichiometric vacancy proxy.
- `hf_local_support_coverage`: supported fraction of the non-Hf cation sublattice.
- `hf_local_correction_applied`: whether the correction is applied.
- `hf_local_concentration_extrapolation`: whether an available support is used outside its dopant's vacancy-fraction range; inspect the applied flag separately if coverage is insufficient.
- `hf_local_support_details_json`: per-dopant support composition, raw and scaled residuals, vacancy fractions, ratio and extrapolation status.

Within support details, `residual_K` is the scaled residual and `unscaled_residual_K` is the signed support residual. `distance` is the absolute Hf-fraction difference; `vacancy_fraction_distance` is the support-selection distance. `support_x_Hf_low` and `support_x_Hf_high` are null because the rule does not interpolate between supports.

The selected family estimate is not a guarantee of improved accuracy for an individual composition. O(-II) is a formal-valence convention, not a measured atomic charge. Missing properties, unavailable reference states or ambiguous formal-valence assignments prevent prediction.

## Evaluation

The manuscript reports a two-cation chemical-system Group-CV RMSE of 224.13 K and a three-cation zero-shot RMSE of 164.81 K, MAE of 133.37 K and R² of 0.842.

For the six literature Hf-family targets, target and support vacancy fractions are all 0.10, so every scaling factor equals one. The family RMSE is 9.23 K; applying this correction to the six supported targets gives an overall hierarchical RMSE of 143.75 K across the 176 three-cation records. This hierarchical metric is separate from the global zero-shot result. The six targets and their supports originate from the same study and do not independently validate concentration scaling.

`verification/three_cation_predictions.csv` contains the 176 record-level calculations. `verification/validation_summary.json` summarizes their metrics.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests cover global examples, the six literature family targets, synthetic concentration cases, formula parsing, signed residuals, support selection, coverage, domain fallback, CSV export and the GUI prediction callback. The callback test does not open a desktop window.

## Data and citation

`data/` contains the coefficients, elemental properties, simple-oxide references, training-pair labels and seven binary OOF supports. Support valence and nominal vacancy fractions are calculated from stoichiometry.

Please cite the associated manuscript when using OxiMelt. Citation metadata are provided in `CITATION.cff`. The code is MIT licensed; third-party data provenance is described in `DATA_NOTICE.md`.
