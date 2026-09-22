# Data and scientific provenance notice

This repository contains the frozen calculation state used for the publication-specific OxiMelt implementation.

- `data/simple_oxide_reference_states.csv` is the final simple-oxide reference-state table used by the global analytical relation.
- `data/element_properties.csv` contains only the five elemental properties required by the frozen global equation.
- `data/training_pair_landscape.csv` is used only to label whether a pair occurred in the two-cation training chemistry; it does not refit the equation.
- `data/model.json` records the frozen coefficients, validation summary, and predeclared local-transfer thresholds.
- `data/hf_local_residual_support.csv` contains the binary chemical-system **out-of-fold** residual profiles used by the family-specific Hf local-transfer rule. No three-cation target melting-point label is used to construct this support table or the correction applied to a new formula.

The local-transfer support is specific to the HfO2-rich oxygen-vacancy family described in the associated paper and must not be interpreted as a universal defect correction.

Before uploading the repository publicly, confirm that redistribution of the elemental-property subset is permitted by the source dataset's license/terms. The MIT license in this repository applies to the software code, not automatically to third-party data.
