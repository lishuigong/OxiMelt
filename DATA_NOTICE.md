# Data and scientific provenance

OxiMelt distributes the fixed calculation inputs for the analytical relation described in the associated manuscript.

- `simple_oxide_reference_states.csv` contains the 57 simple-oxide reference states.
- `element_properties.csv` contains the five elemental properties used by the global relation.
- `training_pair_landscape.csv` identifies cation pairs represented in the two-cation fitting domain; these labels do not refit the equation.
- `model.json` contains the fixed coefficients, evaluation metrics and family-rule settings.
- `hf_local_residual_support.csv` contains seven binary chemical-system out-of-fold residual supports for Lu, Sc and Yb. The dopant-valence and nominal vacancy-fraction fields are derived from stoichiometry under formal Hf(+4), dopant(+3) and O(-II) assignments.

Three-cation target temperatures are not used to construct the binary support library or supplied as prediction inputs. The family correction applies only within the composition and support conditions stated in the manuscript. Equal-vacancy targets do not independently test concentration scaling.

The elemental properties originate from the Oliynyk elemental-property dataset cited in the manuscript. Reference temperatures and support records originate from the experimental sources documented in the manuscript and Supplementary Information.

The MIT license applies to the software code. Third-party data remain subject to their source terms; the software license does not grant additional rights to those data.
