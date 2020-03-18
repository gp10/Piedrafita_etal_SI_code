# Piedrafita_etal_SI_code :: pCellDynasts
Analyzing proliferating Cell Dynamics in squamous epithelial tissues.

Computational methods used to unravel the paradigm of cell proliferation in squamous epithelial tissues. This set of methods can be potentially extended to disentangle proliferating cell hierarchies and the mode of cell renewal in other tissues.

Check the original source for context details:
  > Piedrafita G, Kostiou V, Wabik A, Colom B, Fernandez-Antoran D, Herms A, Murai K, Hall BA, Jones PH (2020) A single-progenitor model as the unifying paradigm of epidermal and esophageal epithelial maintenance in mice. _Nat. Commun_ 11, 1429. https://doi.org/10.1038/s41467-020-15258-0

### Graphical abstract
![GraphicalAbstract](https://github.com/gp10/Piedrafita_etal_SI_code/blob/master/Graphical_abstract_pCellDynasts.png)

### Overview
This code combines computational methods used for cell-fate model simulations, inference and fitting of experimental data, organized in three categories found as separate folders here:
- **Clonal-Lineage-Inference**: inference on clonal dynamics. Division and differentiation fates as well as other kinetic properties are derived from cell lineage tracing.
- **H2BGFP-ModalityTests**: inference on cell proliferation heterogeneity through modality tests. Statistical analysis of multimodality in H2BGFP dilution patterns across individual cells (hints on division rate homogeneity).
- **H2BGFP-tcc-Inference**: inference on cell cycle time and division rate of proliferating cell populations. Cell division timing properties are resolved from H2BGFP intensity distributions.

Documentation (Readme files) on scripts and code used for specific purpose can be found in each corresponding folder.

### Requirements:
- Matlab (majority of code)
- Python
- R (for modality tests)

### Notes:
Figure references require update to agree with final order in publication. We expect to fix this within the coming few days, i.e. by end of March.
