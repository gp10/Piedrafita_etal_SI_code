# H2BGFP_tcc_Inference
Computational methods to estimate the individual cell-cycle period distribution (tcc) from experimental H2BGFP dilution data.

### Dependencies
- H2BGFPdil-tccDist-inference.m : main script to infer tcc, by comparing simulated to experimental data.
- H2BGFP-Int-OE-dataset.mat : contains the experimental data on H2BGFP intensities in esophagus.
- MonteCarloSimulator-BasalCell-H2BGFPdil.m : runs simulations of basal cell turnover and H2BGFP dilution.
- ABCrejection-tccDist.m : runs ABC method for tcc parameter estimation (based on experimental vs. simulation goodness-of-fit).

### Requirements
Matlab R2016b