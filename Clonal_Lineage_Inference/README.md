# Clonal_Lineage_Inference
Computational methods used for parameter inference to estimate the goodness-of-fit of simulated vs. experimental clone size distributions.

### Dependencies
- BasalCloneSizeDist-paramFit.m : main script used for model goodness-of-fit (gets log-Likelihood value and plots clone size distribution fits of specific parameter conditions).
- Lrig1-LineageTracing-dataset.mat : contains the experimental data on experimental clone sizes in esophagus.
- MonteCarloSimulator-SP-BasalCloneDynamics.m : runs simulations of basal-layer clone sizes over time.
- logLike-calc.m : computes the log-Likelihood match of simulated vs. experimental clone size distributions.
- size2freq.m : calculates the frequency histogram (distribution) of clone sizes from their individual sizes.
- size2freqbinned.m : calculates the frequency histogram (distribution) binned in categories increasing in size in powers of 2.

### Requirements
Matlab R2016b