# Measurement of the cancer effects of somatic double-base substitutions

This repository contains the complete analysis for _Testing of custom evolutionary models of oncogenesis : Extension to double-base substitutions_, a poster to be presented at SMBE2023 (Society for Molecular Biology & Evolution). The analysis relies on our **cancereffectsizeR** package ([website](https://townsend-lab-yale.github.io/cancereffectsizeR/index.html) | [code](https://github.com/Townsend-Lab-Yale/cancereffectsizeR/) | [publication](https://aacrjournals.org/cancerres/article/83/4/500/716429/Estimation-of-Neutral-Mutation-Rates-and)).

To completely reproduce the analysis, install the exact versions of cancereffectsizeR and the reference data package (ces.refset.hg38) used for the poster. Alternatively, future versions of these packages (3.0.0 or later for cancereffectsizeR, 1.3.0 or later for the data package) will likely work, too. These updates will be released sometime in late 2023.

```r
options(timeout = 600)
install.packages('remotes')
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@SMBE23", dependencies = TRUE)
remotes::install_github("Townsend-Lab-Yale/ces.refset.hg38@SMBE23")

# Install additional packages required for the analysis:
install.packages('ggstatsplot')
install.pacakges('irr')
```

After the installations, restart your R session, and then run through [dbs_analysis.R](dbs_analysis.R) to acquire data and reproduce the analysis.

We welcome any feedback!
