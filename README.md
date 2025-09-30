# Mass-Luminosity Relations for Stripped Stars

Python scripts to estimate minimum, maximum, pure-helium mass-luminosity relations (MLRs) for fully and partially stripped stellar structure configurations.
Links:
https://arxiv.org/pdf/2508.14161
https://ui.adsabs.harvard.edu/abs/2025arXiv250814161S/abstract

# Overview

This repository provides tools to compute mass-luminosity relations for:
- **Fully stripped stars:** pure-helium composition
- **Partially stripped stars:** helium core + hydrogen shell with varying composition slope giving minimum mass for given luminosity or conversely maximum luminosity for given mass

- **Fully homogeneous stars:** for various amounts of hydrogen giving maximum mass for given luminosity or conversely minimum luminosity for given mass

The MLR estimates are derived from detailed structure model grids. For a comprehensive description of the underlying stellar models and methodology, see [Sabhahit et al. (2025b)](https://arxiv.org/pdf/2508.14161).

Additional script is also provided to compare Fully homogeneous and pure-He MLRs of Gräfener et al. (2011)

# Usage

Please select the required calculator and enter either stellar mass or luminosity, hydrogen and metal abundances as mass fractions. Running the code will provide the minimum, maximum, and pure-He values for the user input parameters. For more details regarding the structure model grid, see the text description at https://gautham-sabhahit.github.io/ml-calculator/

Grid parameter range:

1. For mass M, the chemically homogeneous structures (H profile slope of 0) and pure-He structures (H profile slope of inf) have the range 1 ≤ M/Msun ≤ 40, while the structures with H profile slope in between these two extremes have the range 1 ≤ M/Msun ≤ 18.
2. For surface H mass fraction, the range is 0 ≤ X ≤ 0.7
3. For surface metal mass fraction, the values are Z = 0.008 (LMC-like, 0.4Zsun) and Z = 0.004 (SMC-like, 0.2Zsun) where Zsun = 0.02.

Warnings and Errors:

1. Errors are displayed if the inputs are not valid numbers, or if the mass is zero or negative, or if X or Z is negative. X = 0 and Z = 0 are allowed.
2. A set of warnings is printed based on the parameter range of the synthetic model grid. If the inputs fall outside the grid’s tested parameter range, a general warning is shown. If the inputs are significantly beyond the grid range such that the minimum or maximum value of M or L is not truly a minimum or maximum, then a warning is issued indicating that the ML fits may be unreliable. If a calculation fails, especially in the mass calculator, an error is issued.
3. The model grid was computed for Z = 0.008 and Z = 0.004. For any Z value other than 0.008 or 0.004, interpolation or extrapolation is performed, and a corresponding warning is provided.

## Citation

If you use these scripts in your research, please cite:
