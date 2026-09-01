# NGC–Úcar Failure Criterion Toolkit

*Open-source implementation accompanying the article published in the* **Journal of Structural Geology**.

## 🚀 Launch the online application

The **NGC-Ucar Toolkit** can be used directly from any modern web browser without installing R.

**Online application**

https://luisarlegui-ngc-ucar-toolkit.share.connect.posit.cloud/

---

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-stable-brightgreen)](#)
[![Technical Appendix](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.21869784-blue)](https://doi.org/10.5281/zenodo.21869784)
[![Paper DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.jsg.2026.105782-blue)](https://doi.org/10.1016/j.jsg.2026.105782)
[![Online App](https://img.shields.io/badge/Launch-NGC--Ucar%20Toolkit-success)](https://luisarlegui-ngc-ucar-toolkit.share.connect.posit.cloud/)


---

## Overview

This repository contains the R scripts, experimental datasets, and interactive Shiny application accompanying the article:

**A century after Griffith: towards a combined non-linear criterion for brittle failure**

**Authors:**

- Roberto Úcar
- Luis Arlegui
- Norly Belandria
- José Luis Simón

**Journal:** *Journal of Structural Geology*  
**DOI:** https://doi.org/10.1016/j.jsg.2026.105782

The repository implements a combined rock-failure criterion integrating:

- The **New Griffith Criterion (NGC)** in the tensile and transitional domains.
- The **Úcar Criterion (UC)** in the compressive domain.

The resulting envelope provides a continuous representation of rock failure from tensile conditions to high confinement.

In addition to the reproducible R scripts, the repository includes the NGC-Ucar Toolkit, an interactive Shiny web application that allows users to generate the combined failure envelope, import experimental triaxial datasets, compare laboratory data with theoretical predictions, and export the calculated results.

---

## Main capabilities

The repository allows users to:

- Generate complete NGC–Úcar failure envelopes.
- Visualise the criterion in both principal stress (σ3–σ1) and Mohr (σn–τ) spaces.
- Import experimental triaxial datasets.
- Compare laboratory data with theoretical failure envelopes.
- Export curves and parameter tables.
- Analyse depth-dependent fracture domains.
- Simulate geothermal reservoir stimulation by thermal cooling.
- Calculate thermal failure boundaries.
- Explore all calculations interactively through the online NGC-Ucar Toolkit.
---


## Software requirements

The scripts were developed and tested using:

```text
R ≥ 4.2
```

Main packages:

```r
ggplot2
dplyr
readr
tidyr
purrr
```

Missing packages are automatically installed when required.

---


## Interactive Shiny application

An interactive implementation of the proposed criterion is freely available through the NGC-Ucar Toolkit, an online Shiny application accessible at [https://…](https://luisarlegui-ngc-ucar-toolkit.share.connect.posit.cloud/) . The source code is maintained in the accompanying GitHub repository.

### Input parameters

The left-hand panel contains all user-defined parameters.

**Rock strength**

The first section defines the uniaxial compressive strength ($C_0$) and the uniaxial tensile strength ($T_0$). These two parameters constitute the primary input required to construct the combined NGC–Úcar failure envelope.

**Curve density**

The second section controls the numerical resolution used to generate the NGC and Úcar envelopes. Increasing the number of points produces smoother curves and also increases the resolution of the CSV files that can be exported from the application.

### Experimental datasets

Experimental triaxial datasets can be imported directly into the application.

Each CSV file follows the NGC-Ucar Toolkit format:

```text
Dataset identifier
Units
sigma3,sigma1
Experimental values
```



**UC truncation**

The third section allows the user to activate or deactivate the ductile-transition truncation applied to the Úcar criterion and to modify the corresponding $\sigma_1 / \sigma_3$ limiting ratio.

**Plot appearance**

The final section provides options for modifying the colours used to display the NGC and Úcar curves.

### Output

The application generates:

- the combined failure envelope in Mohr space ($\sigma_n$ - $\tau$);
- the corresponding envelope in principal stress space ($\sigma_1$ - $\sigma_3$);
- a table containing the calculated coordinates of the generated curves, which can be exported as CSV files.


---
## Example datasets

The repository includes the four experimental triaxial datasets used in the manuscript examples.

These datasets can be:

- imported directly into the online app NGC-Ucar Toolkit,
- analysed using the accompanying R scripts,
- used to reproduce the figures presented in the paper.

They are located in the `data/` directory.

---

# Recommended workflow

## Step 1. Generate the combined NGC–Úcar curve

Run:

```r
01_generate_combined_NGC_UC_curve.R
```

User inputs:
The input values used in this section may be obtained directly from uniaxial compression ($C_0$) and direct or indirect tensile ($T_0$) tests.  however, when triaxial laboratory tests are available, the recommended workflow is:

1. Fit the Úcar criterion to the triaxial dataset.
2. Obtain adjusted values of:

```text
C0
T0
k1
k2
```

3. Use those adjusted values as inputs for Script 01.

The recommended fitting procedure is the orthogonal non-linear regression methodology developed in the scripts described in the publication:

> Úcar, R., Arlegui, L., Belandria, N., & Torrijo, F. (2025). Estimating rock strength parameters across varied failure criteria: Application of spreadsheet and R-based orthogonal regression to triaxial test data. Journal of Rock Mechanics and Geotechnical Engineering, 17(8), 4685–4699.
DOI: https://doi.org/10.1016/j.jrmge.2024.11.024



This approach minimises orthogonal distances to the experimental data and generally provides more robust estimates than ordinary least squares.
```r
sigma_t_adj
sigma_c_adj
```

Outputs:

```text
results/tables/curve_NGC_rigorous_new.csv
results/tables/curve_UC.csv
results/tables/curve_combined_NGC_UC.csv
results/tables/parameters_NGC_UC.csv
```

and the corresponding figures.

This script constitutes the core of the repository and generates all parameters subsequently used by the remaining examples.

---

## Step 2. Compare the criterion with triaxial tests

Run:

```r
02_NGC_UC_triaxial_tests.R
```

This script:

- imports triaxial test pairs ($\sigma_1$ - $\sigma_3$);
- plots them against the combined envelope;
- generates Mohr circles;
- exports comparison figures and tables.

This script is useful both for educational purposes and for validation against laboratory datasets.

---

## Step 3. Example 01: Depth-dependent fracture domains

Run:

```r
03_example_01_depth_stress_fracture_domains.R
```

The script:

- reads the NGC–Úcar envelope;
- computes stress states as a function of depth;
- incorporates lithostatic loading;
- incorporates tectonic stress;
- classifies stress states as:

```text
No failure
Pure extensional fracture
Transitional / hybrid shear fracture
```

Outputs include fracture intervals, critical depths and transition points.

---

## Step 4. Example 02: Cooling-induced geothermal failure

Run:

```r
04_example_02_geothermal_cooling_failure.R
```

The script evaluates:

- thermal contraction;
- thermoelastic stress changes;
- cooling required to reach failure.

Outputs include:

```text
critical cooling $\Delta T$
failure mechanism
stress path
Mohr-circle evolution
```

This example reproduces the geothermal-reservoir application discussed in the manuscript.

---

## Step 5. Example 03: Thermal failure boundary

Run:

```r
05_example_03_thermal_failure_boundary.R
```

This script calculates:

```text
ΔTcrit = f(E, αT)
```

where:

- E = Young's modulus
- $\alphaT$ = thermal expansion coefficient

The resulting diagrams provide a rapid assessment of geothermal stimulation potential.



---

# Scientific framework

The combined criterion merges:

## Tensile domain

New Griffith Criterion (NGC)

including:

- pure extensional fractures;
- transitional fractures;
- hybrid tensile-shear failure.

## Compressive domain

Úcar Criterion (UC)

including:

- low confinement;
- intermediate confinement;
- high confinement.

The formulation provides a continuous failure envelope across tensile and compressive regimes.

---
## Repository structure

app/                 Interactive Shiny application

R/                   Computational functions and example scripts

data/                Experimental datasets

examples/            Stand-alone application examples

results/             Generated figures and tables

---

# Citation

If you use this repository or the NGC-Ucar Toolkit, please cite:

Úcar, R., Arlegui, L., Belandria, N., & Simón, J.L. (2026).  
A century after Griffith: towards a combined non-linear criterion for brittle failure.  
*Journal of Structural Geology*.  
https://doi.org/10.1016/j.jsg.2026.105782

### Related resources

**Technical Appendix**  
https://doi.org/10.5281/zenodo.21869784

**GitHub repository**  
https://github.com/LuisArlegui/new_griffith_ucar_combined_failure_criterion

**Online NGC-Ucar Toolkit**  
https://luisarlegui-ngc-ucar-toolkit.share.connect.posit.cloud/

---

# License

This project is distributed under the MIT License.

See:

```text
LICENSE
```

for details.

---


