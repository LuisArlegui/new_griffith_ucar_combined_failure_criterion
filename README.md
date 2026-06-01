# New Griffith–Úcar Combined Failure Criterion

![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/status-pre--release-orange)
![DOI](https://img.shields.io/badge/DOI-TO_BE_ADDED-lightgrey)

---

## Overview

This repository contains the R scripts accompanying the manuscript:

**A century after Griffith: towards a combined non-linear criterion for brittle failure**

**Authors:**


- Roberto Úcar
- Luis Arlegui
- Norly Belandria
- José Luis Simón

The repository implements a combined rock-failure criterion integrating:

- The **New Griffith Criterion (NGC)** in the tensile and transitional domains.
- The **Úcar Criterion (UC)** in the compressive domain.

The resulting envelope provides a continuous representation of rock failure from tensile conditions to high confinement.

---

## Main capabilities

The scripts allow the user to:

- Generate complete NGC–UC failure envelopes.
- Export failure curves in both:
  - Principal stress space (σ₃–σ₁).
  - Mohr space (σₙ–τ).
- Compare the criterion with laboratory triaxial tests.
- Analyse depth-dependent fracture domains.
- Simulate geothermal reservoir stimulation by thermal cooling.
- Calculate thermal failure boundaries for different thermoelastic parameters.

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

# Recommended workflow

## Step 1. Generate the combined NGC–UC curve

Run:

```r
01_generate_combined_NGC_UC_curve.R
```

User inputs:
Los valores que introduciremos en esta sección pueden provenir directamente de ensayos de compresión uniaxial C0 y de tracción indirecta o directa, T0, however, when triaxial laboratory tests are available, the recommended workflow is:

1. Fit the Úcar criterion to the triaxial dataset.
2. Obtain adjusted values of:

```text
C0
T0
k1
k2
```

3. Use those adjusted values as inputs for Script 01.

The recommended fitting procedure is the orthogonal non-linear regression methodology developed in the scripts accompany the publication:

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

- imports triaxial test pairs (σ₃, σ₁);
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

- reads the NGC–UC envelope;
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
critical cooling ΔT
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
- αT = thermal expansion coefficient

The resulting diagrams provide a rapid assessment of geothermal stimulation potential.



---

# Scientific background

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

# Citation

If you use these scripts, please cite:

```text
TO_BE_UPDATED_AFTER_ACCEPTANCE
```

and the associated Zenodo release.

---

# License

This project is distributed under the MIT License.

See:

```text
LICENSE
```

for details.

---

# Status

Current repository status:

```text
Private pre-release version
Prepared during peer review
Public release planned after acceptance
Zenodo DOI to be added after publication
```
