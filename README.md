# XSBFEM-XGBoost Crack Frequency Prediction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-Pending-blue.svg)](https://doi.org/...) 

This repository contains the dataset, pre-trained models, and source code associated with the research paper: **"Data-Driven Framework of Vibrations in Bidirectional Gradient Plates: An Extended Scaled Boundary Element Method."**

Our work proposes a novel hybrid framework that integrates the numerical rigor of the Extended Scaled Boundary Finite Element Method (XSBFEM) with the rapid predictive capabilities of Machine Learning (ML), specifically a BCMO-tuned XGBoost model. 

## 🚀 Key Features
* **High-Fidelity Dataset:** A comprehensive dataset of 3,000 vibration scenarios for cracked Bidirectional Functionally Graded Material (BDFGM) plates, generated using our highly stable XSBFEM solver.
* **Pre-trained ML Models:** Includes the highly accurate BCMO-XGBR model (R² = 0.999), capable of real-time natural frequency prediction.
* **Efficiency Verification Data:** Additional datasets generated via Latin Hypercube Sampling (LHS) specifically designed to test and verify the computational efficiency and generalization of the models.
* **SHAP Interpretability:** Tools to quantify how boundary conditions, grading indices, and crack length influence dynamic responses.
* **XSBFEM Benchmark Code:** Open-source XSBFEM code for a classic fracture mechanics benchmark (an inclined center crack under uniform tensile loads), including the computation of Mode I and Mode II Stress Intensity Factors (SIFs).

## 📁 Repository Structure

```text
├── Data/
│   ├── Calculated_3000_LHS_Dataset.xlsx       # Original dataset generated via XSBFEM
│   └── New_3000_LHS_Dataset.xlsx      # LHS generated data for efficiency testing Best_BCMO_XGBoost_Model
├── Models/
│   ├── Best_BCMO_XGBoost_Model.mat       # The trained BCMO-XGBR model
│   └── ML_predic.m               # Load a new dataset for testing
│   └── v5_model.mat               # The trained BCMO-XGBR model Shape
│   └── ML_Mode_prediction.m               # Load a new dataset for testing
├── xsbfem/
├── main_inclined_crack.m          % Main driver (reproduces Fig. 18)
├── test_sbfem_pure_mode.m         % Standalone validation test
├── get_material_matrix.m          % 2D constitutive matrix
├── sbfem_coeff_matrices.m         % E0,E1,E2 for closed polygon
├── sbfem_coeff_matrices_crack.m   % E0,E1,E2 for crack (open boundary)
├── sbfem_stiffness.m              % SBFEM stiffness (closed polygon)
├── sbfem_crack_tip.m              % Crack tip stiffness + SIF extraction
├── create_graded_mesh.m           % Graded Q4 mesh generation
├── classify_elements.m            % Element classification using level sets
├── heaviside_element_stiffness.m  % Enriched stiffness for cut elements
├── Results/
│   └── SIF_comparisons/             # Analytical vs. Numerical SIF plots
├── requirements.txt
└── README.md

## How to Run

### Quick Validation Test
```matlab
cd xsbfem
test_sbfem_pure_mode
```
This validates the SBFEM SIF computation for pure mode I, II, and mixed mode 
cracks using analytical boundary conditions.

### Full Simulation
```matlab
cd xsbfem
main_inclined_crack
```
This runs the complete simulation:
1. **Part A:** Direct SBFEM validation with Williams expansion BCs
2. **Part B:** Full FEM + enriched SBFEM simulation of the plate

## Method Overview

### SBFEM Core
1. Boundary nodes of a polygon define the SBFEM subdomain
2. Scaling centre placed inside (or at crack tip for fracture)
3. Coefficient matrices $E_0$, $E_1$, $E_2$ computed from boundary integration
4. Hamiltonian matrix $Z$ formed
5. Schur decomposition gives modal quantities
6. Stiffness: 
   $$K = V_q V_u^{-1}$$

### Enrichment Strategy 
- **Standard FE elements:** Computed using conventional Q4 formulation
- **Heaviside enriched elements:** Elements completely cut by crack are split 
  into two SBFEM subdomains with virtual nodes at intersection points
- **Crack tip super element:** Multiple element layers around crack tip form 
  a single SBFEM subdomain with scaling centre at the tip

### SIF Computation 
Two methods implemented:
- **Displacement-based:** Uses crack opening displacement from 
  singular eigenvalue modes
- **Stress-based:** Interpolates singular stress modes to crack front

### Key Equations

- **Singular eigenvalues:** $-1 < \text{Re}(\lambda) < 0$ (corresponding to $r^{-1/2}$ singularity)
  

- **SIF from stress:** $\{K_I;K_{II}\}=\sqrt{2\pi L_0}\times\{\Psi_{yy}^s;\Psi_{xy}^s\}\times c^s$
  

- **SIF from displacement:** $\{K_I;K_{II}\}=\frac{G}{\kappa+1}\times\sqrt{\frac{2\pi}{r_0}}\times\{\Delta u_y;\Delta u_x\}$



![untitled1](https://github.com/user-attachments/assets/e7b87a86-3062-4d24-912f-856ee2d937c6)

<img width="1062" height="1134" alt="image" src="https://github.com/user-attachments/assets/dbf8b060-91dd-4ace-ad6d-4e6a1d23083c" />


<img width="1321" height="947" alt="image" src="https://github.com/user-attachments/assets/a28f3979-9a4b-4d18-a579-84c8f4577fb5" />


## Notes

- The mesh grading parameter (grade_power) concentrates elements near the 
  crack. Increase for better resolution.
- Super element requires 3-5 layers of elements around the crack tip.
- The Schur decomposition may require regularization for ill-conditioned 
  coefficient matrices.
- For best accuracy, increase the number of boundary nodes and Gauss 
  integration points.

## References

[1] Jiang S-Y, Du C-B, Ooi ET. Engineering Fracture Mechanics 222 (2019) 106734  
[2] Song C, Ooi ET, Natarajan S. Engineering Fracture Mechanics 187 (2018) 45-73  
[3] Wolf J, Song C. Comp Methods Appl Mech Engrg 190 (2001) 5551-5568
