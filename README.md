# Satellite Accretion in the Local Group


This repository contains a **Jupyter notebook** analyzing the accretion histories and properties of satellite galaxies around Milky Way (MW) and M31 analogs in the **HESTIA Local Group simulations**.

**Dupuy et al., 2022 – _“Anisotropic satellite accretion onto the Local Group with HESTIA”_** ([arXiv:2208.14648](https://arxiv.org/abs/2208.14648))

![Satellite Accretion Visualization](accretion.png)

---

## Overview

The notebook explores:

- **Satellite infall histories** onto MW and M31 analogs.
- **Radial vs. tangential infall paths** via angles between position and velocity vectors.
- **Distance traveled** from birth location to host at first virial crossing.
- **Alignment with the cosmic web** using tidal and shear eigenvectors.
- **Dependence of infall alignment** on host galaxy, satellite mass, and infall redshift.
- **Accretion maps in host eigenframes**, visualized as octants and full-sky Aitoff projections.

> **Note:** Original datasets are proprietary and not included.

---

## Skills Highlighted

This project demonstrates:

- Handling **structured astrophysical datasets** from multiple simulations
- Computing **satellite infall statistics** and derived quantities
- Applying **vector algebra** to study radial infall and cosmic-web alignments
- Visualizing **histograms, scatter plots, and density maps** with Matplotlib
- Working with **host-centered and cosmic-web reference frames**
- Comparing **observed distributions to isotropic random expectations**

---

## Notebook Structure

1. **Loading HESTIA Satellite Infall Data**
   - Reads infall catalogs for three simulations (`09_18`, `17_11`, `37_11`)
   - Splits satellites by host galaxy (MW vs M31) and baryonic content
   - Allows selection at different virial radii (R200, 2×R200)
   - Prints summary statistics of satellite types

2. **Accretion Histories**
   - Plots **infall redshift distributions** per simulation and combined
   - Examines **radiality of infall paths** using cos(angle between position and velocity)
   - Computes **distance traveled** from birth to infall
   - Scatter plots compare **distance vs. infall redshift**, distinguishing DM-only vs baryonic satellites

3. **Cosmic Web Alignment**
   - Computes **|cos(theta)| between satellite infall directions and local eigenvectors**
     (tidal or shear tensors)
   - Compares to **random isotropic distributions** with 1σ and 2σ bands
   - Studies **dependence on host galaxy** (MW vs M31)
   - Explores **mass dependence** and **redshift dependence** of alignment
   - Visualizes results in histograms with quantitative significance

4. **Accretion Maps in Host Eigenframe**
   - Rotates satellite positions into host-centered eigenframe
   - Generates **octant projections** showing preferential infall directions
   - Creates **full-sky Aitoff projections** of satellite accretion
   - Highlights anisotropies and filamentary infall patterns

---

## Running the Notebook

1. Clone the repository:

```bash
git clone https://github.com/aadupuy/satellite-accretion.git
cd satellite-accretion
```

2. Install required packages

```bash
pip install numpy pandas matplotlib astropy
```

3. Launch the notebook

```bash
jupyter notebook satellite-accretion.ipynb
```

---

## Data Preprocessing & Tensor Eigenvector Extraction

### Satellite Infall Catalogs

Satellite infall catalogs were generated with the Python script `extract_satellite_infall.py`, which:
- Identifies host halos (MW/M31) and their satellites
- Detects satellite infall at a specified multiple of the host virial radius
- Outputs satellite positions, velocities, and masses at infall and at z=0

> **Note:** Original simulation outputs and infall datasets are proprietary and not included.

### Tensor Eigenvector Extraction

Tensor eigenvectors (shear and tidal) were extracted from HESTIA simulation grids using two scripts:

1. **`extract_tensor_eigenvectors.py`**  
   - Extracts eigenvectors **at satellite infall positions only**  
   - Focuses on MW/M31 host halos  
   - Output shape: `(N_sat, 13)`  

2. **`extract_tensor_eigenvectors_birth.py`**  
   - Extracts eigenvectors at **both infall and birth positions**  
   - Handles **multiple rvir selections** (0.5–2.0 × rvir)  
   - Processes all satellites in the infall catalog  
   - Output shape: `(N_sat, 19)`  

> **Note:** Original tensor grid files are proprietary and not included.