# NMR Parameters Analysis from VASP
---
- 📄 Author: **Ouail Zakary**  
- 📧 Email: [Ouail.Zakary@oulu.fi](mailto:Ouail.Zakary@oulu.fi)  
- 🔗 ORCID: [0000-0002-7793-3306](https://orcid.org/0000-0002-7793-3306)  
- 🌐 Website: [Personal Webpage](https://cc.oulu.fi/~nmrwww/members/Ouail_Zakary.html)  
- 📁 Portfolio: [Academic Portfolio](https://ozakary.github.io/)
---
A Python tool to extract and calculate Nuclear Magnetic Resonance (NMR) parameters from VASP output files.

![Version](https://img.shields.io/badge/version-0.2.0-blue)
![Python](https://img.shields.io/badge/python-3.6%2B-blue)
![VASP](https://img.shields.io/badge/VASP-6.4.1-green)
![License](https://img.shields.io/badge/license-MIT-orange)

## Overview

The VASP NMR Parameter Calculator extracts and processes NMR data from VASP output files (OUTCAR) that the standard VASP processing does not provide directly. This tool is particularly useful for researchers working with NMR computations in materials science.

The calculator determines:
- Non-diagonalized total magnetic shielding tensor for each atom
- Principal axis system (PAS) tensor and principal components (σ<sub>11</sub>, σ<sub>22</sub>, σ<sub>33</sub>)
- Isotropic magnetic shielding (σ<sub>iso</sub>)
- Chemical shift anisotropy (σ<sub>CSA</sub>)
- Asymmetry parameter (η<sub>CSA</sub>)
- Span (Ω) and skew (κ)
- Quadrupolar parameters (C<sub>Q</sub> and η<sub>Q</sub>) if the EFG tensor is computed and present in the OUTCAR file

## Conventions

All parameters follow established solid-state NMR conventions. The principal components are ordered in the **shielding convention**:

**σ<sub>11</sub> ≤ σ<sub>22</sub> ≤ σ<sub>33</sub>**

where σ<sub>11</sub> is the least shielded and σ<sub>33</sub> is the most shielded direction. This is the shielding analogue of the standard IUPAC shift convention (δ<sub>11</sub> ≥ δ<sub>22</sub> ≥ δ<sub>33</sub>).

### Haeberlen convention
The chemical shift anisotropy σ<sub>CSA</sub> is the **Haeberlen reduced anisotropy**:

**σ<sub>CSA</sub> = σ<sub>zz</sub> − σ<sub>iso</sub>**

where σ<sub>zz</sub> is the principal component furthest from σ<sub>iso</sub>. σ<sub>CSA</sub> is positive when σ<sub>33</sub> dominates and negative when σ<sub>11</sub> dominates. The full anisotropy is Δσ = (3/2) σ<sub>CSA</sub>.

The asymmetry parameter η<sub>CSA</sub> ∈ [0, 1] is:

**η<sub>CSA</sub> = (σ<sub>yy</sub> − σ<sub>xx</sub>) / σ<sub>CSA</sub>**

where σ<sub>xx</sub> and σ<sub>yy</sub> are the two remaining principal components, ordered by proximity to σ<sub>iso</sub>.

### Herzfeld-Berger convention
The span and skew are defined as:

**Ω = σ<sub>33</sub> − σ<sub>11</sub>** (Ω ≥ 0)

**κ = 3(σ<sub>iso</sub> − σ<sub>22</sub>) / Ω** (−1 ≤ κ ≤ +1)

## Tested With
- VASP version 6.4.1
- Python 3.6+

## Dependencies
- NumPy (>= 1.20.0)
- Python Standard Library (`argparse`, `re`, `os`, `sys`, `csv`)

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/ozakary/NMR-VASP.git
cd NMR-VASP
```

### 2. Install dependencies
```bash
pip install numpy
```

### 3. Make the script executable (optional)
```bash
chmod +x nmr_vasp.py
```

## Usage

Run the calculator with your POSCAR and OUTCAR files:
```bash
python nmr_vasp.py --poscar /path/to/POSCAR --outcar /path/to/OUTCAR
```

Or if you made it executable:
```bash
./nmr_vasp.py --poscar /path/to/POSCAR --outcar /path/to/OUTCAR
```

### Command-line Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--poscar` | Yes | | Path to the POSCAR file |
| `--outcar` | Yes | | Path to the OUTCAR file |
| `--report` | No | `nmr_report.txt` | Output file for detailed report |
| `--organized_report` | No | `nmr_organized_report.md` | Output file for organized markdown report |
| `--csv` | No | `nmr_data.csv` | Output file for CSV data |
| `--verbosity` | No | `1` | Verbosity level: 0=silent, 1=normal, 2=verbose |
| `--sanity_check` | No | `False` | Compare computed parameters against VASP's own reported values |
| `--tolerance` | No | `0.0001` | Tolerance in ppm for the sanity check (default: 4 decimal places) |

### Examples

Basic run:
```bash
python nmr_vasp.py --poscar POSCAR --outcar OUTCAR
```

With custom output files and verbosity:
```bash
python nmr_vasp.py --poscar POSCAR --outcar OUTCAR --report my_nmr_report.txt --verbosity 2
```

With sanity check:
```bash
python nmr_vasp.py --poscar POSCAR --outcar OUTCAR --sanity_check
```

## Output Files

### 1. Detailed Report (`nmr_report.txt`)
Contains comprehensive per-atom information:
- Ion index and element symbol
- Cartesian coordinates
- Non-diagonalized and diagonalized total shielding tensors
- Principal components (σ<sub>11</sub>, σ<sub>22</sub>, σ<sub>33</sub>)
- Magnetic shielding parameters (σ<sub>iso</sub>, σ<sub>CSA</sub>, η<sub>CSA</sub>, Ω, κ)
- Quadrupolar parameters (C<sub>Q</sub>, η<sub>Q</sub>) if available
- Convention notes explaining all parameter definitions

### 2. Organized Report (`nmr_organized_report.md`)
A tabular Markdown report with:
- Convention summary
- Table of all magnetic shielding parameters
- Table of quadrupolar parameters (if available)

### 3. CSV Data (`nmr_data.csv`)
Comma-separated values file containing all parameters for easy import into plotting or data analysis software.

## Sanity Check

The `--sanity_check` flag activates a validation routine that compares the computed σ<sub>iso</sub>, span, and skew against the values VASP itself reports in the `CSA tensor` block of the OUTCAR (specifically the *absolute, valence and core, including G=0 contribution* section). All three quantities should be **identical to 4 decimal places**.

```
==========================================================================================
SANITY CHECK: Computed vs. VASP-reported (absolute, valence+core, incl. G=0)
  Note: VASP ISO_SHIFT sign is reversed; reference values shown here as -ISO_SHIFT.
==========================================================================================
 Ion  Elem    σ_iso calc   σ_iso VASP    Δσ_iso   Span calc  Span VASP   ΔSpan  ...  Status
------------------------------------------------------------------------------------------
   1     C       69.6759      69.6759    0.0000    214.8361   214.8361   0.0000  ...      OK
 ...
==========================================================================================
✓ All atoms are identical to 4 decimal places.
==========================================================================================
```

If mismatches are found, the output flags the affected atoms and suggests likely causes (tensor format issue, wrong conversion factor, or sign error).

## Data Extraction Details

The calculator extracts the following data from VASP output files:

**From POSCAR:**
- Element list, order, and counts (used as the authoritative atom count)
- Atom coordinates (Direct or Cartesian, with scaling factor applied)

**From OUTCAR:**
- Symmetrized shielding tensors (correctly distinguished from the unsymmetrized tensors that precede them in the OUTCAR)
- G=0 contribution to the chemical shift
- Core susceptibility and cell volume
- Core NMR properties (element-specific core shielding contributions)
- Quadrupolar parameters C<sub>Q</sub> and η<sub>Q</sub> (if present)
- VASP's own reported σ<sub>iso</sub>, span, and skew (for sanity checking)

## Algorithm Overview

1. The POSCAR file is read to determine atom types, counts, and coordinates.
2. All necessary data are extracted from the OUTCAR file. The symmetrized tensor block is carefully distinguished from the unsymmetrized tensor block that shares a similar keyword.
3. For each atom, the total shielding tensor is assembled as:

   **σ<sub>total</sub> = σ<sub>symm</sub> + σ<sub>G=0</sub> + σ<sub>core,NMR</sub> + σ<sub>core,suscept</sub>**

   A sign correction is applied to account for VASP's internal sign convention.

4. The tensor is symmetrized and diagonalized. Principal components are sorted as σ<sub>11</sub> ≤ σ<sub>22</sub> ≤ σ<sub>33</sub>.
5. Derived parameters (σ<sub>iso</sub>, σ<sub>CSA</sub>, η<sub>CSA</sub>, Ω, κ) are computed following the Haeberlen and Herzfeld-Berger conventions.
6. Quadrupolar parameters are extracted if available.
7. Results are written to the specified output files.

## What's New in v0.2.0

- **Sanity check** (`--sanity_check`, `--tolerance`): validates computed parameters directly against VASP's own reported values, atom by atom.
- **Fixed σ<sub>CSA</sub>**: now correctly computed as the Haeberlen reduced anisotropy σ<sub>zz</sub> − σ<sub>iso</sub>, where σ<sub>zz</sub> is the component furthest from σ<sub>iso</sub>. Can be positive or negative. The previous version always used σ<sub>33</sub> and was off by a factor of 3/2.
- **Fixed η<sub>CSA</sub>**: denominator now always pairs correctly with the σ<sub>CSA</sub> sign branch, guaranteeing η<sub>CSA</sub> ∈ [0, 1].
- **Fixed tensor block parsing**: the symmetrized tensor block is now correctly identified, avoiding a substring match against the unsymmetrized tensor block that caused atom counts to be doubled.
- **Convention notes** added to all output files.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
