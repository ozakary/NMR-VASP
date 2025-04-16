# NMR Parameters Calculator from VASP

**Author:** Ouail Zakary  
**ORCID:** [0000-0002-7793-3306](https://orcid.org/0000-0002-7793-3306)  
**E-mail:** [Ouail.Zakary@oulu.fi](mailto:Ouail.Zakary@oulu.fi)  
**Website:** [Ouail Zakary - webpage](https://cc.oulu.fi/~nmrwww/members/Ouail_Zakary.html)  
**Academic Portfolio:** [Ouail Zakary - academic portfolio](https://ozakary.github.io/)  

A Python tool to extract and calculate Nuclear Magnetic Resonance (NMR) parameters from VASP output files.

![Version](https://img.shields.io/badge/version-0.1.0-blue)
![Python](https://img.shields.io/badge/python-3.6%2B-blue)
![VASP](https://img.shields.io/badge/VASP-6.4.1-green)
![License](https://img.shields.io/badge/license-MIT-orange)

## Overview

The VASP NMR Parameter Calculator extracts and processes NMR data from VASP output files (OUTCAR) that the standard VASP processing doesn't provide directly. This tool is particularly useful for researchers working with NMR computations in material science.

The calculator determines:
- Non-diagonalized total magnetic shielding tensor for each atom
- Principal axis system (PAS) tensor and principal components (σ<sub>11</sub>, σ<sub>22</sub>, σ<sub>33</sub>)
- Isotropic magnetic shielding (σ<sub>iso</sub>)
- Chemical shift anisotropy (σ<sub>CSA</sub>)
- Asymmetry parameter (η<sub>CSA</sub>)
- Span and skew
- Quadrupolar parameters (C<sub>Q</sub> and η<sub>Q</sub>) if the EFG tensor is computed and present in the OUTCAR file

## Tested With

- VASP version 6.4.1
- Python 3.6+

## Dependencies

- NumPy (>= 1.20.0)
- Python Standard Library (argparse, re, os, sys, csv)

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ozakary/NMR-VASP.git
cd vasp-nmr-calculator
```

### 2. Install dependencies

```bash
pip install numpy
```

### 3. Make the script executable (optional)

```bash
chmod +x nmr_calculator.py
```

## Usage

Run the calculator with your POSCAR and OUTCAR files:

```bash
python nmr_calculator.py --poscar /path/to/POSCAR --outcar /path/to/OUTCAR
```

Or if you made it executable:

```bash
./nmr_calculator.py --poscar /path/to/POSCAR --outcar /path/to/OUTCAR
```

### Command-line Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--poscar` | Yes | | Path to the POSCAR file |
| `--outcar` | Yes | | Path to the OUTCAR file |
| `--report` | No | nmr_report.txt | Output file for detailed report |
| `--organized_report` | No | nmr_organized_report.md | Output file for organized markdown report |
| `--csv` | No | nmr_data.csv | Output file for CSV data |
| `--verbosity` | No | 1 | Verbosity level: 0=silent, 1=normal, 2=verbose |

### Example

```bash
python nmr_calculator.py --poscar POSCAR --outcar OUTCAR --report my_nmr_report.txt --verbosity 2
```

## Output Files

### 1. Detailed Report (nmr_report.txt)

Contains comprehensive information for each atom:
- Ion index and element symbol
- Coordinates
- Non-diagonalized and diagonalized tensors
- Principal components (σ<sub>11</sub>, σ<sub>22</sub>, σ<sub>33</sub>)
- Magnetic shielding parameters (σ<sub>iso</sub>, σ<sub>CSA</sub>, η<sub>CSA</sub>, span, skew)
- Quadrupolar parameters (C<sub>Q</sub>, η<sub>Q</sub>) if available

### 2. Organized Report (nmr_organized_report.md)

A tabular Markdown report with:
- Table of magnetic shielding parameters
- Table of quadrupolar parameters (if available)

### 3. CSV Data (nmr_data.csv)

Comma-separated values file containing all parameters for easy import into plotting or data analysis software.

## Data Extraction Details

The calculator extracts the following data from VASP output files:

- From POSCAR:
  - Element list and order
  - Atom coordinates

- From OUTCAR:
  - Symmetrized tensors
  - G=0 contribution to chemical shift
  - Core susceptibility and cell volume
  - Core NMR properties
  - Quadrupolar parameters (if present)

## Algorithm Overview

1. The script reads the POSCAR file to determine atom types and coordinates
2. It extracts all necessary data from the OUTCAR file
3. For each atom, it:
   - Calculates the total tensor including all contributions
   - Diagonalizes the tensor to obtain principal components
   - Computes derived parameters (σᵢₛₒ, σCSA, ηCSA, span, skew)
   - Extracts quadrupolar parameters if available
4. Results are written to the specified output files

## Acknowledgments

This code is based on work by O. Zakary (ouail.zakary@oulu.fi).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

For questions or issues, please contact:
- Original code author: Ouail Zakary (ouail.zakary@oulu.fi)
- GitHub issues: https://github.com/yourusername/vasp-nmr-calculator/issues
