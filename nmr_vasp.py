#!/usr/bin/env python3
###################################################################################################################################################################################################
########################################################< VASP NMR Parameter Calculator >##########################################################################################################
###################################################################################################################################################################################################
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This Python package calculates NMR parameters from VASP output files, including:
# - Non-diagonalized total magnetic shielding tensor for each atom
# - Principal axis system (PAS) tensor and principal components (sigma_11, sigma_22, sigma_33)
# - Isotropic magnetic shielding (sigma_iso)
# - Chemical shift anisotropy (sigma_csa)
# - Asymmetry parameter (eta_csa)
# - Span and skew
# - Quadrupolar parameters (C_Q and eta_Q) if present in the OUTCAR file
#
# Usage: nmr_calculator.py --poscar <POSCAR_FILE> --outcar <OUTCAR_FILE> --report <REPORT_FILE> --organized_report <ORG_REPORT_FILE> --verbosity <1 or 2>
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

import re
import sys
import numpy as np
import argparse
import os
import csv

class NMRCalculator:
    def __init__(self, outcar_file, poscar_file, verbosity=1):
        self.outcar_file = outcar_file
        self.poscar_file = poscar_file
        self.verbosity = verbosity
        self.results = []
        self.atom_coordinates = []
        
    def extract_symmetrized_tensors(self):
        """Extract symmetrized tensors from OUTCAR file."""
        symm_tensors = []

        with open(self.outcar_file, 'r') as file:
            lines = file.readlines()

        # Locate the section containing the symmetrized tensors
        # NOTE: 'SYMMETRIZED TENSORS' must be matched carefully because
        # 'UNSYMMETRIZED TENSORS' contains it as a substring.
        start_line_idx = None
        for idx, line in enumerate(lines):
            if 'SYMMETRIZED TENSORS' in line and 'UNSYMMETRIZED' not in line:
                start_line_idx = idx + 1
                break

        if start_line_idx is None:
            if self.verbosity > 0:
                print("SYMMETRIZED TENSORS section not found.")
            return []

        lines = lines[start_line_idx:]
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith('ion'):
                i += 1
                tensor_data = []

                for _ in range(3):
                    if i < len(lines):
                        tensor_row = [float(x) for x in lines[i].strip().split()]
                        tensor_data.append(tensor_row)
                        i += 1
                    else:
                        break

                if len(tensor_data) == 3:
                    symm_tensors.append(np.array(tensor_data, dtype=float))
            else:
                i += 1

        return symm_tensors

    def extract_nmr_values(self):
        """Extract core NMR properties from OUTCAR file."""
        nmr_values = {}
        found_block = False
        start_processing = False
        
        with open(self.outcar_file, 'r') as file:
            for line in file:
                line = line.strip()
                
                if 'Core NMR properties' in line:
                    found_block = True
                    continue
                
                if found_block:
                    if '----------------------------' in line:
                        if not start_processing:
                            start_processing = True
                        else:
                            break
                    elif start_processing:
                        parts = line.split()
                        if len(parts) >= 3:
                            element = parts[1]
                            try:
                                value = float(parts[2])
                                nmr_values[element] = np.diag([value, value, value])
                            except ValueError:
                                if self.verbosity > 1:
                                    print(f"Skipping line due to ValueError: {line}")
        
        if not found_block and self.verbosity > 0:
            print("'Core NMR properties' block not found in the file.")
        
        return nmr_values

    def extract_quadrupolar_parameters(self):
        """Extract quadrupolar parameters (C_Q and eta_Q) from OUTCAR file."""
        quadrupolar_params = []
        found_block = False
        start_processing = False
        
        with open(self.outcar_file, 'r') as file:
            for line in file:
                line = line.strip()
                
                if 'NMR quadrupolar parameters' in line:
                    found_block = True
                    continue
                
                if found_block:
                    if '----------------------------------------------------------------------' in line:
                        if not start_processing:
                            start_processing = True
                        elif start_processing and len(quadrupolar_params) > 0:
                            # End of the block
                            break
                    elif start_processing and line and not line.startswith('ion'):
                        parts = line.split()
                        if len(parts) >= 4:
                            try:
                                ion_idx = int(parts[0])
                                c_q = float(parts[1])
                                eta_q = float(parts[2])
                                q_mb = float(parts[3])
                                quadrupolar_params.append({
                                    'ion': ion_idx,
                                    'Cq': c_q,
                                    'eta_Q': eta_q,
                                    'Q_mb': q_mb
                                })
                            except (ValueError, IndexError):
                                if self.verbosity > 1:
                                    print(f"Skipping quadrupolar line due to parsing error: {line}")

        if not found_block and self.verbosity > 0:
            print("'NMR quadrupolar parameters' block not found in the file.")
        
        return quadrupolar_params

    def read_poscar(self):
        """Read POSCAR file to get the order, count of elements, and atom coordinates."""
        with open(self.poscar_file, 'r') as file:
            lines = file.readlines()

        scaling_factor = float(lines[1].strip())
        
        # Read lattice vectors
        a = np.array([float(x) for x in lines[2].strip().split()]) * scaling_factor
        b = np.array([float(x) for x in lines[3].strip().split()]) * scaling_factor
        c = np.array([float(x) for x in lines[4].strip().split()]) * scaling_factor
        
        element_symbols = lines[5].strip().split()
        element_counts = [int(x) for x in lines[6].strip().split()]
        
        # Determine if coordinates are Direct or Cartesian
        coord_type = lines[7].strip()[0].upper()
        
        # Read atom coordinates
        coord_start_line = 8
        atom_coords = []
        atom_index = 0
        
        for elem_idx, count in enumerate(element_counts):
            for _ in range(count):
                if atom_index + coord_start_line < len(lines):
                    coords = np.array([float(x) for x in lines[atom_index + coord_start_line].strip().split()[:3]])
                    # Convert to Cartesian if needed
                    if coord_type == 'D' or coord_type == 'C':
                        if coord_type == 'D':
                            # Convert from Direct to Cartesian
                            coords = coords[0]*a + coords[1]*b + coords[2]*c
                        atom_coords.append({
                            'element': element_symbols[elem_idx],
                            'coords': coords
                        })
                atom_index += 1
        
        # Create element list with the correct order
        elements = [elem for elem, count in zip(element_symbols, element_counts) for _ in range(count)]
        
        return elements, atom_coords

    def extract_g0_contribution(self):
        """Extract G=0 contribution to chemical shift from OUTCAR."""
        g0_contrib = None
        
        with open(self.outcar_file, 'r') as file:
            lines = file.readlines()

        for i, line in enumerate(lines):
            if "G=0 CONTRIBUTION TO CHEMICAL SHIFT" in line:
                start_index = i + 5
                if start_index + 2 < len(lines):
                    try:
                        g0_contrib = np.array([
                            [float(x) for x in lines[start_index].split()[1:]],
                            [float(x) for x in lines[start_index + 1].split()[1:]],
                            [float(x) for x in lines[start_index + 2].split()[1:]]
                        ])
                    except (ValueError, IndexError):
                        if self.verbosity > 0:
                            print("Error extracting G=0 contribution.")
                        g0_contrib = None
                break

        return g0_contrib

    def extract_core_susceptibility_and_volume(self):
        """Extract core susceptibility and volume from OUTCAR file."""
        core_suscept = None
        volume = None
        
        with open(self.outcar_file, 'r') as file:
            for line in file:
                if " Core contribution to magnetic susceptibility:" in line:
                    try:
                        match = re.search(r"Core contribution to magnetic susceptibility:\s+([-+]?[0-9]*\.?[0-9]+)\s+10\^-6", line)
                        if match:
                            core_suscept = float(match.group(1))
                    except ValueError:
                        core_suscept = None

                if "volume of cell :" in line:
                    try:
                        volume = float(line.split()[-1])
                    except ValueError:
                        volume = None
                        
                if core_suscept is not None and volume is not None:
                    break
                        
        return core_suscept, volume

    def extract_vasp_csa_reference(self):
        """
        Extract VASP's own reported sigma_iso, span, and skew from the CSA tensor block
        in the OUTCAR file. Specifically reads the 'absolute, valence and core' section
        under 'INCLUDING G=0 CONTRIBUTION'.

        The OUTCAR block has this structure:
          CSA tensor (J. Mason, ...)
          ...
            EXCLUDING G=0 CONTRIBUTION    INCLUDING G=0 CONTRIBUTION
          ATOM    ISO_SHIFT  SPAN  SKEW    ISO_SHIFT  SPAN  SKEW
          (absolute, valence only)
             1   ...
          (absolute, valence and core)
             1   excl_iso  excl_span  excl_skew   incl_iso  incl_span  incl_skew
             2   ...

        Returns a dict: {ion_index (1-based): {'sigma_iso': ..., 'span': ..., 'skew': ...}}
        """
        reference = {}

        with open(self.outcar_file, 'r') as f:
            lines = f.readlines()

        # Locate the CSA tensor block
        csa_start = None
        for i, line in enumerate(lines):
            if 'CSA tensor' in line:
                csa_start = i
                break

        if csa_start is None:
            if self.verbosity > 0:
                print("Warning: 'CSA tensor' block not found in OUTCAR. Sanity check skipped.")
            return {}

        # Within the block, find '(absolute, valence and core)'
        core_section_start = None
        for i in range(csa_start, len(lines)):
            if 'absolute' in lines[i].lower() and 'core' in lines[i].lower() and 'valence' in lines[i].lower():
                core_section_start = i + 1
                if self.verbosity > 1:
                    print(f"  Found core section at line {i}: {lines[i].rstrip()}")
                break

        if core_section_start is None:
            if self.verbosity > 0:
                print("Warning: '(absolute, valence and core)' sub-block not found. Sanity check skipped.")
            return {}

        # Read atom lines until the closing dashes
        for line in lines[core_section_start:]:
            stripped = line.strip()
            if stripped.startswith('---'):
                break
            if not stripped:
                continue
            parts = stripped.split()
            # Expected format: ion  excl_iso  excl_span  excl_skew  incl_iso  incl_span  incl_skew
            if len(parts) == 7:
                try:
                    ion_idx   = int(parts[0])
                    incl_iso  = float(parts[4])
                    incl_span = float(parts[5])
                    incl_skew = float(parts[6])
                    reference[ion_idx] = {
                        'sigma_iso': incl_iso,
                        'span':      incl_span,
                        'skew':      incl_skew,
                    }
                except (ValueError, IndexError):
                    if self.verbosity > 1:
                        print(f"Skipping CSA reference line: {line.rstrip()}")

        return reference

    def run_sanity_check(self, tolerance=0.0001):
        """
        Compare the computed sigma_iso, span, and skew against VASP's own reported
        values from the 'CSA tensor / absolute, valence and core / including G=0' block.

        Values should be IDENTICAL to 4 decimal places. Any discrepancy flags a
        bug in the tensor assembly, conversion factor, or format parsing.

        Parameters
        ----------
        tolerance : float
            Maximum allowed absolute difference (ppm). Default is 0.0001 (4 decimals).

        Returns
        -------
        bool
            True if all atoms match identically, False otherwise.
        """
        if not self.results:
            print("Sanity check: no computed results found. Run calculate_chemical_shifts() first.")
            return False

        reference = self.extract_vasp_csa_reference()
        if not reference:
            print("Sanity check: could not extract VASP reference values.")
            return False

        print("\n" + "=" * 90)
        print("SANITY CHECK: Computed vs. VASP-reported (absolute, valence+core, incl. G=0)")
        print("  Note: VASP ISO_SHIFT sign is reversed; reference values shown here as -ISO_SHIFT.")
        print("=" * 90)
        print(f"{'Ion':>4}  {'Elem':>4}  "
              f"{'σ_iso calc':>12} {'σ_iso VASP':>12} {'Δσ_iso':>9}  "
              f"{'Span calc':>10} {'Span VASP':>10} {'ΔSpan':>8}  "
              f"{'Skew calc':>10} {'Skew VASP':>10} {'ΔSkew':>8}  {'Status':>8}")
        print("-" * 90)

        all_pass = True
        for result in self.results:
            ion = result['ion']
            ref = reference.get(ion)
            if ref is None:
                print(f"{ion:>4}  {'?':>4}  -- no VASP reference entry --")
                continue

            d_iso  = abs(result['sigma_iso'] - (-ref['sigma_iso']))
            d_span = abs(result['span']      - ref['span'])
            d_skew = abs(result['skew']      - ref['skew'])

            ok = (d_iso <= tolerance) and (d_span <= tolerance) and (d_skew <= tolerance)
            status = "OK" if ok else "MISMATCH"
            if not ok:
                all_pass = False

            print(f"{ion:>4}  {result['element']:>4}  "
                  f"{result['sigma_iso']:>12.4f} {-ref['sigma_iso']:>12.4f} {d_iso:>9.4f}  "
                  f"{result['span']:>10.4f} {ref['span']:>10.4f} {d_span:>8.4f}  "
                  f"{result['skew']:>10.4f} {ref['skew']:>10.4f} {d_skew:>8.4f}  {status:>8}")

        print("=" * 90)
        if all_pass:
            print(f"✓ All atoms are identical to 4 decimal places.")
        else:
            print(f"✗ Some atoms differ — the computed values do not reproduce VASP's output exactly.")
            print(f"  Likely causes: wrong tensor format (3×3 vs Voigt), wrong conversion factor, or sign issue.")
        print("=" * 90 + "\n")

        return all_pass

    def read_outcar(self):
        """
        Read and extract all data from OUTCAR file.
        """
        symm_tensors = self.extract_symmetrized_tensors()
        g0_contrib = self.extract_g0_contribution()
        core_suscept, volume = self.extract_core_susceptibility_and_volume()
        core_nmr_properties = self.extract_nmr_values()
        quadrupolar_params = self.extract_quadrupolar_parameters()
        
        return symm_tensors, g0_contrib, core_suscept, volume, core_nmr_properties, quadrupolar_params

    def calculate_chemical_shifts(self):
        """Calculate magnetic shielding properties for all atoms."""
        element_list, atom_coordinates = self.read_poscar()
        symm_tensors, g0_contrib, core_suscept, volume, core_nmr, quadrupolar_params = self.read_outcar()
        
        if self.verbosity > 0:
            print(f"Processing NMR parameters for {len(element_list)} atoms...")

        # Cap tensors to actual atom count — OUTCAR may contain multiple tensor blocks
        if len(symm_tensors) != len(element_list):
            if self.verbosity > 0:
                print(f"  Note: found {len(symm_tensors)} tensors in OUTCAR but {len(element_list)} atoms in POSCAR. "
                      f"Using first {len(element_list)} tensors.")
            symm_tensors = symm_tensors[:len(element_list)]
        
        # Organize quadrupolar parameters by ion index for quick lookup
        quad_lookup = {param['ion']: param for param in quadrupolar_params}
        
        # Check if we have all required data
        if not symm_tensors:
            print("Error: No symmetrized tensors found.")
            return []
            
        if g0_contrib is None:
            print("Warning: No G=0 contribution found, using zero matrix.")
            g0_contrib = np.zeros((3, 3))
            
        if core_suscept is None:
            print("Warning: No core susceptibility found, using zero.")
            core_suscept = 0.0
            
        if volume is None:
            print("Warning: No cell volume found, using 1.0.")
            volume = 1.0
            
        # Store the coordinates
        self.atom_coordinates = atom_coordinates
        
        # Calculate conversion factor for core susceptibility
        conversion_factor = 1 / (volume * 0.6022142 * 3 / (8 * np.pi))
        
        results = []
        for i, tensor in enumerate(symm_tensors):
            element = element_list[i] if i < len(element_list) else "Unknown"
            coords = atom_coordinates[i]['coords'] if i < len(atom_coordinates) else np.array([0.0, 0.0, 0.0])
            
            # Calculate total tensor
            total_tensor = tensor + g0_contrib
            
            # Add core NMR contribution if available
            core_nmr_contrib = core_nmr.get(element, np.zeros((3, 3)))
            if core_nmr_contrib.shape == (3, 3):
                total_tensor += core_nmr_contrib
            
            # Add core susceptibility contribution
            core_suscept_contrib = np.eye(3) * (core_suscept * conversion_factor)
            total_tensor += core_suscept_contrib
            
            # Apply sign convention for magnetic shielding
            total_tensor *= -1
            
            # Ensure symmetry of the tensor
            total_tensor = (total_tensor + total_tensor.T) / 2
            
            # Diagonalize the tensor
            eigvals, eigvecs = np.linalg.eigh(total_tensor)
            
            # Sort eigenvalues by convention: sigma_11 ≤ sigma_22 ≤ sigma_33
            sigma_11, sigma_22, sigma_33 = sorted(eigvals)
            
            # Calculate tensor properties
            sigma_iso = (sigma_11 + sigma_22 + sigma_33) / 3
            span = sigma_33 - sigma_11
            skew = 3 * (sigma_iso - sigma_22) / span if span != 0 else 0
            
            # Calculate CSA and asymmetry parameter (Haeberlen convention)
            # sigma_zz is the component furthest from sigma_iso
            if abs(sigma_33 - sigma_iso) >= abs(sigma_11 - sigma_iso):
                sigma_csa = sigma_33 - sigma_iso   # Δσ > 0
                eta_csa = (sigma_22 - sigma_11) / sigma_csa if sigma_csa != 0 else 0
            else:
                sigma_csa = sigma_11 - sigma_iso   # Δσ < 0
                eta_csa = (sigma_22 - sigma_33) / sigma_csa if sigma_csa != 0 else 0
            
            # Get quadrupolar parameters if available
            c_q = quad_lookup.get(i+1, {}).get('Cq', None)
            eta_q = quad_lookup.get(i+1, {}).get('eta_Q', None)
            
            # Create result dictionary
            result = {
                "ion": i + 1,
                "element": element,
                "coords": coords,
                "non_diag_tensor": total_tensor,
                "diag_tensor": np.diag([sigma_11, sigma_22, sigma_33]),
                "principal_components": (sigma_11, sigma_22, sigma_33),
                "sigma_iso": sigma_iso,
                "sigma_csa": sigma_csa,
                "eta_csa": eta_csa,
                "span": span,
                "skew": skew
            }
            
            # Add quadrupolar parameters if available
            if c_q is not None:
                result["C_Q"] = c_q
            if eta_q is not None:
                result["eta_Q"] = eta_q
                
            results.append(result)
            
        self.results = results
        return results

    def write_detailed_report(self, output_file):
        """Write the final results to a detailed report file."""
        with open(output_file, 'w') as file:
            file.write(f"####################################< There are {len(self.results)} atoms in this system >#############################################\n")
            file.write("#\n")
            file.write("# Convention: Principal components ordered as σ_11 ≤ σ_22 ≤ σ_33 (shielding convention,\n")
            file.write("#   σ_11 = least shielded, σ_33 = most shielded).\n")
            file.write("#   σ_csa = σ_zz - σ_iso (Haeberlen reduced anisotropy / CSA); σ_zz is the principal\n")
            file.write("#   component furthest from σ_iso. σ_csa > 0 if σ_33 dominates, σ_csa < 0 if σ_11 dominates.\n")
            file.write("#   η_csa = (σ_yy - σ_xx) / σ_csa ∈ [0, 1] (Haeberlen asymmetry parameter).\n")
            file.write("#   Span Ω = σ_33 - σ_11, Skew κ = 3(σ_iso - σ_22) / Ω (Herzfeld-Berger convention).\n")
            file.write("#   Note: the full anisotropy Δσ = (3/2) σ_csa if needed.\n")
            file.write("#\n")
            for result in self.results:
                file.write("\n")
                file.write(f"Ion Index.: {result['ion']}\n")
                file.write(f"Ion Symbol: {result['element']}\n")
                file.write(f"Coordinates: [{result['coords'][0]:.6f}, {result['coords'][1]:.6f}, {result['coords'][2]:.6f}]\n")
                file.write("\n")
                file.write("Non-diagonalized Total Tensor:\n")
                for row in result['non_diag_tensor']:
                    file.write(f"{row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}\n")
                file.write("\n")
                file.write("Diagonalized Total Tensor:\n")
                for i in range(3):
                    row = [0, 0, 0]
                    row[i] = result['diag_tensor'][i, i]
                    file.write(f"{row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}\n")
                file.write("\n")
                file.write(f"Principal Components [ppm]:\n")
                file.write(f"  sigma_11 = {result['principal_components'][0]:12.6f}\n")
                file.write(f"  sigma_22 = {result['principal_components'][1]:12.6f}\n")
                file.write(f"  sigma_33 = {result['principal_components'][2]:12.6f}\n")
                file.write("\n")
                file.write(f"Magnetic Shielding Parameters:\n")
                file.write(f"  sigma_iso [ppm] = {result['sigma_iso']:12.6f}\n")
                file.write(f"  sigma_csa [ppm] = {result['sigma_csa']:12.6f}\n")
                file.write(f"  eta_csa        = {result['eta_csa']:12.6f}\n")
                file.write(f"  span [ppm]     = {result['span']:12.6f}\n")
                file.write(f"  skew           = {result['skew']:12.6f}\n")
                file.write("\n")
                
                # Add quadrupolar parameters if available
                if 'C_Q' in result or 'eta_Q' in result:
                    file.write("Quadrupolar Parameters:\n")
                    if 'C_Q' in result:
                        file.write(f"  C_Q [MHz]      = {result['C_Q']:12.6f}\n")
                    if 'eta_Q' in result:
                        file.write(f"  eta_Q          = {result['eta_Q']:12.6f}\n")
                    file.write("\n")
                    
                file.write("#######################################################################################################################")
                file.write("\n")

    def write_organized_report(self, output_file):
        """Write the results to an organized, tabular report file."""
        with open(output_file, 'w') as file:
            # Write header
            file.write("# VASP NMR Parameters - Organized Report\n\n")
            file.write(f"Total atoms: {len(self.results)}\n\n")
            file.write("## Convention\n\n")
            file.write("- Principal components: σ_11 ≤ σ_22 ≤ σ_33 (shielding convention; σ_11 = least shielded, σ_33 = most shielded)\n")
            file.write("- σ_csa = σ_zz − σ_iso (Haeberlen reduced anisotropy / CSA); σ_zz is the component furthest from σ_iso; σ_csa > 0 if σ_33 dominates, σ_csa < 0 if σ_11 dominates\n")
            file.write("- η_csa ∈ [0, 1]: Haeberlen asymmetry parameter, η = (σ_yy − σ_xx) / σ_csa\n")
            file.write("- Span Ω = σ_33 − σ_11, Skew κ = 3(σ_iso − σ_22) / Ω (Herzfeld-Berger convention)\n")
            file.write("- Note: full anisotropy Δσ = (3/2) σ_csa if needed\n\n")
            
            # Write magnetic shielding parameters table
            file.write("## Magnetic Shielding Parameters\n\n")
            file.write("| Ion | Element | X | Y | Z | σ_iso | σ_csa | η_csa | Span | Skew | σ_11 | σ_22 | σ_33 |\n")
            file.write("|-----|---------|---|---|---|-------|-------|-------|------|------|------|------|------|\n")
            
            for result in self.results:
                x, y, z = result['coords']
                file.write(f"| {result['ion']:3d} | {result['element']:7s} | ")
                file.write(f"{x:6.3f} | {y:6.3f} | {z:6.3f} | ")
                file.write(f"{result['sigma_iso']:7.3f} | {result['sigma_csa']:7.3f} | {result['eta_csa']:7.3f} | ")
                file.write(f"{result['span']:6.3f} | {result['skew']:6.3f} | ")
                file.write(f"{result['principal_components'][0]:6.3f} | {result['principal_components'][1]:6.3f} | {result['principal_components'][2]:6.3f} |\n")
            
            # Write quadrupolar parameters table if available
            has_quadrupolar = any('C_Q' in result for result in self.results)
            if has_quadrupolar:
                file.write("\n## Quadrupolar Parameters\n\n")
                file.write("| Ion | Element | C_Q [MHz] | η_Q |\n")
                file.write("|-----|---------|-----------|-----|\n")
                
                for result in self.results:
                    if 'C_Q' in result:
                        file.write(f"| {result['ion']:3d} | {result['element']:7s} | ")
                        file.write(f"{result.get('C_Q', 'N/A'):9.3f} | {result.get('eta_Q', 'N/A'):5.3f} |\n")

    def write_csv_report(self, output_file):
        """Write the results to a CSV file."""
        with open(output_file, 'w', newline='') as csvfile:
            # Determine all possible fields
            all_fields = set()
            for result in self.results:
                all_fields.update(result.keys())
            
            # Remove fields that need special handling
            special_fields = ['non_diag_tensor', 'diag_tensor', 'principal_components', 'coords']
            for field in special_fields:
                if field in all_fields:
                    all_fields.remove(field)
            
            # Create fieldnames list with special fields added back appropriately
            fieldnames = ['ion', 'element', 'x', 'y', 'z']
            fieldnames.extend([f for f in sorted(all_fields) if f not in ['ion', 'element']])
            fieldnames.extend(['sigma_11', 'sigma_22', 'sigma_33'])
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in self.results:
                row = {k: result[k] for k in all_fields if k in result}
                
                # Handle coordinates
                if 'coords' in result:
                    row['x'] = result['coords'][0]
                    row['y'] = result['coords'][1]
                    row['z'] = result['coords'][2]
                
                # Handle principal components
                if 'principal_components' in result:
                    row['sigma_11'] = result['principal_components'][0]
                    row['sigma_22'] = result['principal_components'][1]
                    row['sigma_33'] = result['principal_components'][2]
                
                writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description='VASP NMR Parameter Calculator')
    parser.add_argument('--poscar', required=True, help='Path to the POSCAR file')
    parser.add_argument('--outcar', required=True, help='Path to the OUTCAR file')
    parser.add_argument('--report', required=False, default='nmr_report.txt', help='Output file for detailed report')
    parser.add_argument('--organized_report', required=False, default='nmr_organized_report.md', help='Output file for organized report')
    parser.add_argument('--csv', required=False, default='nmr_data.csv', help='Output file for CSV data')
    parser.add_argument('--verbosity', type=int, choices=[0, 1, 2], default=1, help='Verbosity level: 0=silent, 1=normal, 2=verbose')
    parser.add_argument('--sanity_check', action='store_true', help='Compare computed parameters against VASP-reported CSA tensor values')
    parser.add_argument('--tolerance', type=float, default=0.0001, help='Tolerance in ppm for sanity check comparison (default: 0.0001, i.e. 4-decimal agreement)')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.isfile(args.poscar):
        print(f"Error: POSCAR file '{args.poscar}' not found.")
        sys.exit(1)
    
    if not os.path.isfile(args.outcar):
        print(f"Error: OUTCAR file '{args.outcar}' not found.")
        sys.exit(1)
    
    # Initialize calculator and process data
    calculator = NMRCalculator(args.outcar, args.poscar, args.verbosity)
    calculator.calculate_chemical_shifts()

    # Optional sanity check against VASP's own CSA tensor block
    if args.sanity_check:
        calculator.run_sanity_check(tolerance=args.tolerance)
    
    # Write output files
    calculator.write_detailed_report(args.report)
    calculator.write_organized_report(args.organized_report)
    calculator.write_csv_report(args.csv)
    
    if args.verbosity > 0:
        print(f"\nProcessing complete!")
        print(f"Detailed report written to: {args.report}")
        print(f"Organized report written to: {args.organized_report}")
        print(f"CSV data written to: {args.csv}")

if __name__ == "__main__":
    main()
