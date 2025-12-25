"""
Full Hartree Potential Calculator for Electron Cloud Visualization
=====================================================================
This script computes the effective potential V_eff(r) for all atoms (Z=1 to 54)
using the Koga high-precision Slater-Type Orbital basis sets.

The Hartree potential is computed via spherical integration:
V_H(r) = (1/r) * ∫₀^r ρ(r') 4πr'² dr' + ∫_r^∞ ρ(r') 4πr' dr'

Where ρ(r) is the total electron density summed over ALL occupied orbitals.

Output: hartree_potential_tables.json containing pre-computed V_eff(r) for each atom.
"""

import json
import numpy as np
from math import factorial as math_factorial
from scipy import integrate
import os

# Configuration
R_MAX = 20.0  # Maximum radius in Bohr
N_POINTS = 500  # Number of radial grid points
OUTPUT_FILE = "hartree_potential_tables.json"

def factorial(n):
    """Safe factorial for large n."""
    if n < 0:
        return 0
    return math_factorial(int(n))

def slater_radial_normalized(r, basis):
    """
    Compute the normalized radial wavefunction R(r) for a Slater-type orbital.
    
    The STO radial function is: R(r) = Σᵢ cᵢ * Nᵢ * r^(nᵢ-1) * exp(-ζᵢ * r)
    
    Normalization: ∫₀^∞ R(r)² r² dr = 1
    """
    val = 0.0
    for term in basis:
        n = term['nStar']
        zeta = term['zeta']
        c = term['coeff']
        
        # Normalization factor for STO: N = sqrt((2ζ)^(2n+1) / (2n)!)
        N = np.sqrt((2 * zeta) ** (2 * n + 1) / factorial(2 * n))
        
        val += c * N * (r ** (n - 1)) * np.exp(-zeta * r)
    
    return val

def compute_radial_probability_density(r_grid, basis):
    """
    Compute P(r) = R(r)² * r² (radial probability density).
    ∫₀^∞ P(r) dr = 1 for a normalized orbital.
    """
    R_vals = np.array([slater_radial_normalized(r, basis) for r in r_grid])
    P_vals = (R_vals ** 2) * (r_grid ** 2)
    return R_vals, P_vals

def compute_hartree_potential_for_density(r_grid, rho_total):
    """
    Compute the Hartree potential V_H(r) for a given spherical charge density rho(r).
    
    Using the classical electrostatics formula for spherical charge distribution:
    V_H(r) = (1/r) * ∫₀^r ρ(r') 4πr'² dr' + ∫_r^∞ (ρ(r') 4πr'² / r') dr'
    
    This simplifies to (with 4π factor absorbed into rho being probability density):
    V_H(r) = Q_in(r)/r + integral of outer shell contribution
    
    Returns: V_H(r) array (positive values, representing electron repulsion)
    """
    dr = r_grid[1] - r_grid[0]
    n_points = len(r_grid)
    
    # rho_total here is the radial probability density P(r) = R²r²
    # Q(r) = ∫₀^r P(r') dr' = enclosed probability (electron charge in a.u.)
    
    V_hartree = np.zeros(n_points)
    
    # Precompute cumulative integrals for efficiency
    # Q_in[i] = ∫₀^r[i] P(r') dr' (enclosed charge)
    Q_in = np.cumsum(rho_total) * dr
    
    # For outer integral: ∫_r^∞ P(r')/r' dr'
    # We compute this as reverse cumsum of P(r)/r
    integrand_outer = np.divide(rho_total, r_grid, out=np.zeros_like(rho_total), where=r_grid > 1e-10)
    # Reverse cumulative sum (from infinity to r)
    Q_out_contribution = np.flip(np.cumsum(np.flip(integrand_outer))) * dr
    
    for i in range(n_points):
        r = r_grid[i]
        if r < 1e-10:
            # At r=0, we need to be careful. V_H(0) is finite.
            # Use L'Hopital or set to the outer integral value at 0.
            V_hartree[i] = Q_out_contribution[i] if i < len(Q_out_contribution) else 0
        else:
            # V_H(r) = Q_in(r)/r + ∫_r^∞ P(r')/r' dr'
            V_hartree[i] = Q_in[i] / r + Q_out_contribution[i]
    
    return V_hartree

def parse_ground_state(ground_state_str):
    """
    Parse ground state string like "1s²2s²2p⁶" into orbital occupancies.
    Returns dict: {'1s': 2, '2s': 2, '2p': 6, ...}
    """
    import re
    
    # Superscript to digit mapping
    superscripts = {'⁰': 0, '¹': 1, '²': 2, '³': 3, '⁴': 4, 
                    '⁵': 5, '⁶': 6, '⁷': 7, '⁸': 8, '⁹': 9, '¹⁰': 10}
    
    occupancies = {}
    
    # Pattern: orbital name followed by superscript number(s)
    # e.g., "1s²" or "3d¹⁰"
    pattern = r'(\d[spdfg])([⁰¹²³⁴⁵⁶⁷⁸⁹]+)'
    
    for match in re.finditer(pattern, ground_state_str):
        orbital = match.group(1)
        superscript_str = match.group(2)
        
        # Convert superscript string to number
        num = 0
        i = 0
        while i < len(superscript_str):
            # Check for two-character superscripts like ¹⁰
            if i + 1 < len(superscript_str) and superscript_str[i:i+2] in superscripts:
                num = num * 10 + superscripts[superscript_str[i:i+2]]
                i += 2
            elif superscript_str[i] in superscripts:
                num = num * 10 + superscripts[superscript_str[i]]
                i += 1
            else:
                i += 1
        
        occupancies[orbital] = num
    
    return occupancies

def compute_total_electron_density(r_grid, atom_data):
    """
    Compute total electron density by summing over all occupied orbitals.
    
    rho_total(r) = Σᵢ nᵢ * |Rᵢ(r)|² * r²
    
    where nᵢ is the occupancy of orbital i.
    """
    orbitals = atom_data['orbitals']
    ground_state = atom_data.get('groundState', '')
    
    # Parse ground state to get occupancies
    occupancies = parse_ground_state(ground_state)
    
    # If no ground state info, assume all available orbitals are filled
    # according to standard filling order
    if not occupancies:
        # Fallback: use aufbau principle
        Z = atom_data['Z']
        # Standard filling order
        filling_order = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p']
        max_electrons = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
        
        remaining = Z
        for orb in filling_order:
            if orb in orbitals and remaining > 0:
                orbital_type = orb[-1]
                max_e = max_electrons[orbital_type]
                occupancies[orb] = min(remaining, max_e)
                remaining -= occupancies[orb]
    
    # Compute total density
    rho_total = np.zeros_like(r_grid)
    
    for orbital_key, orbital_data in orbitals.items():
        if orbital_key in occupancies:
            n_electrons = occupancies[orbital_key]
            # Extract basis array from orbital data (structure is {energy, basis})
            basis = orbital_data['basis'] if isinstance(orbital_data, dict) else orbital_data
            R_vals, P_vals = compute_radial_probability_density(r_grid, basis)
            # Each electron in this orbital contributes P(r) to the density
            rho_total += n_electrons * P_vals
    
    return rho_total

def compute_effective_potential(r_grid, Z, rho_total):
    """
    Compute the effective potential:
    V_eff(r) = V_nuclear(r) + V_Hartree(r)
             = -Z/r + V_H(r)
    
    The V_Hartree is positive (repulsion from other electrons).
    """
    V_nuclear = -Z / np.maximum(r_grid, 1e-10)
    V_hartree = compute_hartree_potential_for_density(r_grid, rho_total)
    V_eff = V_nuclear + V_hartree
    
    return V_eff, V_nuclear, V_hartree

def main():
    print("=" * 70)
    print("Full Hartree Potential Calculator")
    print("Using Koga High-Precision Slater-Type Orbital Basis Sets")
    print("=" * 70)
    
    # Load Koga database
    db_path = 'koga_basis_database.json'
    if not os.path.exists(db_path):
        print(f"ERROR: {db_path} not found!")
        return
    
    with open(db_path, 'r', encoding='utf-8') as f:
        koga_db = json.load(f)
    
    # Create radial grid
    r_grid = np.linspace(0.01, R_MAX, N_POINTS)
    
    # Results storage
    results = {
        "metadata": {
            "r_min": float(r_grid[0]),
            "r_max": float(r_grid[-1]),
            "n_points": N_POINTS,
            "method": "Hartree Potential from Koga STO Basis (1999/2000)"
        },
        "atoms": {}
    }
    
    # Process each atom (Z=1 to 54, excluding special entries like Au_R)
    for symbol, atom_data in koga_db.items():
        if symbol in ['Au_R']:  # Skip relativistic variants
            continue
        
        Z = atom_data['Z']
        if Z > 54:  # Only process up to Xe
            continue
        
        print(f"\nProcessing {symbol} (Z={Z})...")
        
        # Compute total electron density
        rho_total = compute_total_electron_density(r_grid, atom_data)
        
        # Verify normalization (should equal Z for neutral atom)
        dr = r_grid[1] - r_grid[0]
        total_electrons = np.sum(rho_total) * dr
        print(f"  Total electrons (integral check): {total_electrons:.4f} (expected: {Z})")
        
        # Compute effective potential
        V_eff, V_nuclear, V_hartree = compute_effective_potential(r_grid, Z, rho_total)
        
        # Get orbital energy from database for validation
        E_orb = {}
        if 'orbitals' in atom_data:
            for orb_key in atom_data['orbitals']:
                if orb_key in atom_data.get('energies', {}):
                    E_orb[orb_key] = atom_data['energies'][orb_key]
        
        # Store results (downsample for JSON efficiency)
        # Store every 5th point to reduce file size
        step = 5
        results["atoms"][symbol] = {
            "Z": Z,
            "r": r_grid[::step].tolist(),
            "V_eff": V_eff[::step].tolist(),
            "V_nuc": V_nuclear[::step].tolist(),
            "V_hartree": V_hartree[::step].tolist(),
            "E_orb": E_orb
        }
        
        # Print some diagnostic info
        print(f"  V_eff at r=1 a₀: {V_eff[np.argmin(np.abs(r_grid - 1.0))]:.4f} Hartree")
        print(f"  V_hartree at r=1 a₀: {V_hartree[np.argmin(np.abs(r_grid - 1.0))]:.4f} Hartree")
    
    # Save results
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'=' * 70}")
    print(f"Results saved to {OUTPUT_FILE}")
    print(f"Processed {len(results['atoms'])} atoms")
    print("=" * 70)
    
    # Generate verification plot for Helium
    try:
        import matplotlib.pyplot as plt
        
        he_data = results["atoms"]["He"]
        r = np.array(he_data["r"])
        V_eff = np.array(he_data["V_eff"])
        V_nuc = np.array(he_data["V_nuc"])
        V_h = np.array(he_data["V_hartree"])
        
        plt.figure(figsize=(12, 8))
        
        plt.subplot(2, 1, 1)
        plt.plot(r, V_nuc, 'b--', label=f'Bare Nucleus (-Z/r), Z={he_data["Z"]}')
        plt.plot(r, V_h, 'orange', linewidth=2, label='Hartree Shielding V_H(r)')
        plt.plot(r, V_eff, 'r-', linewidth=2, label='Effective Potential V_eff(r)')
        plt.axhline(y=0, color='gray', linestyle=':')
        plt.xlim(0, 10)
        plt.ylim(-5, 2)
        plt.xlabel('r (a₀)')
        plt.ylabel('Potential (Hartree)')
        plt.title('Helium: Effective Potential with Hartree Shielding')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Also plot for a heavier atom (e.g., Ti)
        if "Ti" in results["atoms"]:
            ti_data = results["atoms"]["Ti"]
            r_ti = np.array(ti_data["r"])
            V_eff_ti = np.array(ti_data["V_eff"])
            V_nuc_ti = np.array(ti_data["V_nuc"])
            V_h_ti = np.array(ti_data["V_hartree"])
            
            plt.subplot(2, 1, 2)
            plt.plot(r_ti, V_nuc_ti, 'b--', label=f'Bare Nucleus (-Z/r), Z={ti_data["Z"]}')
            plt.plot(r_ti, V_h_ti, 'orange', linewidth=2, label='Hartree Shielding V_H(r)')
            plt.plot(r_ti, V_eff_ti, 'r-', linewidth=2, label='Effective Potential V_eff(r)')
            plt.axhline(y=0, color='gray', linestyle=':')
            plt.xlim(0, 5)
            plt.ylim(-50, 25)
            plt.xlabel('r (a₀)')
            plt.ylabel('Potential (Hartree)')
            plt.title('Titanium: Effective Potential with Hartree Shielding')
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('hartree_verification.png', dpi=150)
        print(f"\nVerification plot saved to hartree_verification.png")
        
    except ImportError:
        print("matplotlib not available, skipping verification plot")

if __name__ == "__main__":
    main()
