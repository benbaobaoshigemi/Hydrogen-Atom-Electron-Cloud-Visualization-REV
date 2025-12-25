"""
Schrödinger Equation Inversion Method for Extracting Effective Potential
=========================================================================
Given a known wavefunction R(r) and its energy eigenvalue ε, we can INVERT
the radial Schrödinger equation to find the exact effective potential V_eff(r).

The radial Schrödinger equation is:
    [-ℏ²/2m d²/dr² - (ℏ²/m)(1/r)d/dr + ℏ²l(l+1)/(2mr²) + V_eff(r)] R(r) = ε R(r)

Rearranging (in atomic units where ℏ=m=1):
    V_eff(r) = ε - [-1/2 R''(r)/R(r) - (1/r) R'(r)/R(r) + l(l+1)/(2r²)]
    V_eff(r) = ε + (1/2)[R''(r)/R(r)] + (1/r)[R'(r)/R(r)] - l(l+1)/(2r²)

This method is EXACT because it extracts the potential that makes R(r) a solution.
"""

import json
import numpy as np
from math import factorial as math_factorial
import os

# Configuration
R_MAX = 20.0
N_POINTS = 1000  # High resolution for derivatives
OUTPUT_FILE = "veff_inversion_tables.json"

def factorial(n):
    if n < 0:
        return 0
    return math_factorial(int(n))

def slater_radial(r, basis):
    """Compute R(r) from STO basis."""
    val = 0.0
    for term in basis:
        n = term['nStar']
        zeta = term['zeta']
        c = term['coeff']
        N = np.sqrt((2 * zeta) ** (2 * n + 1) / factorial(2 * n))
        val += c * N * (r ** (n - 1)) * np.exp(-zeta * r)
    return val

def slater_radial_derivative(r, basis):
    """Compute R'(r) = dR/dr analytically."""
    val = 0.0
    for term in basis:
        n = term['nStar']
        zeta = term['zeta']
        c = term['coeff']
        N = np.sqrt((2 * zeta) ** (2 * n + 1) / factorial(2 * n))
        # d/dr [r^(n-1) * exp(-ζr)] = (n-1)r^(n-2) * exp(-ζr) - ζ * r^(n-1) * exp(-ζr)
        #                          = [(n-1)/r - ζ] * r^(n-1) * exp(-ζr)
        if r > 1e-10:
            factor = (n - 1) / r - zeta
        else:
            factor = -zeta  # At r=0, the (n-1)/r term needs care
        val += c * N * factor * (r ** (n - 1)) * np.exp(-zeta * r)
    return val

def slater_radial_second_derivative(r, basis):
    """Compute R''(r) = d²R/dr² analytically."""
    val = 0.0
    for term in basis:
        n = term['nStar']
        zeta = term['zeta']
        c = term['coeff']
        N = np.sqrt((2 * zeta) ** (2 * n + 1) / factorial(2 * n))
        
        # f = r^(n-1) * exp(-ζr)
        # f' = [(n-1)/r - ζ] * f
        # f'' = d/dr[(n-1)/r - ζ] * f + [(n-1)/r - ζ] * f'
        #     = -(n-1)/r² * f + [(n-1)/r - ζ]² * f
        #     = [-(n-1)/r² + ((n-1)/r - ζ)²] * f
        
        if r > 1e-10:
            factor1 = (n - 1) / r - zeta
            factor2 = -(n - 1) / (r * r) + factor1 * factor1
        else:
            factor2 = zeta * zeta  # Limiting behavior
        
        val += c * N * factor2 * (r ** (n - 1)) * np.exp(-zeta * r)
    return val

def extract_effective_potential(r_grid, basis, epsilon, l):
    """
    Extract V_eff(r) by inverting the Schrödinger equation.
    
    V_eff(r) = ε + (1/2)[R''(r)/R(r)] + (1/r)[R'(r)/R(r)] - l(l+1)/(2r²)
    """
    V_eff = np.zeros_like(r_grid)
    
    for i, r in enumerate(r_grid):
        if r < 0.05:  # Skip very small r to avoid numerical issues
            V_eff[i] = np.nan
            continue
        
        R = slater_radial(r, basis)
        R_prime = slater_radial_derivative(r, basis)
        R_double_prime = slater_radial_second_derivative(r, basis)
        
        if abs(R) < 1e-15:  # Avoid division by zero
            V_eff[i] = np.nan
            continue
        
        # Kinetic energy contributions
        kinetic_radial = 0.5 * R_double_prime / R
        kinetic_angular_1 = (1 / r) * R_prime / R
        kinetic_angular_2 = l * (l + 1) / (2 * r * r)
        
        V_eff[i] = epsilon + kinetic_radial + kinetic_angular_1 - kinetic_angular_2
    
    # Interpolate NaN values at small r
    valid_mask = ~np.isnan(V_eff)
    if np.any(valid_mask):
        first_valid_idx = np.argmax(valid_mask)
        if first_valid_idx > 0:
            V_eff[:first_valid_idx] = V_eff[first_valid_idx]
    
    return V_eff

def get_l_from_orbital(orbital_key):
    """Extract angular momentum l from orbital name like '1s', '2p', '3d'."""
    orbital_type = orbital_key[-1].lower()
    l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}
    return l_map.get(orbital_type, 0)

def main():
    print("=" * 70)
    print("Schrödinger Equation Inversion Method")
    print("Extracting Exact V_eff(r) from Koga Wavefunctions")
    print("=" * 70)
    
    db_path = 'koga_basis_database.json'
    if not os.path.exists(db_path):
        print(f"ERROR: {db_path} not found!")
        return
    
    with open(db_path, 'r', encoding='utf-8') as f:
        koga_db = json.load(f)
    
    r_grid = np.linspace(0.01, R_MAX, N_POINTS)
    
    results = {
        "metadata": {
            "method": "Schrödinger Equation Inversion",
            "description": "Exact V_eff(r) extracted from Koga HF wavefunctions",
            "r_min": float(r_grid[0]),
            "r_max": float(r_grid[-1]),
            "n_points": N_POINTS
        },
        "atoms": {}
    }
    
    for symbol, atom_data in koga_db.items():
        if symbol in ['Au_R']:
            continue
        
        Z = atom_data['Z']
        if Z > 54:
            continue
        
        print(f"\nProcessing {symbol} (Z={Z})...")
        
        atom_results = {
            "Z": Z,
            "orbitals": {}
        }
        
        for orbital_key, orbital_data in atom_data['orbitals'].items():
            basis = orbital_data['basis']
            epsilon = orbital_data.get('energy', 0)
            l = get_l_from_orbital(orbital_key)
            
            print(f"  {orbital_key}: ε = {epsilon:.6f} Hartree, l = {l}")
            
            V_eff = extract_effective_potential(r_grid, basis, epsilon, l)
            
            # Validate: at large r, V_eff should approach -Z_eff/r
            # For outer electrons, Z_eff ≈ 1-10
            r_test = 5.0
            idx_test = np.argmin(np.abs(r_grid - r_test))
            if not np.isnan(V_eff[idx_test]):
                Z_eff_apparent = -V_eff[idx_test] * r_test
                print(f"    V_eff(r=5) = {V_eff[idx_test]:.4f}, apparent Z_eff = {Z_eff_apparent:.2f}")
            
            # Downsample for JSON
            step = 10
            atom_results["orbitals"][orbital_key] = {
                "epsilon": epsilon,
                "l": l,
                "r": r_grid[::step].tolist(),
                "V_eff": V_eff[::step].tolist()
            }
        
        results["atoms"][symbol] = atom_results
    
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'=' * 70}")
    print(f"Results saved to {OUTPUT_FILE}")
    print("=" * 70)
    
    # Generate verification plot
    try:
        import matplotlib.pyplot as plt
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Plot 1: Hydrogen 1s (should give -1/r exactly)
        h_data = results["atoms"]["H"]["orbitals"]["1s"]
        r_h = np.array(h_data["r"])
        V_h = np.array(h_data["V_eff"])
        V_exact_h = -1 / r_h  # Hydrogen exact
        
        axes[0, 0].plot(r_h, V_h, 'b-', linewidth=2, label='V_eff (inverted)')
        axes[0, 0].plot(r_h, V_exact_h, 'r--', label='-1/r (exact)')
        axes[0, 0].set_xlim(0, 10)
        axes[0, 0].set_ylim(-2, 0)
        axes[0, 0].set_xlabel('r (a₀)')
        axes[0, 0].set_ylabel('V_eff (Hartree)')
        axes[0, 0].set_title('Hydrogen 1s (Exact Test)')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Helium 1s
        he_data = results["atoms"]["He"]["orbitals"]["1s"]
        r_he = np.array(he_data["r"])
        V_he = np.array(he_data["V_eff"])
        V_nuc_he = -2 / r_he
        
        axes[0, 1].plot(r_he, V_nuc_he, 'b--', label='-Z/r (bare nucleus)')
        axes[0, 1].plot(r_he, V_he, 'r-', linewidth=2, label='V_eff (inverted)')
        axes[0, 1].set_xlim(0, 5)
        axes[0, 1].set_ylim(-5, 0)
        axes[0, 1].set_xlabel('r (a₀)')
        axes[0, 1].set_ylabel('V_eff (Hartree)')
        axes[0, 1].set_title('Helium 1s: Effective Potential')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Carbon orbitals
        if "C" in results["atoms"]:
            c_data = results["atoms"]["C"]
            for orb_key, orb_data in c_data["orbitals"].items():
                r_c = np.array(orb_data["r"])
                V_c = np.array(orb_data["V_eff"])
                axes[1, 0].plot(r_c, V_c, linewidth=2, label=f'{orb_key}')
            
            axes[1, 0].plot(r_c, -6/r_c, 'k--', alpha=0.5, label='-Z/r')
            axes[1, 0].set_xlim(0, 5)
            axes[1, 0].set_ylim(-20, 0)
            axes[1, 0].set_xlabel('r (a₀)')
            axes[1, 0].set_ylabel('V_eff (Hartree)')
            axes[1, 0].set_title('Carbon: Orbital-Dependent V_eff')
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Titanium orbitals (transition metal)
        if "Ti" in results["atoms"]:
            ti_data = results["atoms"]["Ti"]
            for orb_key, orb_data in ti_data["orbitals"].items():
                r_ti = np.array(orb_data["r"])
                V_ti = np.array(orb_data["V_eff"])
                axes[1, 1].plot(r_ti, V_ti, linewidth=2, label=f'{orb_key}')
            
            axes[1, 1].plot(r_ti, -22/r_ti, 'k--', alpha=0.5, label='-Z/r')
            axes[1, 1].set_xlim(0, 3)
            axes[1, 1].set_ylim(-50, 0)
            axes[1, 1].set_xlabel('r (a₀)')
            axes[1, 1].set_ylabel('V_eff (Hartree)')
            axes[1, 1].set_title('Titanium: Orbital-Dependent V_eff')
            axes[1, 1].legend()
            axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('veff_inversion_verification.png', dpi=150)
        print(f"\nVerification plot saved to veff_inversion_verification.png")
        
    except ImportError:
        print("matplotlib not available")

if __name__ == "__main__":
    main()
