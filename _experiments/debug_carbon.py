# -*- coding: utf-8 -*-
"""
CARBON (Z=6) DETAILED ENERGY ANALYSIS
=====================================

Debug script to understand why Carbon has ~1.5% error.

Reference values from Koga:
- E_tot = -37.68861896 Ha
- ε_1s = -11.3255187 Ha
- ε_2s = -0.7056273 Ha  
- ε_2p = -0.4333405 Ha

Carbon: 1s² 2s² 2p² (³P)
"""

import numpy as np
from pathlib import Path
import json
import re

# Import functions from validate_full_hf
import sys
sys.path.insert(0, str(Path(__file__).parent))
from validate_full_hf import (
    create_log_grid, 
    compute_orbital_radial, 
    compute_Yk,
    compute_kinetic_density,
    compute_nuclear_density,
    compute_J_integral,
    compute_K_integral,
    parse_slater_basis,
    Orbital,
    STOTerm,
    EXCHANGE_COEFFS,
)


def main():
    print("=" * 70)
    print("CARBON (Z=6) DETAILED ENERGY ANALYSIS")
    print("=" * 70)
    
    # Load data
    basis_file = Path(__file__).parent.parent / "slater_basis.js"
    data = parse_slater_basis(basis_file)
    C_data = data['C']
    
    r = create_log_grid(r_min=1e-6, r_max=50.0, N=2000)
    Z = 6
    
    # Build orbitals
    orbitals = {}
    for orb_name in ['1s', '2s', '2p']:
        terms = [
            STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff'])
            for t in C_data['orbitals'][orb_name]
        ]
        l = 0 if 's' in orb_name else 1
        orb = Orbital(name=orb_name, l=l, terms=terms)
        R = compute_orbital_radial(orb, r)
        orbitals[orb_name] = (orb, R, l)
        
        # Check normalization
        norm = np.trapz(R ** 2 * r ** 2, r)
        print(f"{orb_name}: norm = {norm:.8f}")
    
    # Occupation numbers
    n = {'1s': 2, '2s': 2, '2p': 2}
    
    print("\n" + "=" * 70)
    print("ONE-ELECTRON INTEGRALS")
    print("=" * 70)
    
    E_1e = 0.0
    for orb_name in ['1s', '2s', '2p']:
        orb, R, l = orbitals[orb_name]
        T = np.trapz(compute_kinetic_density(orb, R, r) * r ** 2, r)
        V = np.trapz(compute_nuclear_density(R, r, Z) * r ** 2, r)
        h = T + V
        E_1e += n[orb_name] * h
        print(f"{orb_name}: T = {T:10.6f}, V = {V:10.6f}, h = {h:10.6f}, n*h = {n[orb_name]*h:10.6f}")
    
    print(f"\nE_1e (total) = {E_1e:.8f} Ha")
    
    print("\n" + "=" * 70)
    print("TWO-ELECTRON INTEGRALS (J and K matrix)")
    print("=" * 70)
    
    J_mat = {}
    K_mat = {}
    
    for name_i in ['1s', '2s', '2p']:
        orb_i, R_i, l_i = orbitals[name_i]
        for name_j in ['1s', '2s', '2p']:
            orb_j, R_j, l_j = orbitals[name_j]
            
            J_ij = compute_J_integral(R_i, R_j, r)
            K_ij = compute_K_integral(R_i, R_j, l_i, l_j, r)
            
            J_mat[(name_i, name_j)] = J_ij
            K_mat[(name_i, name_j)] = K_ij
    
    print("\nJ matrix:")
    print(f"{'':>6} {'1s':>10} {'2s':>10} {'2p':>10}")
    for name_i in ['1s', '2s', '2p']:
        row = f"{name_i:>6}"
        for name_j in ['1s', '2s', '2p']:
            row += f" {J_mat[(name_i, name_j)]:>10.6f}"
        print(row)
    
    print("\nK matrix:")
    print(f"{'':>6} {'1s':>10} {'2s':>10} {'2p':>10}")
    for name_i in ['1s', '2s', '2p']:
        row = f"{name_i:>6}"
        for name_j in ['1s', '2s', '2p']:
            row += f" {K_mat[(name_i, name_j)]:>10.6f}"
        print(row)
    
    print("\n" + "=" * 70)
    print("E_2e CALCULATION METHODS")
    print("=" * 70)
    
    # Method 1: Standard closed-shell formula (treating all as closed)
    E_2e_closed = 0.0
    for name_i in ['1s', '2s', '2p']:
        for name_j in ['1s', '2s', '2p']:
            E_2e_closed += 0.5 * n[name_i] * n[name_j] * (J_mat[(name_i, name_j)] - 0.5 * K_mat[(name_i, name_j)])
    
    print(f"\nMethod 1 (all closed-shell formula):")
    print(f"  E_2e = {E_2e_closed:.8f} Ha")
    print(f"  E_tot = {E_1e + E_2e_closed:.8f} Ha")
    print(f"  E_ref = {C_data['E_tot']:.8f} Ha")
    print(f"  Error = {100 * (E_1e + E_2e_closed - C_data['E_tot']) / abs(C_data['E_tot']):.4f}%")
    
    # Method 2: Roothaan open-shell for 2p
    # For p², the 2p-2p self-interaction uses modified coefficients
    # Standard: 0.5 * n * n * (J - 0.5K) = 0.5 * 2 * 2 * (J - 0.5K) = 2 * (J - 0.5K)
    # Roothaan: pairs * f * (2A*J - B*K) = 1 * (1/3) * (J - 1.5K) for p² ³P
    
    f_2p = 1/3  # from advisor
    A_2p = 0.5  # 2A = 1
    B_2p = 1.5  # from advisor
    
    E_2e_roothaan = 0.0
    for name_i in ['1s', '2s', '2p']:
        for name_j in ['1s', '2s', '2p']:
            if name_i == '2p' and name_j == '2p':
                # Open-shell self-interaction
                # n_open = 2, pairs = 1
                pairs = 1
                E_2e_roothaan += 0.5 * pairs * f_2p * (2 * A_2p * J_mat[(name_i, name_j)] - B_2p * K_mat[(name_i, name_j)])
            else:
                # Standard interaction
                E_2e_roothaan += 0.5 * n[name_i] * n[name_j] * (J_mat[(name_i, name_j)] - 0.5 * K_mat[(name_i, name_j)])
    
    print(f"\nMethod 2 (Roothaan 2p-2p only):")
    print(f"  E_2e = {E_2e_roothaan:.8f} Ha")
    print(f"  E_tot = {E_1e + E_2e_roothaan:.8f} Ha")
    print(f"  Error = {100 * (E_1e + E_2e_roothaan - C_data['E_tot']) / abs(C_data['E_tot']):.4f}%")
    
    # Method 3: Try different f values
    print(f"\nMethod 3 (scanning f values):")
    for test_f in [0.0, 0.1, 0.2, 0.25, 1/3, 0.4, 0.5, 0.6, 2/3, 0.75, 1.0]:
        E_2e_test = 0.0
        for name_i in ['1s', '2s', '2p']:
            for name_j in ['1s', '2s', '2p']:
                if name_i == '2p' and name_j == '2p':
                    pairs = 1
                    E_2e_test += 0.5 * pairs * test_f * (J_mat[(name_i, name_j)] - 1.5 * K_mat[(name_i, name_j)])
                else:
                    E_2e_test += 0.5 * n[name_i] * n[name_j] * (J_mat[(name_i, name_j)] - 0.5 * K_mat[(name_i, name_j)])
        
        E_tot_test = E_1e + E_2e_test
        error = 100 * (E_tot_test - C_data['E_tot']) / abs(C_data['E_tot'])
        print(f"  f = {test_f:.4f}: E_tot = {E_tot_test:.6f}, error = {error:+.4f}%")
    
    # Method 4: What if we DON'T use Roothaan at all for 2p-2p?
    print(f"\nMethod 4 (completely standard formula):")
    E_2e_std = 0.0
    for name_i in ['1s', '2s', '2p']:
        n_i = n[name_i]
        for name_j in ['1s', '2s', '2p']:
            n_j = n[name_j]
            # Standard: 0.5 * n_i * n_j * (J - 0.5K)
            E_2e_std += 0.5 * n_i * n_j * (J_mat[(name_i, name_j)] - 0.5 * K_mat[(name_i, name_j)])
    
    print(f"  E_2e = {E_2e_std:.8f} Ha")
    print(f"  E_tot = {E_1e + E_2e_std:.8f} Ha")
    
    # Reference breakdown
    print("\n" + "=" * 70)
    print("REFERENCE VALUES")
    print("=" * 70)
    print(f"E_tot (ref)  = {C_data['E_tot']:.8f} Ha")
    print(f"ε_1s (ref)   = {C_data['energies']['1s']:.8f} Ha")
    print(f"ε_2s (ref)   = {C_data['energies']['2s']:.8f} Ha")
    print(f"ε_2p (ref)   = {C_data['energies']['2p']:.8f} Ha")


if __name__ == "__main__":
    main()
