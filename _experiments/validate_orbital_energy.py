# -*- coding: utf-8 -*-
"""
ORBITAL ENERGY DENSITY VALIDATION
==================================

This script validates that the integral of orbital energy density E_i(r)
converges to the reference orbital eigenvalue ε_i from Koga basis.

Target: ∫ E_i(r) r² dr → ε_i

This is a simpler validation than total energy, focusing on individual orbitals.
"""

import numpy as np
import json
import re
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict

# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

@dataclass
class STOTerm:
    n: int
    zeta: float
    coeff: float


@dataclass
class Orbital:
    name: str
    l: int
    terms: List[STOTerm]


# ==============================================================================
# GRID AND STO FUNCTIONS
# ==============================================================================

def create_log_grid(r_min=1e-6, r_max=50.0, N=2000):
    x = np.linspace(0, np.log(r_max / r_min), N)
    return r_min * np.exp(x)


def factorial(n):
    if n <= 1:
        return 1.0
    result = 1.0
    for i in range(2, n + 1):
        result *= i
    return result


def sto_normalization(n, zeta):
    return (2 * zeta) ** n * np.sqrt(2 * zeta / factorial(2 * n))


def compute_orbital_radial(orbital, r):
    R = np.zeros_like(r)
    for term in orbital.terms:
        N = sto_normalization(term.n, term.zeta)
        R += term.coeff * N * np.power(r, term.n - 1) * np.exp(-term.zeta * r)
    return R


# ==============================================================================
# Y_k POTENTIAL (k=0, 1, 2)
# ==============================================================================

# Slater-Condon angular coefficients for exchange
EXCHANGE_COEFFS = {
    (0, 0): [(0, 1.0)],           # s-s
    (0, 1): [(1, 1/3)],           # s-p
    (1, 0): [(1, 1/3)],           # p-s
    (1, 1): [(0, 1.0), (2, 2/5)], # p-p: k=0 (coeff 1) + k=2 (coeff 2/5)
    (0, 2): [(2, 1/5)],           # s-d
    (2, 0): [(2, 1/5)],           # d-s
    (1, 2): [(1, 2/15), (3, 3/35)], # p-d
    (2, 1): [(1, 2/15), (3, 3/35)], # d-p
    (2, 2): [(0, 1.0), (2, 2/7), (4, 2/7)], # d-d
}


def compute_Yk(R_i, R_j, r, k=0):
    """
    Compute Y_k(r) for any k.
    
    Y_0(r) = ∫_0^r ρ r'^2 dr' + r ∫_r^∞ ρ r' dr'
    Y_k(r) = (1/r^k) ∫_0^r ρ r'^(k+2) dr' + r^(k+1) ∫_r^∞ ρ r'^(1-k) dr'
    """
    N = len(r)
    rho = R_i * R_j
    
    if k == 0:
        integrand_inner = rho * r ** 2
        cumsum_inner = np.zeros(N)
        for idx in range(1, N):
            cumsum_inner[idx] = np.trapz(integrand_inner[:idx+1], r[:idx+1])
        
        integrand_outer = rho * r
        total_outer = np.trapz(integrand_outer, r)
        cumsum_outer = np.zeros(N)
        for idx in range(N):
            if idx == 0:
                cumsum_outer[idx] = total_outer
            else:
                cumsum_outer[idx] = total_outer - np.trapz(integrand_outer[:idx+1], r[:idx+1])
        
        return cumsum_inner + r * cumsum_outer
    else:
        Y = np.zeros(N)
        for idx in range(N):
            r_val = r[idx]
            if r_val < 1e-15:
                continue
            
            inner_integrand = rho[:idx+1] * r[:idx+1] ** (k + 2)
            I_inner = np.trapz(inner_integrand, r[:idx+1]) if idx > 0 else 0.0
            
            if k <= 1:
                outer_integrand = rho[idx:] * r[idx:] ** (1 - k)
            else:
                outer_integrand = rho[idx:] / (r[idx:] ** (k - 1))
            
            I_outer = np.trapz(outer_integrand, r[idx:]) if len(outer_integrand) > 1 else 0.0
            
            Y[idx] = I_inner / (r_val ** k) + (r_val ** (k + 1)) * I_outer
        
        return Y


# ==============================================================================
# ENERGY DENSITY COMPONENTS
# ==============================================================================

def compute_kinetic_density(orbital, R, r):
    """T density = R * (-1/2) * ∇²R"""
    l = orbital.l
    
    R_val = np.zeros_like(r)
    R_prime = np.zeros_like(r)
    R_double_prime = np.zeros_like(r)
    
    for term in orbital.terms:
        n = term.n
        zeta = term.zeta
        c = term.coeff
        N = sto_normalization(n, zeta)
        
        safe_r = np.where(r > 1e-15, r, 1e-15)
        exp_term = np.exp(-zeta * r)
        
        r_nm1 = np.power(safe_r, n - 1)
        r_nm2 = np.power(safe_r, n - 2) if n >= 2 else np.zeros_like(r)
        r_nm3 = np.power(safe_r, n - 3) if n >= 3 else np.zeros_like(r)
        
        chi = N * r_nm1 * exp_term
        chi_prime = N * ((n - 1) * r_nm2 - zeta * r_nm1) * exp_term
        chi_double_prime = N * (
            (n - 1) * (n - 2) * r_nm3
            - 2 * zeta * (n - 1) * r_nm2
            + zeta ** 2 * r_nm1
        ) * exp_term
        
        R_val += c * chi
        R_prime += c * chi_prime
        R_double_prime += c * chi_double_prime
    
    safe_r = np.where(r > 1e-15, r, 1e-15)
    laplacian = R_double_prime + (2 / safe_r) * R_prime - l * (l + 1) / (safe_r ** 2) * R_val
    
    return R_val * (-0.5) * laplacian


def compute_nuclear_density(R, r, Z):
    safe_r = np.where(r > 1e-15, r, 1e-15)
    return R ** 2 * (-Z / safe_r)


def compute_hartree_density(R_i, R_j, r, n_j):
    """Coulomb potential density from orbital j acting on orbital i"""
    Y0 = compute_Yk(R_j, R_j, r, k=0)
    V = np.zeros_like(r)
    mask = r > 1e-10
    V[mask] = Y0[mask] / r[mask]
    return R_i ** 2 * V * n_j


def compute_exchange_density(R_i, R_j, l_i, l_j, r, n_j):
    """
    Exchange density with Slater-Condon angular coefficients.
    
    K = Σ_k c_k * (R_i R_j) * Y_k(i,j) / r
    """
    coeffs = EXCHANGE_COEFFS.get((l_i, l_j), [(0, 1.0)])
    
    E_exch = np.zeros_like(r)
    mask = r > 1e-10
    
    for k, c_k in coeffs:
        Yk = compute_Yk(R_i, R_j, r, k=k)
        V_k = np.zeros_like(r)
        V_k[mask] = Yk[mask] / r[mask]
        E_exch += c_k * R_i * R_j * V_k
    
    return E_exch * n_j * 0.5  # 0.5 for spin average


# ==============================================================================
# PARSING
# ==============================================================================

def parse_slater_basis(filepath):
    content = filepath.read_text(encoding='utf-8')
    start = content.find('globalScope.SlaterBasis = ')
    json_start = start + len('globalScope.SlaterBasis = ')
    
    depth, json_end, in_str, escape = 0, json_start, False, False
    for i, c in enumerate(content[json_start:]):
        if escape:
            escape = False
            continue
        if c == '\\':
            escape = True
            continue
        if c == '"':
            in_str = not in_str
            continue
        if in_str:
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                json_end = json_start + i + 1
                break
    
    json_str = content[json_start:json_end]
    json_str = re.sub(r',\s*}', '}', json_str)
    json_str = re.sub(r',\s*]', ']', json_str)
    
    return json.loads(json_str)


def get_l(name):
    if 's' in name: return 0
    if 'p' in name: return 1
    if 'd' in name: return 2
    return 0


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 70)
    print("ORBITAL ENERGY DENSITY VALIDATION")
    print("∫ E_i(r) r² dr → ε_i (reference orbital energy)")
    print("=" * 70)
    
    basis_file = Path(__file__).parent.parent / "slater_basis.js"
    data = parse_slater_basis(basis_file)
    r = create_log_grid(r_min=1e-6, r_max=50.0, N=2000)
    
    # Test closed-shell atoms only
    test_cases = [
        # He (Z=2): 1s²
        ('He', '1s', {'1s': 2}),
        
        # Be (Z=4): 1s² 2s²
        ('Be', '1s', {'1s': 2, '2s': 2}),
        ('Be', '2s', {'1s': 2, '2s': 2}),
        
        # Ne (Z=10): 1s² 2s² 2p⁶
        ('Ne', '1s', {'1s': 2, '2s': 2, '2p': 6}),
        ('Ne', '2s', {'1s': 2, '2s': 2, '2p': 6}),
        ('Ne', '2p', {'1s': 2, '2s': 2, '2p': 6}),
        
        # Mg (Z=12): 1s² 2s² 2p⁶ 3s²
        ('Mg', '1s', {'1s': 2, '2s': 2, '2p': 6, '3s': 2}),
        ('Mg', '2s', {'1s': 2, '2s': 2, '2p': 6, '3s': 2}),
        ('Mg', '2p', {'1s': 2, '2s': 2, '2p': 6, '3s': 2}),
        ('Mg', '3s', {'1s': 2, '2s': 2, '2p': 6, '3s': 2}),
        
        # Ar (Z=18): 1s² 2s² 2p⁶ 3s² 3p⁶
        ('Ar', '1s', {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6}),
        ('Ar', '2s', {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6}),
        ('Ar', '2p', {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6}),
        ('Ar', '3s', {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6}),
        ('Ar', '3p', {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6}),
    ]
    
    print(f"\n{'Atom':>4} {'Orb':>4} {'ε_calc':>12} {'ε_ref':>12} {'Error%':>10}")
    print("-" * 50)
    
    for symbol, target_orb, config in test_cases:
        atom_data = data[symbol]
        Z = atom_data['Z']
        
        # Build all orbitals
        orbitals = {}
        for orb_name, occ in config.items():
            terms = [
                STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff'])
                for t in atom_data['orbitals'][orb_name]
            ]
            orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
            R = compute_orbital_radial(orb, r)
            orbitals[orb_name] = (orb, R, occ)
        
        # Compute energy density for target orbital
        target, R_target, n_target = orbitals[target_orb]
        
        E_kin = compute_kinetic_density(target, R_target, r)
        E_nuc = compute_nuclear_density(R_target, r, Z)
        
        E_hartree = np.zeros_like(r)
        E_exchange = np.zeros_like(r)
        
        l_target = target.l
        for orb_name, (orb, R_j, n_j) in orbitals.items():
            E_hartree += compute_hartree_density(R_target, R_j, r, n_j)
            E_exchange += compute_exchange_density(R_target, R_j, l_target, orb.l, r, n_j)
        
        E_total = E_kin + E_nuc + E_hartree - E_exchange
        
        # Integrate
        eps_calc = np.trapz(E_total * r ** 2, r)
        eps_ref = atom_data['energies'][target_orb]
        error = 100 * (eps_calc - eps_ref) / abs(eps_ref)
        
        print(f"{symbol:>4} {target_orb:>4} {eps_calc:>12.6f} {eps_ref:>12.6f} {error:>10.4f}")
    
    print("-" * 50)
    print("\nNote: Errors should be < 1% for validated implementation.")


if __name__ == "__main__":
    main()
