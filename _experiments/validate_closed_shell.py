# -*- coding: utf-8 -*-
"""
FULL CLOSED-SHELL VALIDATION (H to Kr)
======================================

Verify orbital energy density integration for ALL closed-shell atoms 
from Z=1 to Z=36.

Closed-shell atoms:
- He (Z=2): 1s²
- Be (Z=4): 1s² 2s²
- Ne (Z=10): 1s² 2s² 2p⁶
- Mg (Z=12): 1s² 2s² 2p⁶ 3s²
- Ar (Z=18): 1s² 2s² 2p⁶ 3s² 3p⁶
- Ca (Z=20): [Ar] 4s²
- Zn (Z=30): [Ar] 3d¹⁰ 4s²
- Kr (Z=36): [Ar] 3d¹⁰ 4s² 4p⁶
"""

import numpy as np
import json
import re
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict

# ==============================================================================
# CLOSED-SHELL ELECTRON CONFIGURATIONS
# ==============================================================================

CLOSED_SHELL_CONFIGS = {
    'He': {'1s': 2},
    'Be': {'1s': 2, '2s': 2},
    'Ne': {'1s': 2, '2s': 2, '2p': 6},
    'Mg': {'1s': 2, '2s': 2, '2p': 6, '3s': 2},
    'Ar': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6},
    'Ca': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2},
    'Zn': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10},
    'Kr': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10, '4p': 6},
}

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


def get_l(name):
    if 's' in name: return 0
    if 'p' in name: return 1
    if 'd' in name: return 2
    if 'f' in name: return 3
    return 0


# BARE 3j-squared coefficients - these are the correct values when using 0.5*n_j
# in the exchange formula (which implicitly provides the (2l'+1) factor)
EXCHANGE_COEFFS = {
    # s-s: bare = full (since 2l+1=1)
    (0, 0): [(0, 1.0)],
    
    # s-p and p-s: bare values
    (0, 1): [(1, 1/3)],
    (1, 0): [(1, 1/3)],
    
    # p-p: BARE = FULL / 3 (since 2l+1=3 for p)
    # Full: k=0 coeff=1, k=2 coeff=2/5
    # Bare: k=0 coeff=1/3, k=2 coeff=2/15
    (1, 1): [(0, 1/3), (2, 2/15)],
    
    # s-d and d-s: bare values  
    (0, 2): [(2, 1/5)],
    (2, 0): [(2, 1/5)],
    
    # p-d and d-p: bare values
    (1, 2): [(1, 2/15), (3, 3/35)],
    (2, 1): [(1, 2/15), (3, 3/35)],
    
    # d-d: BARE = FULL / 5 (since 2l+1=5 for d)
    # Full: k=0 coeff=1, k=2 coeff=2/7, k=4 coeff=2/7
    # Bare: k=0 coeff=1/5, k=2 coeff=2/35, k=4 coeff=2/35
    (2, 2): [(0, 1/5), (2, 2/35), (4, 2/35)],
}


def compute_Yk(R_i, R_j, r, k=0):
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
    Y0 = compute_Yk(R_j, R_j, r, k=0)
    V = np.zeros_like(r)
    mask = r > 1e-10
    V[mask] = Y0[mask] / r[mask]
    return R_i ** 2 * V * n_j


def compute_exchange_density(R_i, R_j, l_i, l_j, r, n_j):
    """
    Exchange density with Slater-Condon coefficients.
    Using bare 3j-squared coefficients without additional degeneracy factors.
    """
    coeffs = EXCHANGE_COEFFS.get((l_i, l_j), [(0, 1.0)])
    
    E_exch = np.zeros_like(r)
    mask = r > 1e-10
    
    for k, c_k in coeffs:
        Yk = compute_Yk(R_i, R_j, r, k=k)
        V_k = np.zeros_like(r)
        V_k[mask] = Yk[mask] / r[mask]
        E_exch += c_k * R_i * R_j * V_k
    
    return E_exch * n_j * 0.5


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


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 80)
    print("CLOSED-SHELL ORBITAL ENERGY VALIDATION (H to Kr)")
    print("=" * 80)
    
    basis_file = Path(__file__).parent.parent / "slater_basis.js"
    data = parse_slater_basis(basis_file)
    r = create_log_grid(r_min=1e-6, r_max=80.0, N=3000)  # Extended for heavier atoms
    
    all_results = []
    
    for symbol, config in CLOSED_SHELL_CONFIGS.items():
        if symbol not in data:
            print(f"\n{symbol}: NOT IN BASIS FILE")
            continue
        
        atom_data = data[symbol]
        Z = atom_data['Z']
        
        print(f"\n{'='*80}")
        print(f"{symbol} (Z={Z})")
        print(f"{'='*80}")
        
        # Check if all required orbitals exist
        missing = [orb for orb in config if orb not in atom_data['orbitals']]
        if missing:
            print(f"  Missing orbitals: {missing}")
            continue
        
        # Build all orbitals
        orbitals = {}
        for orb_name, occ in config.items():
            terms = [
                STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff'])
                for t in atom_data['orbitals'][orb_name]
            ]
            orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
            R = compute_orbital_radial(orb, r)
            
            # Check normalization
            norm = np.trapz(R ** 2 * r ** 2, r)
            orbitals[orb_name] = (orb, R, occ, norm)
        
        print(f"\n{'Orbital':>6} {'Occ':>4} {'Norm':>10} {'ε_calc':>12} {'ε_ref':>12} {'Error%':>10}")
        print("-" * 70)
        
        for orb_name in config.keys():
            target, R_target, n_target, norm = orbitals[orb_name]
            
            E_kin = compute_kinetic_density(target, R_target, r)
            E_nuc = compute_nuclear_density(R_target, r, Z)
            
            E_hartree = np.zeros_like(r)
            E_exchange = np.zeros_like(r)
            
            l_target = target.l
            for other_name, (orb, R_j, n_j, _) in orbitals.items():
                E_hartree += compute_hartree_density(R_target, R_j, r, n_j)
                E_exchange += compute_exchange_density(R_target, R_j, l_target, orb.l, r, n_j)
            
            E_total = E_kin + E_nuc + E_hartree - E_exchange
            
            eps_calc = np.trapz(E_total * r ** 2, r)
            eps_ref = atom_data['energies'].get(orb_name, 0)
            
            if eps_ref != 0:
                error = 100 * (eps_calc - eps_ref) / abs(eps_ref)
            else:
                error = float('nan')
            
            status = "✓" if abs(error) < 1.0 else ("~" if abs(error) < 5.0 else "✗")
            
            print(f"{orb_name:>6} {n_target:>4} {norm:>10.6f} {eps_calc:>12.6f} {eps_ref:>12.6f} {error:>10.4f} {status}")
            
            all_results.append({
                'atom': symbol,
                'orbital': orb_name,
                'error': error
            })
        
        # Also compute total energy
        E_tot_ref = atom_data['E_tot']
        print(f"\n  E_tot (ref) = {E_tot_ref:.8f} Ha")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    for symbol in CLOSED_SHELL_CONFIGS.keys():
        atom_results = [r for r in all_results if r['atom'] == symbol]
        if not atom_results:
            continue
        
        max_error = max(abs(r['error']) for r in atom_results if not np.isnan(r['error']))
        status = "✓ PASS" if max_error < 1.0 else ("~ MARGINAL" if max_error < 5.0 else "✗ FAIL")
        print(f"  {symbol:>3}: max |error| = {max_error:>8.4f}% {status}")


if __name__ == "__main__":
    main()
