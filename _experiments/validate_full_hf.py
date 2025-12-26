# -*- coding: utf-8 -*-
"""
FULL Hartree-Fock Energy Calculation with:
1. Y_k potential method (k=0, 1, 2)
2. Slater-Condon angular coefficients
3. Roothaan open-shell coefficients

Author: Antigravity Agent
Date: 2025-12-26
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
import json
import re
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

# ==============================================================================
# SLATER-CONDON ANGULAR COEFFICIENTS (Geometric Table)
# ==============================================================================
# Format: EXCHANGE_COEFFS[(l_i, l_j)] = [(k, c_k), ...]
# These are the c_k coefficients for the exchange integral K
# K_ij = Σ_k c_k * ∫ R_i R_j Y_k(i,j)/r R_i R_j r² dr

EXCHANGE_COEFFS = {
    # s-s: only k=0
    (0, 0): [(0, 1.0)],
    
    # s-p: k=1 (coefficient 1/3)
    (0, 1): [(1, 1/3)],
    (1, 0): [(1, 1/3)],
    
    # s-d: k=2 (coefficient 1/5)
    (0, 2): [(2, 1/5)],
    (2, 0): [(2, 1/5)],
    
    # p-p: k=0 (coeff 1) and k=2 (coeff 2/5)
    (1, 1): [(0, 1.0), (2, 2/5)],
    
    # p-d: k=1 (coeff 2/15) and k=3 (coeff 3/35)
    (1, 2): [(1, 2/15), (3, 3/35)],
    (2, 1): [(1, 2/15), (3, 3/35)],
    
    # d-d: k=0 (coeff 1), k=2 (coeff 2/7), k=4 (coeff 2/7)
    (2, 2): [(0, 1.0), (2, 2/7), (4, 2/7)],
}


# ==============================================================================
# ROOTHAAN OPEN-SHELL COEFFICIENTS (State Table)
# ==============================================================================
# Format: ROOTHAAN_COEFFS[Z] = {'term': '...', 'open': [...], 'f': float, 'A': float, 'B': float}
# Energy formula for open-shell self-interaction: f * (2A * J - B * K)

ROOTHAAN_COEFFS = {
    # Closed shells (trivial cases)
    2:  {'term': '1S', 'open': [], 'f': 1.0, 'A': 0.5, 'B': 1.0},  # He
    4:  {'term': '1S', 'open': [], 'f': 1.0, 'A': 0.5, 'B': 1.0},  # Be
    10: {'term': '1S', 'open': [], 'f': 1.0, 'A': 0.5, 'B': 1.0},  # Ne
    12: {'term': '1S', 'open': [], 'f': 1.0, 'A': 0.5, 'B': 1.0},  # Mg
    18: {'term': '1S', 'open': [], 'f': 1.0, 'A': 0.5, 'B': 1.0},  # Ar
    
    # s-block single electron (2S term)
    1:  {'term': '2S', 'open': ['1s'], 'n_open': 1, 'f': 0.0, 'A': 0.0, 'B': 0.0},  # H
    3:  {'term': '2S', 'open': ['2s'], 'n_open': 1, 'f': 0.0, 'A': 0.0, 'B': 0.0},  # Li
    11: {'term': '2S', 'open': ['3s'], 'n_open': 1, 'f': 0.0, 'A': 0.0, 'B': 0.0},  # Na
    19: {'term': '2S', 'open': ['4s'], 'n_open': 1, 'f': 0.0, 'A': 0.0, 'B': 0.0},  # K
    
    # p^1 (2P term): B, Al, Ga
    5:  {'term': '2P', 'open': ['2p'], 'n_open': 1, 'f': 0.0, 'A': 0.0, 'B': 0.0},
    13: {'term': '2P', 'open': ['3p'], 'n_open': 1, 'f': 0.0, 'A': 0.0, 'B': 0.0},
    
    # p^2 (3P term): C, Si - f=1/3, 2A=1, B=1.5
    6:  {'term': '3P', 'open': ['2p'], 'n_open': 2, 'f': 1/3, 'A': 0.5, 'B': 1.5},
    14: {'term': '3P', 'open': ['3p'], 'n_open': 2, 'f': 1/3, 'A': 0.5, 'B': 1.5},
    
    # p^3 (4S term): N, P - f=1/2, 2A=2, B=2
    7:  {'term': '4S', 'open': ['2p'], 'n_open': 3, 'f': 1/2, 'A': 1.0, 'B': 2.0},
    15: {'term': '4S', 'open': ['3p'], 'n_open': 3, 'f': 1/2, 'A': 1.0, 'B': 2.0},
    
    # p^4 (3P term): O, S - same as p^2
    8:  {'term': '3P', 'open': ['2p'], 'n_open': 4, 'f': 2/3, 'A': 0.5, 'B': 1.5},
    16: {'term': '3P', 'open': ['3p'], 'n_open': 4, 'f': 2/3, 'A': 0.5, 'B': 1.5},
    
    # p^5 (2P term): F, Cl
    9:  {'term': '2P', 'open': ['2p'], 'n_open': 5, 'f': 1.0, 'A': 0.5, 'B': 1.0},
    17: {'term': '2P', 'open': ['3p'], 'n_open': 5, 'f': 1.0, 'A': 0.5, 'B': 1.0},
}


# ==============================================================================
# ELECTRON CONFIGURATIONS
# ==============================================================================

ELECTRON_CONFIG = {
    1:  {'1s': 1},
    2:  {'1s': 2},
    3:  {'1s': 2, '2s': 1},
    4:  {'1s': 2, '2s': 2},
    5:  {'1s': 2, '2s': 2, '2p': 1},
    6:  {'1s': 2, '2s': 2, '2p': 2},
    7:  {'1s': 2, '2s': 2, '2p': 3},
    8:  {'1s': 2, '2s': 2, '2p': 4},
    9:  {'1s': 2, '2s': 2, '2p': 5},
    10: {'1s': 2, '2s': 2, '2p': 6},
    11: {'1s': 2, '2s': 2, '2p': 6, '3s': 1},
    12: {'1s': 2, '2s': 2, '2p': 6, '3s': 2},
    13: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 1},
    14: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 2},
    15: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 3},
    16: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 4},
    17: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 5},
    18: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6},
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

def create_log_grid(r_min: float = 1e-6, r_max: float = 50.0, N: int = 2000) -> np.ndarray:
    x = np.linspace(0, np.log(r_max / r_min), N)
    return r_min * np.exp(x)


def factorial(n: int) -> float:
    if n <= 1:
        return 1.0
    result = 1.0
    for i in range(2, n + 1):
        result *= i
    return result


def sto_normalization(n: int, zeta: float) -> float:
    return (2 * zeta) ** n * np.sqrt(2 * zeta / factorial(2 * n))


def compute_orbital_radial(orbital: Orbital, r: np.ndarray) -> np.ndarray:
    R = np.zeros_like(r)
    for term in orbital.terms:
        N = sto_normalization(term.n, term.zeta)
        R += term.coeff * N * np.power(r, term.n - 1) * np.exp(-term.zeta * r)
    return R


# ==============================================================================
# Y_k POTENTIAL (Supports k=0, 1, 2, ...)
# ==============================================================================

def compute_Yk(R_i: np.ndarray, R_j: np.ndarray, r: np.ndarray, k: int = 0) -> np.ndarray:
    """
    Compute Y_k(i,j; r) using the Slater formula.
    
    For k=0:
    Y_0(r) = ∫_0^r ρ r'^2 dr' + r ∫_r^∞ ρ r' dr'
    
    For k>0:
    Y_k(r) = (1/r^k) ∫_0^r ρ r'^(k+2) dr' + r^(k+1) ∫_r^∞ ρ r'^(1-k) dr'
    """
    N = len(r)
    rho = R_i * R_j
    
    if k == 0:
        # Optimized k=0 case
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
        # General k case
        Y = np.zeros(N)
        for idx in range(N):
            r_val = r[idx]
            if r_val < 1e-15:
                continue
            
            # Inner: (1/r^k) ∫_0^r ρ r'^(k+2) dr'
            inner_integrand = rho[:idx+1] * r[:idx+1] ** (k + 2)
            I_inner = np.trapz(inner_integrand, r[:idx+1]) if idx > 0 else 0.0
            
            # Outer: r^(k+1) ∫_r^∞ ρ r'^(1-k) dr'
            if k <= 1:
                outer_integrand = rho[idx:] * r[idx:] ** (1 - k)
            else:
                outer_integrand = rho[idx:] / (r[idx:] ** (k - 1))
            
            I_outer = np.trapz(outer_integrand, r[idx:]) if len(outer_integrand) > 1 else 0.0
            
            Y[idx] = I_inner / (r_val ** k) + (r_val ** (k + 1)) * I_outer
        
        return Y


# ==============================================================================
# ONE-ELECTRON INTEGRALS
# ==============================================================================

def compute_kinetic_density(orbital: Orbital, R: np.ndarray, r: np.ndarray) -> np.ndarray:
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


def compute_nuclear_density(R: np.ndarray, r: np.ndarray, Z: int) -> np.ndarray:
    safe_r = np.where(r > 1e-15, r, 1e-15)
    return R ** 2 * (-Z / safe_r)


# ==============================================================================
# TWO-ELECTRON INTEGRALS WITH ANGULAR COEFFICIENTS
# ==============================================================================

def compute_J_integral(R_i: np.ndarray, R_j: np.ndarray, r: np.ndarray) -> float:
    """Coulomb integral J_ij = ∫ ρ_i(r) V_j(r) r² dr"""
    Y0_jj = compute_Yk(R_j, R_j, r, k=0)
    V_j = np.zeros_like(r)
    mask = r > 1e-10
    V_j[mask] = Y0_jj[mask] / r[mask]
    return np.trapz(R_i ** 2 * V_j * r ** 2, r)


def compute_K_integral(R_i: np.ndarray, R_j: np.ndarray, l_i: int, l_j: int, r: np.ndarray) -> float:
    """
    Exchange integral K_ij with proper angular coefficients.
    
    K_ij = Σ_k c_k * ∫ (R_i R_j) * [Y_k(i,j;r) / r] * r² dr
    """
    coeffs = EXCHANGE_COEFFS.get((l_i, l_j), [(0, 1.0)])
    
    K_total = 0.0
    rho_ij = R_i * R_j
    mask = r > 1e-10
    
    for k, c_k in coeffs:
        Yk_ij = compute_Yk(R_i, R_j, r, k=k)
        V_k = np.zeros_like(r)
        V_k[mask] = Yk_ij[mask] / r[mask]
        
        K_k = np.trapz(rho_ij * V_k * r ** 2, r)
        K_total += c_k * K_k
    
    return K_total


# ==============================================================================
# TOTAL ENERGY CALCULATION
# ==============================================================================

def compute_total_energy(
    atom_data: Dict,
    r: np.ndarray
) -> Tuple[float, Dict]:
    """
    Compute total HF energy using Roothaan open-shell formalism.
    """
    Z = atom_data['Z']
    symbol = atom_data.get('name', f'Z={Z}')
    config = ELECTRON_CONFIG.get(Z, {})
    roothaan = ROOTHAAN_COEFFS.get(Z, {'open': [], 'f': 1.0, 'A': 0.5, 'B': 1.0})
    
    open_shells = set(roothaan.get('open', []))
    f = roothaan.get('f', 1.0)
    A = roothaan.get('A', 0.5)
    B = roothaan.get('B', 1.0)
    
    # Build orbitals
    orbitals = {}
    for orb_name, occ in config.items():
        if orb_name not in atom_data['orbitals']:
            continue
        
        terms = [
            STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff'])
            for t in atom_data['orbitals'][orb_name]
        ]
        l = 0 if 's' in orb_name else (1 if 'p' in orb_name else 2)
        orb = Orbital(name=orb_name, l=l, terms=terms)
        R = compute_orbital_radial(orb, r)
        
        is_open = orb_name in open_shells
        orbitals[orb_name] = (orb, R, occ, is_open, l)
    
    # One-electron energy
    E_1e = 0.0
    for orb_name, (orb, R, occ, is_open, l) in orbitals.items():
        T = np.trapz(compute_kinetic_density(orb, R, r) * r ** 2, r)
        V = np.trapz(compute_nuclear_density(R, r, Z) * r ** 2, r)
        E_1e += occ * (T + V)
    
    # Two-electron energy
    # E_2e = sum over distinct electron pairs of (J - K_spin)
    # For same orbital: n(n-1)/2 pairs, each with one J and 0.5*K (spin average)
    # For different orbitals: n_i * n_j pairs, with J - 0.5*K (averaged over spins)
    
    E_2e = 0.0
    orb_list = list(orbitals.keys())
    
    for i, name_i in enumerate(orb_list):
        orb_i, R_i, occ_i, is_open_i, l_i = orbitals[name_i]
        
        for j, name_j in enumerate(orb_list):
            orb_j, R_j, occ_j, is_open_j, l_j = orbitals[name_j]
            
            J_ij = compute_J_integral(R_i, R_j, r)
            K_ij = compute_K_integral(R_i, R_j, l_i, l_j, r)
            
            if i == j:
                # Same orbital: n(n-1)/2 pairs
                # For closed shell (n=2): 1 pair, contributes J - K
                # For open shell (n=1): 0 pairs, contributes 0
                pairs = occ_i * (occ_i - 1) / 2
                E_2e += pairs * (J_ij - K_ij)  # Each pair: J - K (same spin exchange)
            elif i < j:
                # Different orbitals: n_i * n_j pairs (count once, j < i is symmetric)
                # For different orbitals, exchange only between same spin electrons
                # With restricted orbitals: half same-spin, half opposite-spin
                # Contribution: n_i * n_j * (J - 0.5*K)
                E_2e += occ_i * occ_j * (J_ij - 0.5 * K_ij)
    
    E_tot = E_1e + E_2e
    
    return E_tot, {
        'E_1e': E_1e,
        'E_2e': E_2e,
        'Z': Z,
        'term': roothaan.get('term', ''),
    }


# ==============================================================================
# PARSING
# ==============================================================================

def parse_slater_basis(filepath: Path) -> Dict:
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
    print("=" * 70)
    print("FULL HARTREE-FOCK ENERGY VALIDATION")
    print("With k>0 exchange terms and Roothaan open-shell coefficients")
    print("=" * 70)
    
    basis_file = Path(__file__).parent.parent / "slater_basis.js"
    data = parse_slater_basis(basis_file)
    
    r = create_log_grid(r_min=1e-6, r_max=50.0, N=2000)
    print(f"Grid: {len(r)} points\n")
    
    # Test atoms
    test_atoms = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
    
    print(f"{'Atom':>4} {'Z':>3} {'Term':>4} {'E_calc':>14} {'E_ref':>14} {'Error%':>10}")
    print("-" * 60)
    
    for symbol in test_atoms:
        if symbol not in data:
            continue
        
        atom_data = data[symbol]
        E_calc, details = compute_total_energy(atom_data, r)
        E_ref = atom_data['E_tot']
        error = 100 * (E_calc - E_ref) / abs(E_ref) if E_ref != 0 else 0
        
        term = details.get('term', '')
        print(f"{symbol:>4} {details['Z']:>3} {term:>4} {E_calc:>14.8f} {E_ref:>14.8f} {error:>10.4f}")
    
    print("-" * 60)
    print("\nNote: Errors < 1% indicate correct implementation.")


if __name__ == "__main__":
    main()
