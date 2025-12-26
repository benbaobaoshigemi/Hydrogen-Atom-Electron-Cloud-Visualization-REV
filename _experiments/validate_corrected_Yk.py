# -*- coding: utf-8 -*-
"""
CORRECTED Orbital Energy Validation
Bug fix: V = Y_k (NOT Y_k / r)
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
import json, re
from pathlib import Path
from dataclasses import dataclass
from typing import List

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

def create_log_grid(r_min=1e-7, r_max=100.0, N=5000):
    x = np.linspace(0, np.log(r_max / r_min), N)
    return r_min * np.exp(x)

def factorial(n):
    if n <= 1: return 1.0
    result = 1.0
    for i in range(2, n + 1): result *= i
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
    return 0

# BARE 3j-squared coefficients
EXCHANGE_COEFFS = {
    (0, 0): [(0, 1.0)],
    (0, 1): [(1, 1/3)],
    (1, 0): [(1, 1/3)],
    (1, 1): [(0, 1/3), (2, 2/15)],
    (0, 2): [(2, 1/5)],
    (2, 0): [(2, 1/5)],
    (1, 2): [(1, 2/15), (3, 3/35)],
    (2, 1): [(1, 2/15), (3, 3/35)],
    (2, 2): [(0, 1/5), (2, 2/35), (4, 2/35)],
}

def compute_Yk(P_i, P_j, r, k=0):
    """
    Compute Y_k screening function (for P = r*R functions)
    Y_k(r) = I_inner(r) / r^(k+1) + r^k * I_outer(r)
    where:
        I_inner(r) = ∫_0^r P_i*P_j * s^k ds
        I_outer(r) = ∫_r^∞ P_i*P_j / s^(k+1) ds
    """
    rho = P_i * P_j
    N = len(r)
    safe_r = np.where(r > 1e-20, r, 1e-20)
    
    # Inner integral: ∫_0^r ρ * s^k ds
    integrand_inner = rho * r**k
    I_inner = np.zeros(N)
    I_inner[1:] = cumulative_trapezoid(integrand_inner, r)
    
    # Outer integral: ∫_r^∞ ρ / s^(k+1) ds
    integrand_outer = rho / safe_r**(k+1)
    total_outer = np.trapezoid(integrand_outer, r)
    I_outer = np.zeros(N)
    I_outer[0] = total_outer
    I_outer[1:] = total_outer - cumulative_trapezoid(integrand_outer, r)
    
    # Y_k = I_inner / r^(k+1) + r^k * I_outer
    Y = I_inner / safe_r**(k+1) + r**k * I_outer
    
    return Y

def compute_kinetic_energy(orbital, R, r):
    """Compute kinetic energy integral: <nl|T|nl>"""
    l = orbital.l
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
        chi_double_prime = N * ((n - 1) * (n - 2) * r_nm3 - 2 * zeta * (n - 1) * r_nm2 + zeta**2 * r_nm1) * exp_term
        R_prime += c * chi_prime
        R_double_prime += c * chi_double_prime
    
    safe_r = np.where(r > 1e-15, r, 1e-15)
    
    # T = -1/2 * <R | ∇² | R> = -1/2 * ∫ R * [R'' + 2R'/r - l(l+1)R/r²] * r² dr
    P = r * R
    P_prime = R + r * R_prime
    P_double_prime = 2 * R_prime + r * R_double_prime
    
    # Using P: T = -1/2 * ∫ P * [P'' - l(l+1)P/r²] dr
    integrand = P * (P_double_prime - l*(l+1)*P/safe_r**2)
    T = -0.5 * np.trapezoid(integrand, r)
    
    return T

def parse_slater_basis(filepath):
    content = filepath.read_text(encoding='utf-8')
    start = content.find('globalScope.SlaterBasis = ')
    json_start = start + len('globalScope.SlaterBasis = ')
    depth, json_end, in_str, escape = 0, json_start, False, False
    for i, c in enumerate(content[json_start:]):
        if escape: escape = False; continue
        if c == '\\': escape = True; continue
        if c == '"': in_str = not in_str; continue
        if in_str: continue
        if c == '{': depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0: json_end = json_start + i + 1; break
    json_str = content[json_start:json_end]
    json_str = re.sub(r',\s*}', '}', json_str)
    json_str = re.sub(r',\s*]', ']', json_str)
    return json.loads(json_str)

# Test
data = parse_slater_basis(Path('../slater_basis.js'))
r = create_log_grid(r_min=1e-7, r_max=100.0, N=5000)

print("CORRECTED FORMULA: V = Y_k (not Y_k/r)")
print("=" * 70)

CONFIGS = {
    'He': {'1s': 2},
    'Be': {'1s': 2, '2s': 2},
    'Ne': {'1s': 2, '2s': 2, '2p': 6},
    'Ar': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6},
    'Zn': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10},
}

results = []
for symbol, config in CONFIGS.items():
    if symbol not in data:
        continue
    
    atom_data = data[symbol]
    Z = atom_data['Z']
    
    orbitals = {}
    for orb_name, occ in config.items():
        terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
                 for t in atom_data['orbitals'][orb_name]]
        orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
        R = compute_orbital_radial(orb, r)
        P = r * R  # P = r*R
        orbitals[orb_name] = (orb, R, P, occ)
    
    print(f"\n{symbol} (Z={Z}):")
    
    for target_name in config.keys():
        target, R_target, P_target, n_target = orbitals[target_name]
        l_target = target.l
        
        # T = kinetic energy
        T = compute_kinetic_energy(target, R_target, r)
        
        # V_nuc = ∫ P² * (-Z/r) dr = ∫ R² * r² * (-Z/r) dr
        safe_r = np.where(r > 1e-15, r, 1e-15)
        V_nuc = np.trapezoid(P_target**2 * (-Z/safe_r), r)
        
        # Hartree (Coulomb): J = n_j * ∫ P_i² * Y_0(P_j, P_j) dr
        E_hartree = 0
        for source_name, (source, R_j, P_j, n_j) in orbitals.items():
            Y0 = compute_Yk(P_j, P_j, r, k=0)
            J_ij = np.trapezoid(P_target**2 * Y0, r)
            E_hartree += n_j * J_ij
        
        # Exchange: K = (n_j/2) * Σ_k c_k * ∫ P_i * P_j * Y_k(P_i, P_j) dr
        E_exchange = 0
        for source_name, (source, R_j, P_j, n_j) in orbitals.items():
            l_source = source.l
            coeffs = EXCHANGE_COEFFS.get((l_target, l_source), [(0, 1.0)])
            for k, c_k in coeffs:
                Yk = compute_Yk(P_target, P_j, r, k=k)
                K_k = np.trapezoid(P_target * P_j * Yk, r)
                E_exchange += c_k * (n_j / 2) * K_k
        
        eps_calc = T + V_nuc + E_hartree - E_exchange
        eps_ref = atom_data['energies'].get(target_name, 0)
        error = 100 * (eps_calc - eps_ref) / abs(eps_ref) if eps_ref != 0 else 0
        
        status = "✓" if abs(error) < 1 else ("~" if abs(error) < 5 else "✗")
        print(f"  {target_name}: eps={eps_calc:.6f}, ref={eps_ref:.6f}, err={error:.2f}% {status}")
        results.append((symbol, target_name, error))

print("\n" + "=" * 70)
print("SUMMARY:")
for s in CONFIGS.keys():
    atom_errs = [abs(e) for (sym, _, e) in results if sym == s]
    if atom_errs:
        status = "PASS" if max(atom_errs) < 1 else ("MARGINAL" if max(atom_errs) < 5 else "FAIL")
        print(f"  {s}: max error = {max(atom_errs):.2f}% {status}")
