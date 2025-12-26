# -*- coding: utf-8 -*-
"""
Completely diagnose Zn 3d orbital energy components.
Calculate T, V_nuc, J_total, and breakdown of K contributions.
Compare with reference energy to find the missing piece.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
import json, re
from pathlib import Path
from dataclasses import dataclass
from typing import List
from scipy.special import factorial

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

def sto_normalization(n, zeta):
    return (2 * zeta) ** n * np.sqrt(2 * zeta / float(factorial(2 * n)))

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

def create_log_grid(r_min=1e-7, r_max=100.0, N=5000):
    x = np.linspace(0, np.log(r_max / r_min), N)
    return r_min * np.exp(x)

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

def compute_Yk(P_i, P_j, r, k=0):
    rho = P_i * P_j
    N = len(r)
    safe_r = np.where(r > 1e-20, r, 1e-20)
    integrand_inner = rho * r**k
    I_inner = np.zeros(N)
    I_inner[1:] = cumulative_trapezoid(integrand_inner, r)
    integrand_outer = rho / safe_r**(k+1)
    total_outer = np.trapezoid(integrand_outer, r)
    I_outer = np.zeros(N)
    I_outer[0] = total_outer
    I_outer[1:] = total_outer - cumulative_trapezoid(integrand_outer, r)
    return I_inner / safe_r**(k+1) + r**k * I_outer

def compute_kinetic_energy(orbital, R, r):
    l = orbital.l
    R_prime = np.gradient(R, r)
    R_double_prime = np.gradient(R_prime, r)
    P = r * R
    P_prime = R + r * R_prime
    P_double_prime = 2 * R_prime + r * R_double_prime
    safe_r = np.where(r > 1e-15, r, 1e-15)
    integrand = P * (P_double_prime - l*(l+1)*P/safe_r**2)
    T = -0.5 * np.trapezoid(integrand, r)
    return T

# BARE coefficients (currently used)
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

# Load Data
data = parse_slater_basis(Path('../slater_basis.js'))
r = create_log_grid(r_min=1e-7, r_max=100.0, N=5000)

symbol = 'Zn'
atom_data = data[symbol]
config = {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10}
Z = atom_data['Z']

print("Calculating Zn 3d components breakdown:")
print("=" * 60)

orbitals = {}
for orb_name, occ in config.items():
    terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
             for t in atom_data['orbitals'][orb_name]]
    orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
    R = compute_orbital_radial(orb, r)
    P = r * R
    orbitals[orb_name] = (orb, R, P, occ)

target_name = '3d'
target, R_target, P_target, n_target = orbitals[target_name]
l_target = target.l

# 1. Kinetic Energy
T = compute_kinetic_energy(target, R_target, r)
print(f"T = {T:.6f} Ha")

# 2. Nuclear Potential
safe_r = np.where(r > 1e-15, r, 1e-15)
V_nuc = np.trapezoid(P_target**2 * (-Z/safe_r), r)
print(f"V_nuc = {V_nuc:.6f} Ha")

# 3. Hartree (J)
print("\nHartree (J) Breakdown:")
E_hartree = 0
for source_name, (source, R_j, P_j, n_j) in orbitals.items():
    Y0 = compute_Yk(P_j, P_j, r, k=0)
    F0 = np.trapezoid(P_target**2 * Y0, r)
    term = n_j * F0
    E_hartree += term
    print(f"  J({target_name}-{source_name}): N={n_j}, F0={F0:.6f} => {term:.6f}")
print(f"Total J = {E_hartree:.6f} Ha")

# 4. Exchange (K)
print("\nExchange (K) Breakdown (using BARE coeffs):")
E_exchange = 0
for source_name, (source, R_j, P_j, n_j) in orbitals.items():
    l_source = source.l
    coeffs = EXCHANGE_COEFFS.get((l_target, l_source), [(0, 1.0)])
    
    term_k = 0
    print(f"  K({target_name}-{source_name}):")
    for k, c_k in coeffs:
        Yk = compute_Yk(P_target, P_j, r, k=k)
        Gk = np.trapezoid(P_target * P_j * Yk, r)
        # Formula: (n_j/2) * c_k * Gk
        contrib = (n_j / 2) * c_k * Gk
        term_k += contrib
        print(f"    k={k}: c_k={c_k:.4f}, G^k={Gk:.6f} => {contrib:.6f}")
    
    E_exchange += term_k
    print(f"    Sum = {term_k:.6f}")

print(f"Total K = {E_exchange:.6f} Ha")

# Total Energy and Ref
eps_calc = T + V_nuc + E_hartree - E_exchange
eps_ref = atom_data['energies'][target_name]
print("-" * 60)
print(f"Calculated epsilon = {eps_calc:.6f} Ha")
print(f"Reference epsilon  = {eps_ref:.6f} Ha")
diff = eps_calc - eps_ref
print(f"Difference (Needs to be removed) = {diff:.6f} Ha")
print(f"Missing Exchange energy = {diff:.6f} Ha (K should be larger)")

# Ratio check
ratio = (E_exchange + diff) / E_exchange
print(f"\nK_required / K_current = {ratio:.4f}")
