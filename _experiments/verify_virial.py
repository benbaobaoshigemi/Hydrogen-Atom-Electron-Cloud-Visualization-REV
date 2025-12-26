# -*- coding: utf-8 -*-
"""
Verify Zn total energy with virial theorem
This checks if the Koga basis set gives correct total energy
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

def create_log_grid(r_min=1e-7, r_max=100.0, N=5000):
    x = np.linspace(0, np.log(r_max / r_min), N)
    return r_min * np.exp(x)

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
    safe_r = np.where(r > 1e-15, r, 1e-15)
    P_double_prime = 2 * R_prime + r * R_double_prime
    integrand = P * (P_double_prime - l*(l+1)*P/safe_r**2)
    return -0.5 * np.trapezoid(integrand, r)

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

data = parse_slater_basis(Path('../slater_basis.js'))
r = create_log_grid(r_min=1e-7, r_max=100.0, N=5000)

print("Verifying Zn total energy with Koga basis")
print("=" * 60)

symbol = 'Zn'
config = {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10}
atom_data = data[symbol]
Z = atom_data['Z']

orbitals = {}
for orb_name, occ in config.items():
    terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
             for t in atom_data['orbitals'][orb_name]]
    orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
    R = compute_orbital_radial(orb, r)
    P = r * R
    norm = np.trapezoid(P**2, r)
    orbitals[orb_name] = (orb, R, P, occ, norm)
    print(f'{orb_name}: norm = {norm:.6f}')

print()

# Calculate total energy components
T_total = 0
V_nuc_total = 0
safe_r = np.where(r > 1e-15, r, 1e-15)

for orb_name, (orb, R, P, n, _) in orbitals.items():
    T = n * compute_kinetic_energy(orb, R, r)
    V_nuc = n * np.trapezoid(P**2 * (-Z/safe_r), r)
    T_total += T
    V_nuc_total += V_nuc

# Coulomb energy (careful not to double count)
V_ee_coulomb = 0
for i, (name_i, (orb_i, R_i, P_i, n_i, _)) in enumerate(orbitals.items()):
    for j, (name_j, (orb_j, R_j, P_j, n_j, _)) in enumerate(orbitals.items()):
        Y0 = compute_Yk(P_j, P_j, r, k=0)
        J = np.trapezoid(P_i**2 * Y0, r)
        if i == j:
            V_ee_coulomb += 0.5 * n_i * (n_i - 1) * J / 2  # Self-interaction correction
        else:
            V_ee_coulomb += 0.5 * n_i * n_j * J  # Different orbitals

# Exchange energy (skip for simplicity)
print(f'T_total = {T_total:.4f} Ha')
print(f'V_nuc_total = {V_nuc_total:.4f} Ha')
print(f'V_ee_coulomb = {V_ee_coulomb:.4f} Ha')
print()
print(f'E_tot (ref) = {atom_data["E_tot"]:.4f} Ha')
print()
print('Virial theorem: E = -T for a bound system')
print(f'-T = {-T_total:.4f} Ha')
