# -*- coding: utf-8 -*-
"""
Test grid convergence for Kinetic Energy of Zn
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

def compute_kinetic_energy(orbital, R, r):
    l = orbital.l
    R_prime = np.gradient(R, r)
    R_double_prime = np.gradient(R_prime, r)
    P = r * R
    safe_r = np.where(r > 1e-15, r, 1e-15)
    # T = -1/2 <P | P'' - l(l+1)/r^2 P>
    P_prime = R + r * R_prime
    P_double_prime = 2 * R_prime + r * R_double_prime
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

data = parse_slater_basis(Path('../slater_basis.js'))

print("Testing Grid Convergence for Zn Total Kinetic Energy")
print("Reference -E_tot = 1777.8 ha (Target T)")
print("=" * 60)

configs = [
    (1e-5, 50.0, 2000),
    (1e-6, 100.0, 5000),
    (1e-7, 100.0, 10000),
    (1e-8, 100.0, 20000),
    (1e-9, 200.0, 50000),
]

def get_l(name):
    if 's' in name: return 0
    if 'p' in name: return 1
    if 'd' in name: return 2
    return 0

config_electrons = {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10}
atom_data = data['Zn']

for r_min, r_max, N in configs:
    x = np.linspace(0, np.log(r_max / r_min), N)
    r = r_min * np.exp(x)
    
    T_total = 0
    for orb_name, n_elec in config_electrons.items():
        terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
                 for t in atom_data['orbitals'][orb_name]]
        orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
        R = compute_orbital_radial(orb, r)
        T = compute_kinetic_energy(orb, R, r)
        T_total += n_elec * T
    
    print(f"Grid: r_min={r_min:.1e}, r_max={r_max}, N={N}")
    print(f"  Total Kinetic Energy: {T_total:.4f} Ha")
