# -*- coding: utf-8 -*-
"""
Calculate accurate F^k integrals for Zn 3d using Koga basis.
Used to reverse engineer the correct exchange coefficients.
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
    terms: List[STOTerm]

def sto_normalization(n, zeta):
    return (2 * zeta) ** n * np.sqrt(2 * zeta / float(factorial(2 * n)))

def compute_orbital_radial(orbital, r):
    R = np.zeros_like(r)
    for term in orbital.terms:
        N = sto_normalization(term.n, term.zeta)
        R += term.coeff * N * np.power(r, term.n - 1) * np.exp(-term.zeta * r)
    return R

def create_log_grid(r_min=1e-7, r_max=100.0, N=10000):
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

# Main calculation
try:
    data = parse_slater_basis(Path('../slater_basis.js'))
    print("Loaded basis data successfully.")
except Exception as e:
    print(f"Error loading basis data: {e}")
    exit(1)

r = create_log_grid(r_min=1e-7, r_max=100.0, N=10000)

print("Calculating Zn 3d integrals...")
terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
         for t in data['Zn']['orbitals']['3d']]
orb_3d = Orbital(name='3d', terms=terms)
R_3d = compute_orbital_radial(orb_3d, r)
P_3d = r * R_3d
norm = np.trapezoid(P_3d**2, r)
print(f'3d Norm: {norm:.6f}')

# Compute F^2 and F^4
# F^k = ∫ P^2 * Y_k dr

Y2 = compute_Yk(P_3d, P_3d, r, k=2)
F2 = np.trapezoid(P_3d**2 * Y2, r)

Y4 = compute_Yk(P_3d, P_3d, r, k=4)
F4 = np.trapezoid(P_3d**2 * Y4, r)

print(f'F^2(3d,3d) = {F2:.6f}')
print(f'F^4(3d,3d) = {F4:.6f}')

# Current BARE coefficient contribution
# BARE for d-d: k=2 -> 2/35, k=4 -> 2/35
# Formula: K = 5 * (BARE_2 * F2 + BARE_4 * F4)

K_calc_BARE_2 = 5 * (2/35) * F2
K_calc_BARE_4 = 5 * (2/35) * F4
K_total_calc = K_calc_BARE_2 + K_calc_BARE_4

print(f'K_calc (BARE, k=2) = {K_calc_BARE_2:.6f}')
print(f'K_calc (BARE, k=4) = {K_calc_BARE_4:.6f}')
print(f'K_total_calc (self K>0) = {K_total_calc:.6f}')
