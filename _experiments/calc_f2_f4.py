# -*- coding: utf-8 -*-
"""
Calculate the energy difference between Average d9 and Ground State 2D term
using calculated Slater Integrals F2, F4 for Zn 3d.
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
symbol = 'Zn'
atom_data = data[symbol]

print(f"Calculating F2, F4 for Zn 3d...")
terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
             for t in atom_data['orbitals']['3d']]
orb = Orbital(name='3d', l=2, terms=terms)
R = compute_orbital_radial(orb, r)
P = r * R

# Calculate Fk = Integral P^2 Yk
Y2 = compute_Yk(P, P, r, k=2)
F2 = np.trapezoid(P**2 * Y2, r)

Y4 = compute_Yk(P, P, r, k=4)
F4 = np.trapezoid(P**2 * Y4, r)

print(f"F2 = {F2:.6f} Ha")
print(f"F4 = {F4:.6f} Ha")

# Calculated Gap from BARE Orbital E (-0.499) to Reference (-0.783)
# Gap = 0.284 Ha.

print("-" * 50)
print("Theoretical Energy Separation d9 (Average) - d9 (2D term)")
print("Standard Multiplet Formula (Condon-Shortley):")
# E(2D) is the ONLY term for d9 (or d1).
# E(d10 -> d9) Ionization Potential = difference in total energy.
# 
# Wait. If d9 has ONLY ONE term (2D), then Average Energy IS 2D Energy?
# NO.
# The "Average of Configuration" includes "phantom" states or is defined over the full set of determinants including those that don't satisfy Pauli? No.
# For d9, there are 10 determinants.
# They form ONE term: 2D (10 states: M_L=-2..2, S=1/2).
# So Average Energy of Configuration MUST EQUAL Energy of 2D Term.
# 
# Unless... the "Average Energy" formula used in RHF
# includes approximations (like discarding off-diagonal L terms).
#
# Let's check "Spin-Polarized Hartree Fock" (UHF) vs RHF?
# Closed Shell RHF is exact for d10.
# The Ion is d9. d9 is Open Shell.
# RHF for Open Shell (ROHF) vs "Average Energy HF".
# 
# If Clementi-Roetti uses ROHF (Restricted Open Shell) for the ion,
# but uses "Average Energy HF" for the orbital energy definition?
#
# Let's just calculate the values.
# 
# Hypothesis: The 'Missing Factor' 2.25 applies to the Exchange term.
# My BARE exchange term was ~0.22 Ha.
# If I scale by 2.25, it becomes ~0.50 Ha.
# Total Energy diff = 0.28 Ha.
#
# Let's verify what 2/7 * F2 + 2/7 * F4 is.
Ex_avg = (2/7) * (F2 + F4) # Simplistic sum
print(f"My BARE Exchange K (approx): {Ex_avg:.6f} (This is roughly 0.2-0.3)")

