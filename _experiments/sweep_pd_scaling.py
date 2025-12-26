# -*- coding: utf-8 -*-
"""
Sweep p-d inter-shell exchange scaling to minimize Kr 3d error.
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

# BARE coefficients 
EXCHANGE_COEFFS_BARE = {
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
    rho = P_i * P_j
    N = len(r)
    safe_r = np.where(r > 1e-20, r, 1e-20)
    # Use vectorized trapezoid if possible or cumulative
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
    P = r * R
    safe_r = np.where(r > 1e-15, r, 1e-15)
    
    # Robust numeric derivatives
    P_prime = np.gradient(P, r, edge_order=2)
    P_double_prime = np.gradient(P_prime, r, edge_order=2)
    
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

print("Sweeping p-d exchange scaling for Kr...")
symbol = 'Kr'
atom_data = data[symbol]
config = {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10, '4p': 6}
Z = atom_data['Z']

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
ref_eps = atom_data['energies'][target_name]

print(f"Target: Kr 3d (Ref: {ref_eps:.6f} Ha)")
print(f"Fixed: d-d intra-shell scaling = sqrt(5) = 2.236")
print("-" * 60)

scales = [1.0, 1.2, 1.4, 1.5, 1.732, 2.0, 2.236]

for pd_scale in scales:
    T = compute_kinetic_energy(target, R_target, r)
    safe_r = np.where(r > 1e-15, r, 1e-15)
    V_nuc = np.trapezoid(P_target**2 * (-Z/safe_r), r)
    
    E_coulomb = 0
    for source_name, (source, R_j, P_j, n_j) in orbitals.items():
        Y0 = compute_Yk(P_j, P_j, r, k=0)
        F0 = np.trapezoid(P_target**2 * Y0, r)
        E_coulomb += n_j * F0
    
    E_exchange = 0
    for source_name, (source, R_j, P_j, n_j) in orbitals.items():
        l_source = source.l
        coeffs = EXCHANGE_COEFFS_BARE.get((l_target, l_source), [(0, 1.0)])
        
        # LOGIC
        scale = 1.0
        # Fixed d-d scaling
        if l_source == l_target and l_source >= 2:
            scale = np.sqrt(2 * l_source + 1)
        
        # Test pd scaling
        if (l_target==1 and l_source==2) or (l_target==2 and l_source==1):
            scale = pd_scale
        
        for k, c_k in coeffs:
            current_scale = scale
            if k == 0 and l_source == l_target:
                current_scale = 1.0 
            
            Yk = compute_Yk(P_target, P_j, r, k=k)
            Gk = np.trapezoid(P_target * P_j * Yk, r)
            E_exchange += (n_j / 2) * current_scale * c_k * Gk

    eps_calc = T + V_nuc + E_coulomb - E_exchange
    error = 100 * (eps_calc - ref_eps) / abs(ref_eps)
    print(f"Scale p-d = {pd_scale:.3f} | Calc = {eps_calc:.6f} | Err = {error:.2f}%")

