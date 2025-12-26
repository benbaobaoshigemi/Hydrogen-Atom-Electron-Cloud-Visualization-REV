# -*- coding: utf-8 -*-
"""
CORRECTED Orbital Energy Validation
====================================

Using the textbook RHF formula:
ε_a = h_aa + Σ_b (2*J_ab - K_ab)

where the sum is over OCCUPIED SPATIAL ORBITALS b (not electrons!)
- 2*J because each spatial orbital has 2 electrons
- K because only same-spin electrons exchange (factor of 1, not 2)

For subshells with degeneracy (p, d):
- We treat each subshell as ONE orbital for the purpose of this formula
- But the degeneracy factor appears in the angular coefficients
"""

import numpy as np
import json
import re
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

def create_log_grid(r_min=1e-6, r_max=50.0, N=2000):
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

# FULL coefficients for exchange (includes 2l'+1 factor for source orbital)
# K_ab = Σ_k c_k * R^k(a,b)
# For same l: c_k = (2l+1) * (3j(l,k,l;0,0,0))^2
EXCHANGE_COEFFS_FULL = {
    (0, 0): [(0, 1.0)],
    (0, 1): [(1, 1.0)],   # (2*1+1) * (1/3) = 1
    (1, 0): [(1, 1.0)],
    (1, 1): [(0, 1.0), (2, 0.4)],  # 3 * (1/3, 2/15) = (1, 2/5)
    (0, 2): [(2, 1.0)],   # (2*2+1) * (1/5) = 1
    (2, 0): [(2, 1.0)],
    (1, 2): [(1, 2/3), (3, 3/7)],  # 5 * (2/15, 3/35)
    (2, 1): [(1, 2/3), (3, 3/7)],
    (2, 2): [(0, 1.0), (2, 2/7), (4, 2/7)],  # 5 * (1/5, 2/35, 2/35)
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
            if idx == 0: cumsum_outer[idx] = total_outer
            else: cumsum_outer[idx] = total_outer - np.trapz(integrand_outer[:idx+1], r[:idx+1])
        return cumsum_inner + r * cumsum_outer
    else:
        Y = np.zeros(N)
        for idx in range(N):
            r_val = r[idx]
            if r_val < 1e-15: continue
            inner_integrand = rho[:idx+1] * r[:idx+1] ** (k + 2)
            I_inner = np.trapz(inner_integrand, r[:idx+1]) if idx > 0 else 0.0
            if k <= 1: outer_integrand = rho[idx:] * r[idx:] ** (1 - k)
            else: outer_integrand = rho[idx:] / (r[idx:] ** (k - 1))
            I_outer = np.trapz(outer_integrand, r[idx:]) if len(outer_integrand) > 1 else 0.0
            Y[idx] = I_inner / (r_val ** k) + (r_val ** (k + 1)) * I_outer
        return Y

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
        chi_double_prime = N * ((n - 1) * (n - 2) * r_nm3 - 2 * zeta * (n - 1) * r_nm2 + zeta ** 2 * r_nm1) * exp_term
        R_val += c * chi
        R_prime += c * chi_prime
        R_double_prime += c * chi_double_prime
    safe_r = np.where(r > 1e-15, r, 1e-15)
    laplacian = R_double_prime + (2 / safe_r) * R_prime - l * (l + 1) / (safe_r ** 2) * R_val
    return R_val * (-0.5) * laplacian

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

# Test on He first
data = parse_slater_basis(Path('../slater_basis.js'))
r = create_log_grid(r_min=1e-6, r_max=50.0, N=3000)

print("CORRECTED FORMULA TEST")
print("=" * 60)
print("Formula: ε_a = h_aa + Σ_b (2*J_ab - K_ab)")
print()

for symbol in ['He', 'Ne', 'Ar']:
    atom_data = data[symbol]
    Z = atom_data['Z']
    
    if symbol == 'He':
        config = {'1s': 2}
    elif symbol == 'Ne':
        config = {'1s': 2, '2s': 2, '2p': 6}
    elif symbol == 'Ar':
        config = {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6}
    
    orbitals = {}
    for orb_name, occ in config.items():
        terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
                 for t in atom_data['orbitals'][orb_name]]
        orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
        R = compute_orbital_radial(orb, r)
        orbitals[orb_name] = (orb, R, occ)
    
    print(f"{symbol} (Z={Z}):")
    
    for target_name in config.keys():
        target, R_target, n_target = orbitals[target_name]
        l_target = target.l
        
        # h_aa = T + V_nuc
        E_kin = compute_kinetic_density(target, R_target, r)
        T = np.trapz(E_kin * r**2, r)
        safe_r = np.where(r > 1e-15, r, 1e-15)
        V_nuc = np.trapz(R_target**2 * (-Z / safe_r) * r**2, r)
        h_aa = T + V_nuc
        
        # Sum over occupied orbitals
        two_J_minus_K = 0
        mask = r > 1e-10
        
        for source_name, (source, R_j, n_j) in orbitals.items():
            l_source = source.l
            
            # J_ab = ∫ R_a^2 * Y_0(b,b) / r * r^2 dr
            Y0 = compute_Yk(R_j, R_j, r, k=0)
            V_J = np.zeros_like(r)
            V_J[mask] = Y0[mask] / r[mask]
            J_ab = np.trapz(R_target**2 * V_J * r**2, r)
            
            # K_ab = Σ_k c_k * ∫ R_a * R_b * Y_k(a,b) / r * r^2 dr
            K_ab = 0
            coeffs = EXCHANGE_COEFFS_FULL.get((l_target, l_source), [(0, 1.0)])
            for k, c_k in coeffs:
                Yk = compute_Yk(R_target, R_j, r, k=k)
                V_k = np.zeros_like(r)
                V_k[mask] = Yk[mask] / r[mask]
                K_ab += c_k * np.trapz(R_target * R_j * V_k * r**2, r)
            
            # For doubly-occupied orbitals: 2J - K
            two_J_minus_K += 2 * J_ab - K_ab
        
        eps_calc = h_aa + two_J_minus_K
        eps_ref = atom_data['energies'].get(target_name, 0)
        error = 100 * (eps_calc - eps_ref) / abs(eps_ref) if eps_ref != 0 else 0
        
        status = "✓" if abs(error) < 1 else ("~" if abs(error) < 5 else "✗")
        print(f"  {target_name}: calc={eps_calc:.6f}, ref={eps_ref:.6f}, err={error:.2f}% {status}")
    
    print()
