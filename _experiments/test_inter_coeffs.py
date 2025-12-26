# -*- coding: utf-8 -*-
"""Test with INTERMEDIATE coefficients: c_k = (2l_source+1) * (3j)^2"""
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

# INTERMEDIATE coefficients: c_k = (2l_source+1) * (3j)^2
EXCHANGE_COEFFS_INTER = {
    (0, 0): [(0, 1.0)],         # 1 * 1.0
    (0, 1): [(1, 1.0)],         # 3 * 1/3
    (1, 0): [(1, 1/3)],         # 1 * 1/3
    (1, 1): [(0, 1.0), (2, 0.4)],  # 3 * (1/3, 2/15)
    (0, 2): [(2, 1.0)],         # 5 * 1/5
    (2, 0): [(2, 0.2)],         # 1 * 1/5
    (1, 2): [(1, 2/3), (3, 3/7)],  # 5 * (2/15, 3/35)
    (2, 1): [(1, 0.4), (3, 9/35)],  # 3 * (2/15, 3/35)
    (2, 2): [(0, 1.0), (2, 2/7), (4, 2/7)],  # 5 * (1/5, 2/35, 2/35)
}

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
    Y = I_inner / safe_r**(k+1) + r**k * I_outer
    return Y

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

print("TESTING INTERMEDIATE COEFFICIENTS: c_k = (2*l_source+1) * (3j)^2")
print("=" * 70)

CONFIGS = {
    'He': {'1s': 2},
    'Be': {'1s': 2, '2s': 2},
    'Ne': {'1s': 2, '2s': 2, '2p': 6},
    'Ar': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6},
    'Zn': {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10},
}

for symbol, config in CONFIGS.items():
    atom_data = data[symbol]
    Z = atom_data['Z']
    
    orbitals = {}
    for orb_name, occ in config.items():
        terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
                 for t in atom_data['orbitals'][orb_name]]
        orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
        R = compute_orbital_radial(orb, r)
        P = r * R
        orbitals[orb_name] = (orb, R, P, occ)
    
    print(f"\n{symbol} (Z={Z}):")
    
    for target_name in config.keys():
        target, R_target, P_target, n_target = orbitals[target_name]
        l_target = target.l
        
        T = compute_kinetic_energy(target, R_target, r)
        safe_r = np.where(r > 1e-15, r, 1e-15)
        V_nuc = np.trapezoid(P_target**2 * (-Z/safe_r), r)
        
        # Coulomb: n_j * F^0
        E_coulomb = 0
        for source_name, (source, R_j, P_j, n_j) in orbitals.items():
            Y0 = compute_Yk(P_j, P_j, r, k=0)
            F0 = np.trapezoid(P_target**2 * Y0, r)
            E_coulomb += n_j * F0
        
        # Exchange with INTERMEDIATE coefficients: (n_j/2) * Σ_k c_k * G^k
        # But note: n_j/2 already has the (2l_source+1) factor
        # So we should use: (n_j / (2*(2l_source+1))) * Σ_k c_k_INTER * G^k
        # = (n_j / (2*(2l_source+1))) * (2l_source+1) * Σ_k c_k_BARE * G^k
        # = (n_j / 2) * Σ_k c_k_BARE * G^k
        # This is my original formula with BARE coefficients!
        
        E_exchange = 0
        for source_name, (source, R_j, P_j, n_j) in orbitals.items():
            l_source = source.l
            # Use (n_j / (2*(2l_source+1))) with INTER coefficients = n_j/2 with BARE
            # Let's try direct: (2l_source+1) spatial orbitals, each contributes K_avg
            # K_total = (2l_source+1) * K_avg = (2l_source+1) * Σ_k c_k_INTER * G^k / (2l_source+1)
            # Wait, this doesn't make sense...
            
            # Let me use: each source orbital pair contributes K_avg
            # Number of source orbitals = (2l_source+1)
            # K_contribution = (2l_source+1) * Σ_k [c_k / (2l_source+1)] * G^k = Σ_k c_k * G^k
            # Hmm, this gives FULL coefficients...
            
            # Actually: K = (# same-spin source electrons) * K_avg_per_electron
            # = (n_j/2) * K_avg
            # K_avg = Σ_k c_k_INTER * G^k  (INTER includes (2l_source+1))
            # So K = (n_j/2) * Σ_k (2l_source+1)*(3j)^2 * G^k
            
            coeffs = EXCHANGE_COEFFS_INTER.get((l_target, l_source), [(0, 1.0)])
            for k, c_k in coeffs:
                Yk = compute_Yk(P_target, P_j, r, k=k)
                Gk = np.trapezoid(P_target * P_j * Yk, r)
                # Divide by (2l_source+1) to get per-electron contribution
                E_exchange += (n_j / 2) * c_k / (2*l_source+1) * Gk
        
        eps_calc = T + V_nuc + E_coulomb - E_exchange
        eps_ref = atom_data['energies'].get(target_name, 0)
        error = 100 * (eps_calc - eps_ref) / abs(eps_ref) if eps_ref != 0 else 0
        
        status = "✓" if abs(error) < 1 else ("~" if abs(error) < 5 else "✗")
        print(f"  {target_name}: calc={eps_calc:.6f}, ref={eps_ref:.6f}, err={error:.2f}% {status}")
