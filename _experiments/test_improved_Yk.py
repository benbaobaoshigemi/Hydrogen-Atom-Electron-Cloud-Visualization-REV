# -*- coding: utf-8 -*-
"""
Test improved Y_k with scipy.integrate for better precision
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
    """Extended and finer grid for heavy atoms"""
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

def compute_Yk_improved(R_i, R_j, r, k=0):
    """Improved Y_k using scipy cumulative_trapezoid"""
    rho = R_i * R_j
    N = len(r)
    
    if k == 0:
        # Y_0(r) = (1/r) * ∫_0^r ρ*r'^2 dr' + ∫_r^∞ ρ*r' dr'
        integrand_inner = rho * r**2
        I_forward = np.zeros(N)
        I_forward[1:] = cumulative_trapezoid(integrand_inner, r)
        
        integrand_outer = rho * r
        I_backward = np.zeros(N)
        # Compute backward integral: ∫_r^∞ = ∫_0^∞ - ∫_0^r  
        total = np.trapezoid(integrand_outer, r)
        I_backward[0] = total
        I_backward[1:] = total - cumulative_trapezoid(integrand_outer, r)
        
        return I_forward + r * I_backward
    else:
        # General Y_k formula
        Y = np.zeros(N)
        
        # Forward integral: ∫_0^r ρ * r'^(k+2) dr'
        integrand_forward = rho * r**(k+2)
        I_forward = np.zeros(N)
        I_forward[1:] = cumulative_trapezoid(integrand_forward, r)
        
        # Backward integral: ∫_r^∞ ρ * r'^(1-k) dr'
        if k <= 1:
            integrand_back = rho * r**(1-k)
        else:
            safe_r = np.where(r > 1e-15, r, 1e-15)
            integrand_back = rho / safe_r**(k-1)
        
        total_back = np.trapezoid(integrand_back, r)
        I_backward = np.zeros(N)
        I_backward[0] = total_back
        I_backward[1:] = total_back - cumulative_trapezoid(integrand_back, r)
        
        # Y_k = I_forward / r^k + r^(k+1) * I_backward
        safe_r = np.where(r > 1e-15, r, 1e-15)
        Y = I_forward / safe_r**k + r**(k+1) * I_backward
        
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
        chi_double_prime = N * ((n - 1) * (n - 2) * r_nm3 - 2 * zeta * (n - 1) * r_nm2 + zeta**2 * r_nm1) * exp_term
        R_val += c * chi
        R_prime += c * chi_prime
        R_double_prime += c * chi_double_prime
    safe_r = np.where(r > 1e-15, r, 1e-15)
    laplacian = R_double_prime + (2/safe_r)*R_prime - l*(l+1)/(safe_r**2)*R_val
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

# Test
data = parse_slater_basis(Path('../slater_basis.js'))
r = create_log_grid(r_min=1e-7, r_max=100.0, N=5000)

print("IMPROVED PRECISION TEST (finer grid, scipy integration)")
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
        orbitals[orb_name] = (orb, R, occ)
    
    print(f"\n{symbol} (Z={Z}):")
    
    for target_name in config.keys():
        target, R_target, n_target = orbitals[target_name]
        l_target = target.l
        
        E_kin = compute_kinetic_density(target, R_target, r)
        T = np.trapezoid(E_kin * r**2, r)
        
        safe_r = np.where(r > 1e-15, r, 1e-15)
        V_nuc = np.trapezoid(R_target**2 * (-Z/safe_r) * r**2, r)
        
        E_hartree = 0
        E_exchange = 0
        mask = r > 1e-10
        
        for source_name, (source, R_j, n_j) in orbitals.items():
            Y0 = compute_Yk_improved(R_j, R_j, r, k=0)
            V = np.zeros_like(r)
            V[mask] = Y0[mask] / r[mask]
            E_hartree += np.trapezoid(R_target**2 * V * n_j * r**2, r)
            
            l_source = source.l
            coeffs = EXCHANGE_COEFFS.get((l_target, l_source), [(0, 1.0)])
            for k, c_k in coeffs:
                Yk = compute_Yk_improved(R_target, R_j, r, k=k)
                V_k = np.zeros_like(r)
                V_k[mask] = Yk[mask] / r[mask]
                E_exchange += c_k * np.trapezoid(R_target * R_j * V_k * n_j * 0.5 * r**2, r)
        
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
        print(f"  {s}: max error = {max(atom_errs):.2f}%")
