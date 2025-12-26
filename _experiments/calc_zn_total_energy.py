# -*- coding: utf-8 -*-
"""
Calculate TOTAL Energy of Zn using BARE coefficients.
Goal: Check if the BARE coefficients reproduce the Total Energy (-1777.8 Ha).
If they do, the orbital energy mismatch is a property of the RHF method (Koopmans' deviation),
not a bug in my code.
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

def get_n(name):
    import re
    match = re.search(r'(\d+)', name)
    if match: return int(match.group(1))
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
    integrand_inner = rho * r**k
    I_inner = np.zeros(N)
    I_inner[1:] = cumulative_trapezoid(integrand_inner, r)
    integrand_outer = rho / safe_r**(k+1)
    total_outer = np.trapezoid(integrand_outer, r)
    I_outer = np.zeros(N)
    I_outer[0] = total_outer
    I_outer[1:] = total_outer - cumulative_trapezoid(integrand_outer, r)
    return I_inner / safe_r**(k+1) + r**k * I_outer

# Analytic Kinetic Energy for STO
# T = <R| -1/2 nabla^2 |R>
# For STO R_nl = N * r^(n-1) * exp(-zeta*r)
# Kinetic Energy integral for normalized STO:
# T = zeta^2 / 2  (for 1s? No, general formula needed)
# General formula for STO(n, zeta):
# T = zeta^2 / 2 * [ 1 - n*(n-1)/ ( (n+l)*(n+l+1/2)? No... ) ]
# Easier to calculate overlap of Laplacian.
# Or better: Use the property that T = zeta^2/2 - Z_eff * zeta ?? No.

def compute_kinetic_energy_analytic_term(t1: STOTerm, t2: STOTerm, l: int):
    # Calculate kinetic energy integral between two STO terms: <t1| -1/2 \nabla^2 |t2>
    # Overlap integral formula:
    # <n1,z1 | n2,z2> = (2*sqrt(z1*z2))^(n1+0.5) * (z1+z2)^-(n1+n2+1) * factorial(n1+n2) ...
    #
    # Analytic Kinetic Energy Matrix Element T_ij:
    # T_ij = -0.5 * zeta_j^2 * S(n_i, z_i; n_j, z_j) 
    #        + zeta_j * n_j * S(n_i, z_i; n_j-1, z_j) ... (This is using the derivative of STO)
    
    # Let's construct it properly.
    # (-)1/2 * nabla^2 (r^(n-1) e^(-zr) Y_lm)
    # Radial part: 
    # (-1/2) * [ R'' + (2/r)R' - l(l+1)/r^2 R ]
    
    # Applying operator to R_j = r^(n_j-1) e^(-z_j r):
    # A = -0.5
    # Term 1: -z_j^2 * r^(n_j-1) e^(-z_j r)  (from second deriv of exp)
    # Term 2: 2 * z_j * (n_j-1) * r^(n_j-2) e^(-z_j r) (cross term)
    # Term 3: - (n_j-1)(n_j-2) * r^(n_j-3) e^(-z_j r) (second deriv of poly)
    # ... This is getting complex to implement inline.
    
    # Alternative: Use simple 1D integration of the EXACT analytic derivative expression.
    # This avoids np.gradient errors but keeps the loop structure.
    pass

def compute_kinetic_energy_analytic_integral(orbital, r):
    # Calculates T by integrating the exact analytic Laplacian over the grid.
    # This avoids finite difference errors.
    
    # 1. Evaluate Laplacian analytically
    l = orbital.l
    Laplacian_R = np.zeros_like(r)
    
    for term in orbital.terms:
        n = term.n
        z = term.zeta
        N = sto_normalization(n, z)
        
        # R(r) = N * r^(n-1) * exp(-z*r)
        # dR/dr = N * [ (n-1)r^(n-2) - z*r^(n-1) ] * exp(-z*r)
        # d^2R/dr^2 = N * [ (n-1)(n-2)r^(n-3) - 2z(n-1)r^(n-2) + z^2 r^(n-1) ] * exp(-z*r)
        
        # Laplacian_radial = d^2R/dr^2 + (2/r)dR/dr - l(l+1)/r^2 R
        
        # Component 1: d^2R/dr^2
        d2R = (n-1)*(n-2) * np.power(r, n-3) - 2*z*(n-1)*np.power(r, n-2) + z**2 * np.power(r, n-1)
        
        # Component 2: (2/r)dR/dr
        # (2/r) * [ (n-1)r^(n-2) - z*r^(n-1) ] = 2(n-1)r^(n-3) - 2z r^(n-2)
        r_inv_dR = 2*(n-1)*np.power(r, n-3) - 2*z*np.power(r, n-2)
        
        # Component 3: -l(l+1)/r^2 * R
        # -l(l+1) * r^(n-3)
        centrifugal = -l*(l+1) * np.power(r, n-3)
        
        # Combine
        # Power n-3 terms: (n-1)(n-2) + 2(n-1) - l(l+1) = n^2 - 3n + 2 + 2n - 2 - l(l+1) = n(n-1) - l(l+1)
        # Power n-2 terms: -2z(n-1) - 2z = -2zn
        # Power n-1 terms: z^2
        
        term_laplacian = N * np.exp(-z*r) * (
            (n*(n-1) - l*(l+1)) * np.power(r, n-3)
            - 2*z*n * np.power(r, n-2)
            + z**2 * np.power(r, n-1)
        )
        
        Laplacian_R += term.coeff * term_laplacian
        
    # 2. Integrate -0.5 * R * Laplacian_R * r^2
    # Re-evaluate R (we already have R vector passed in, but safe to recompute or pass)
    # Actually R is passed in.
    
    # NOTE: r^(n-3) * r^2 = r^(n-1). 
    # For n=1 (1s), n-3 = -2. r^(n-1) = r^0 = 1.
    # We need to handle r->0 carefully.
    # Use the grid r.
    
    integrand = -0.5 * orbital_R(orbital, r) * Laplacian_R * r**2
    return np.trapezoid(integrand, r)

def orbital_R(orbital, r):
    # Helper to reconstruct R if needed
    R = np.zeros_like(r)
    for term in orbital.terms:
        N = sto_normalization(term.n, term.zeta)
        R += term.coeff * N * np.power(r, term.n - 1) * np.exp(-term.zeta * r)
    return R

def compute_kinetic_energy(orbital, R, r):
    # Wrapper to switch to semi-analytic method
    return compute_kinetic_energy_analytic_integral(orbital, r)

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
r = create_log_grid(r_min=1e-9, r_max=100.0, N=50000)

symbol = 'Zn'
atom_data = data[symbol]
Z = atom_data['Z']
config = {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10}

print(f"CALCULATING TOTAL ENERGY FOR {symbol} (Z={Z})")
print("Using Standard BARE coefficients (No Scaling)")
print("-" * 60)

orbitals = {}
for orb_name, occ in config.items():
    terms = [STOTerm(n=t['nStar'], zeta=t['zeta'], coeff=t['coeff']) 
             for t in atom_data['orbitals'][orb_name]]
    orb = Orbital(name=orb_name, l=get_l(orb_name), terms=terms)
    R = compute_orbital_radial(orb, r)
    P = r * R
    orbitals[orb_name] = (orb, R, P, occ)

# 1. Total Kinetic Energy
T_total = 0
for name, (orb, R, P, occ) in orbitals.items():
    t_orb = compute_kinetic_energy(orb, R, r)
    T_total += occ * t_orb
print(f"Total Kinetic Energy (T): {T_total:.6f} Ha")

# 2. Total Nuclear Energy
V_nuc_total = 0
for name, (orb, R, P, occ) in orbitals.items():
    safe_r = np.where(r > 1e-15, r, 1e-15)
    v_orb = np.trapezoid(P**2 * (-Z/safe_r), r)
    V_nuc_total += occ * v_orb
print(f"Total Nuclear Potential (V_nuc): {V_nuc_total:.6f} Ha")

# 3. Total Electron-Electron Energy (J - K)
# E_ee = 0.5 * Sum_{ij} (J_{ij} - K_{ij})
# Sum over all electron pairs.
# My previous orbital loops calculated: epsilon_i = T_i + V_nuc_i + Sum_j (J_ij - K_ij)
# Total E = Sum_i epsilon_i - E_ee (double counting correction)
# Or just Sum_i (T_i + V_nuc_i) + E_ee

J_total = 0
K_total = 0

orb_names = list(orbitals.keys())
for i in range(len(orb_names)):
    name_i = orb_names[i]
    orb_i, R_i, P_i, occ_i = orbitals[name_i]
    l_i = orb_i.l
    
    for j in range(len(orb_names)):
        name_j = orb_names[j]
        orb_j, R_j, P_j, occ_j = orbitals[name_j]
        l_j = orb_j.l
        
        # Coulomb J_ij
        # Interaction between SHELL i and SHELL j
        # E_J = 0.5 * occ_i * occ_j * F0(i, j)
        Y0 = compute_Yk(P_j, P_j, r, k=0)
        F0 = np.trapezoid(P_i**2 * Y0, r)
        J_total += 0.5 * occ_i * occ_j * F0
        
        # Exchange K_ij
        # E_K = 0.5 * occ_i * occ_j * Sum_k c_k G^k ?
        # My formula in orbital calc is: V_ex_i = Sum_j (occ_j/2) * coeffs * Gk
        # E_ex_tot = 0.5 * Sum_i occ_i * V_ex_i
        #          = 0.5 * Sum_{ij} occ_i * (occ_j/2) * coeffs * Gk
        #          = 0.25 * occ_i * occ_j * coeffs * Gk
        
        # Wait, the pair factor.
        # Total Exchange Energy between Shell I and Shell J.
        # E_ex(I,J) = - Sum_{pairs} K_{uv}
        # If I=J (Intra-shell):
        #   E_ex = - Sum_{pairs} K (n(n-1)/2 pairs? No, spin)
        
        # Let's use the explicit Average Energy Coefficient formula from derive_coeffs.py results
        # To be rigorous, we should use the Exact Fk/Gk coefficients for Total Energy
        # But here I want to test MY CODE's BARE coefficients.
        
        # My code uses:
        # V_ex term in Fock = Sum_j (n_j/2) * c_k * G^k
        # Total E contribution = 1/2 * Sum_i n_i * <i|V_ex|i>
        # = 1/2 * Sum_i n_i * [ Sum_j (n_j/2) * c_k * G^k(ij) ]
        # = 1/4 * Sum_{ij} n_i n_j * c_k * G^k
        
        coeffs = EXCHANGE_COEFFS_BARE.get((l_i, l_j), [(0, 1.0)])
        
        # FORCE BARE SCALING (1.0) for Validation
        scale = 1.0 
        # (Removed Heuristic Scaling Logic)
             
        for k, c_k in coeffs:
            # Do not scale F0/G0 if any
            current_scale = scale
            if k == 0 and l_i == l_j:
                current_scale = 1.0

            Yk = compute_Yk(P_i, P_j, r, k=k)
            Gk = np.trapezoid(P_i * P_j * Yk, r)
            K_total += 0.25 * occ_i * occ_j * current_scale * c_k * Gk

def get_n(name):
    import re
    match = re.search(r'(\d+)', name)
    if match: return int(match.group(1))
    return 0


print(f"Total Coulomb (J): {J_total:.6f} Ha")
print(f"Total Exchange (K) [BARE]: {K_total:.6f} Ha")
E_ee = J_total - K_total
print(f"Total Electron-Electron (V_ee): {E_ee:.6f} Ha")

E_total = T_total + V_nuc_total + E_ee
print(f"Total Energy: {E_total:.6f} Ha")
print(f"Reference Total E: -1777.8481 Ha")
print(f"Difference: {E_total - (-1777.8481):.6f} Ha")

print(f"Virial Ratio (-V/T): {-(V_nuc_total + E_ee)/T_total:.6f} (Should be 2.0)")

