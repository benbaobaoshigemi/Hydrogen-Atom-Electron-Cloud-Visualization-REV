# -*- coding: utf-8 -*-
"""
Derive EXACT Average Energy of Configuration Coefficients f_k and g_k.
Method:
1. Generate all microstates (determinants) for a given configuration (l^N).
2. Calculate the exchange energy contribution to the Total Energy for each determinant.
   E_ex(D) = - Sum_{i<j in D, ms_i=ms_j} K_{ij}
   where K_{ij} = Sum_k c^k(lm_i, lm_j) * F^k
3. Average E_ex over all determinants.
4. Extract the coefficient of each F^k term.
   E_avg_ex = - Sum_k (Coeff_k * F^k)

This provides the rigorous scalar coefficients for the "Average Energy of Configuration".
"""

import numpy as np
from sympy.physics.wigner import wigner_3j
from sympy import simplify, Rational
from itertools import combinations

def get_ck_coeff(l1, m1, l2, m2, k):
    """
    Calculate c^k(lm, l'm') coefficient for raw interaction.
    c^k = (2l+1)(2l'+1) [3j(...)]^2 ... Wait, the standard formula is:
    <lm, l'm' | 1/r12 | lm, l'm'> = Sum_k c^k(lm, l'm') F^k
    c^k(lm, l'm') = (-1)^(m+m') * (2l+1)(2l'+1) * 
                    ( (l 0 l' 0 | k 0) * (l -m l' m' | k m'-m) ) ... No.
    
    Let's use the definition:
    <lm, l'm' | g | lm, l'm'> = Sum_k c^k(lm,l'm') F^k(nl, n'l')
    c^k = (2l+1)(2l'+1) * [3j(l k l'; 0 0 0)]^2 * (-1)^(m-m') ...?
    
    Actually, use Gaunt coefficients.
    c^k(lm, l'm') = coeff of F^k in Coulomb expansion.
    """
    # The interaction energy is:
    # J(i,j) = Sum_k a^k(lm_i, lm_j) F^k
    # K(i,j) = Sum_k b^k(lm_i, lm_j) G^k
    
    # a^k(lm, l'm') = c^k(lm, lm) * c^k(l'm', l'm')
    # b^k(lm, l'm') = [c^k(lm, l'm')]^2
    
    # where c^k(lm, l'm') is the coefficient in the expansion of spherical harmonics.
    # c^k(lm, l'm') = sqrt((2l+1)(2l'+1))/(2k+1) * <lm|Yk|l'm'>... NO.
    
    # Using 3j symbols directly:
    # <l m | Y_kq | l' m'> = (-1)^m * sqrt((2l+1)(2k+1)(2l'+1)/(4pi)) * 3j(l k l'; -m q m') * 3j(l k l'; 0 0 0)
    # But we define the radial integral F^k to absorb the 4pi denominators?
    # Standard definition:
    # coeff of F^k in Direct J_ij: c^k(l m, l m) * c^k(l' m', l' m')
    # coeff of G^k in Exchange K_ij: [c^k(l m, l' m')]^2
    
    # The c^k used here is:
    # c^k(l m, l' m') = sqrt((2l+1)(2l'+1)) * (-1)^m * 3j(l k l'; -m (m-m') m') * 3j(l k l'; 0 0 0)
    
    w3j_0 = wigner_3j(l1, k, l2, 0, 0, 0)
    if w3j_0 == 0: return 0
    
    q = m1 - m2
    w3j_m = wigner_3j(l1, k, l2, -m1, q, m2)
    
    # Note: (-1)^m phase
    phase = (-1)**abs(m1) 
    
    # Is it exactly this? 
    # Let's rely on the property: sum_m c^k c^k = (2l+1) 
    
    val = phase * np.sqrt((2*l1+1)*(2*l2+1)) * float(w3j_0) * float(w3j_m)
    return val

def calculate_average_exchange(l, N_electrons, config_name):
    print(f"--- Configuration {config_name} ({l}^{N_electrons}) ---")
    
    # 1. Generate spin-orbitals
    # (m, s) where m in -l..l, s in -0.5, 0.5
    orbitals = []
    for m in range(-l, l+1):
        orbitals.append( (m, 0.5) )
        orbitals.append( (m, -0.5) )
    
    # 2. Generate Determinants
    dets = list(combinations(orbitals, N_electrons))
    n_dets = len(dets)
    print(f"Number of determinants: {n_dets}")
    
    # 3. Sum Exchange Energy coeffs
    # We want final form: E_ex_avg = - Sum_k (Coeff_k * F^k)
    # Note: For equivalent electrons, J and K use the SAME parameters F^k.
    # Total interaction = Sum_{i<j} [ J_ij - delta(ms_i, ms_j) K_ij ]
    # J_ij comes from F^k (k=0, 2, 4...)
    # K_ij comes from F^k (k=0, 2, 4...)
    
    total_coeffs = {} # k -> sum_value
    
    for det in dets:
        # For each pair in the determinant
        for i in range(len(det)):
            for j in range(i+1, len(det)):
                m1, s1 = det[i]
                m2, s2 = det[j]
                
                # Direct Term (J)
                # J_12 = Sum_k c^k(1,1)c^k(2,2) F^k
                for k in range(0, 2*l+1, 2):
                    c1 = get_ck_coeff(l, m1, l, m1, k)
                    c2 = get_ck_coeff(l, m2, l, m2, k)
                    # We are only interested in Exchange average, but actually
                    # Average Energy includes J and K. The Coeff f_k usually combines both.
                    # But if we JUST want exchange contribution to check my K formula:
                    # Let's track J and K separately.
                    pass
                
                # Exchange Term (K) - only if parallel spins
                if s1 == s2:
                    for k in range(0, 2*l+1, 2): # k must be even for equivalent electrons? No, standard parity.
                        # Wait, for equivalent electrons, K_ij involves F^k.
                        # The coefficient is [c^k(l m1, l m2)]^2
                        ck = get_ck_coeff(l, m1, l, m2, k)
                        val = ck**2
                        
                        total_coeffs[k] = total_coeffs.get(k, 0) + val

    # 4. Average
    print("Average Exchange Coefficients (coeff of -F^k):")
    for k in sorted(total_coeffs.keys()):
        avg = total_coeffs[k] / n_dets
        print(f"  k={k}: {avg:.6f}")
        
    return total_coeffs

# Run for d^2, d^5, d^10
calculate_average_exchange(2, 2, "d^2")
calculate_average_exchange(2, 5, "d^5")
calculate_average_exchange(2, 10, "d^10")

# Run for p^2, p^6
calculate_average_exchange(1, 2, "p^2")
calculate_average_exchange(1, 6, "p^6")

def calculate_intershell_exchange(l1, N1, l2, N2, config_name):
    print(f"--- Configuration {config_name} ({l1}^{N1} {l2}^{N2}) ---")
    
    # 1. Generate spin-orbitals for subshell 1
    orb1 = []
    for m in range(-l1, l1+1):
        orb1.append((m, 0.5, 1))
        orb1.append((m, -0.5, 1))
        
    # 2. Generate spin-orbitals for subshell 2
    orb2 = []
    for m in range(-l2, l2+1):
        orb2.append((m, 0.5, 2))
        orb2.append((m, -0.5, 2))
        
    # 3. Generate Determinants
    # We select N1 from orb1 and N2 from orb2
    from itertools import product
    dets1 = list(combinations(orb1, N1))
    dets2 = list(combinations(orb2, N2))
    
    # Combined determinants
    dets = list(product(dets1, dets2))
    n_dets = len(dets)
    print(f"Number of determinants: {n_dets}")
    
    # 4. Calculate Average Exchange between Shell 1 and Shell 2
    # We only care about K_ij where i in Shell 1, j in Shell 2
    
    total_coeffs = {} # k -> sum matrix element
    
    for d1, d2 in dets:
        # Full list of occupied spin-orbitals
        occupied = d1 + d2
        
        # Cross-terms only
        for i in range(len(d1)):
            for j in range(len(d2)):
                m1, s1, shell1 = d1[i]
                m2, s2, shell2 = d2[j]
                
                # Exchange only if parallel spins
                if s1 == s2:
                    # G^k coefficient: [c^k(l1 m1, l2 m2)]^2
                    for k in range(abs(l1-l2), l1+l2+1): 
                        # parity check? l1+l2+k must be even? 
                        # No, for G^k (exchange), parity condition is l1+l2+k is even.
                        if (l1 + l2 + k) % 2 == 0:
                            ck = get_ck_coeff(l1, m1, l2, m2, k)
                            val = ck**2
                            total_coeffs[k] = total_coeffs.get(k, 0) + val
                            
    # 5. Average
    print("Average Inter-shell Exchange Coefficients (coeff of -G^k):")
    for k in sorted(total_coeffs.keys()):
        avg = total_coeffs[k] / n_dets
        print(f"  k={k}: {avg:.6f}")

# Check s-d exchange
calculate_intershell_exchange(0, 2, 2, 10, "s^2 d^10 (Full)")
# Check p-d exchange
calculate_intershell_exchange(1, 6, 2, 10, "p^6 d^10 (Full)")

# KOOPMANS' THEOREM CHECK
# E(d10) - E(d9)
calculate_average_exchange(2, 9, "d^9")


