# -*- coding: utf-8 -*-
"""
Calculate the EXACT sum of exchange coefficients for a SINGLE electron
interacting with a closed shell.

We want to calculate:
Sum_{m'} [ c^k(lm, lm') ]^2

where c^k is the coefficient such that:
Exchange_Integr(lm, lm') = [c^k]^2 * F^k

This will tell us the exact K potential felt by one electron.
"""
import numpy as np
from sympy.physics.wigner import wigner_3j

def get_ck_sq(l, m, lp, mp, k):
    # c^k coefficient for Exchange K:
    # 3j(l k lp; -m (m-mp) mp) * (...)
    # Standard formula for coefficient of G^k/F^k in exchange:
    # (2l+1)(2lp+1) * [3j(l k lp; 0 0 0)]^2 * [3j(l k lp; -m q mp)]^2 * (-1)...?
    
    # Let's use the explicit Gaunt coefficient formula relation.
    # The exchange integral K_{ij} coeff for F^k is:
    # (2l+1)^2 * [3j(l k l; 0 0 0)]^2 * [3j(l k l; -m (m-mp) mp)]^2
    # Note: No phase factor needed as we square it.
    
    prefactor = (2*l+1)*(2*lp+1) * float(wigner_3j(l, k, lp, 0, 0, 0))**2
    geom = float(wigner_3j(l, k, lp, -m, m-mp, mp))**2
    return prefactor * geom

def calculate_potential_sum(l_target, l_source):
    print(f"Calculating Potential Sum for l={l_target} <-> l'={l_source}")
    
    # Average over all target m
    total_avg_k = {}
    
    for m in range(-l_target, l_target+1):
        # Sum over all source m' (full shell)
        sum_k = {}
        for mp in range(-l_source, l_source+1):
            # Same spin exchange only
            # Check allowed k
            # |l-l'| <= k <= l+l'
            # Parity: l+l'+k even
            k_min = abs(l_target - l_source)
            k_max = l_target + l_source
            for k in range(k_min, k_max+1, 2):
                val = get_ck_sq(l_target, m, l_source, mp, k)
                sum_k[k] = sum_k.get(k, 0) + val
        
        # Add to total average
        for k, v in sum_k.items():
            total_avg_k[k] = total_avg_k.get(k, 0) + v
            
    # Divide by (2l_target+1) to get average over m
    deg = 2*l_target + 1
    print(f"Average Potential Coefficients (per electron):")
    for k in sorted(total_avg_k.keys()):
        avg = total_avg_k[k] / deg
        print(f"  k={k}: {avg:.6f}")
        
    return total_avg_k

print("--- d-d Intra-shell (l=2, l'=2) ---")
avg_dd = calculate_potential_sum(2, 2)

print("\n--- p-p Intra-shell (l=1, l'=1) ---")
avg_pp = calculate_potential_sum(1, 1)

print("\n--- s-d Inter-shell (l=0, l'=2) ---")
avg_sd = calculate_potential_sum(0, 2)

print("\n--- p-d Inter-shell (l=1, l'=2) ---")
avg_pd = calculate_potential_sum(1, 2)


print("\n--- Check My BARE Coeffs for comparison ---")
# My BARE coeffs are essentially (2l+1)(2l'+1) [3j]^2 ? NO.
# My BARE code uses: c_k = (2l+1)(2l'+1) [3j(l k l'; 0 0 0)]^2 ??? 
# Wait, let's verify what my BARE code (Step 662) thinks BARE is.
# My code: c_k = (3j(l k l'; 0 0 0))^2  (just the 3j squared)
# And then multiplied by (n_j/2).

def check_bare(l, lp, k):
    w = float(wigner_3j(l, k, lp, 0, 0, 0))**2
    print(f"BARE (3j^2) for l={l}, l'={lp}, k={k}: {w:.6f}")
    
check_bare(2, 2, 2)
check_bare(2, 2, 4)
