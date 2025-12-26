# -*- coding: utf-8 -*-
"""
Calculate correct Slater-Condon angular coefficients using SymPy
================================================================

The c_k coefficients for exchange integrals are derived from 
Gaunt coefficients (integrals of three spherical harmonics).
"""

from sympy.physics.wigner import gaunt, wigner_3j
from sympy import sqrt, pi, simplify, Rational, S
import numpy as np

def compute_ck_coefficient(l1, l2, k):
    """
    Compute the c_k angular coefficient for exchange integral.
    
    The exchange integral involves:
    K = Σ_k c_k(l1, l2) * R_k(l1, l2)
    
    where c_k is the angular factor and R_k is the radial Slater integral.
    
    The angular factor is:
    c_k = (2l1+1)(2l2+1) * [Wigner 3j symbol]^2 * (some factor)
    
    For a filled subshell, we sum over all m values.
    """
    # For exchange between orbitals with angular momenta l1 and l2,
    # the c_k coefficient is given by the sum of squared Gaunt coefficients
    
    # Selection rules: k must satisfy |l1-l2| <= k <= l1+l2 and l1+l2+k even
    if k < abs(l1 - l2) or k > l1 + l2:
        return 0
    if (l1 + l2 + k) % 2 != 0:
        return 0
    
    # The coefficient is related to Wigner 3j symbols:
    # c_k = (4*pi / (2k+1)) * sum_{m1, m2} |<l1 m1 | k (m1-m2) | l2 m2>|^2
    # 
    # For a filled shell (all m values occupied), the formula simplifies to:
    # c_k = (2l1+1)(2l2+1) / (2k+1) * |<l1 0 k 0 | l2 0>|^2
    # where the 3j symbol is <l1 0 k 0 | l2 0>
    
    w3j = wigner_3j(l1, k, l2, 0, 0, 0)
    result = (2*l1 + 1) * (2*l2 + 1) * w3j**2
    
    return simplify(result)


def compute_ck_slater(l1, l2, k):
    """
    Compute c_k using the standard Slater formulation.
    
    From Slater's "Quantum Theory of Atomic Structure" (1960):
    c^k(l, l') = (2l+1)(2l'+1) * [3j(l k l'; 0 0 0)]^2
    
    This is the factor that multiplies the radial Slater integral R^k.
    """
    if k < abs(l1 - l2) or k > l1 + l2:
        return S(0)
    if (l1 + l2 + k) % 2 != 0:
        return S(0)
    
    w3j = wigner_3j(l1, k, l2, 0, 0, 0)
    c_k = (2*l1 + 1) * (2*l2 + 1) * w3j**2
    
    return simplify(c_k)


def main():
    print("=" * 60)
    print("SLATER-CONDON ANGULAR COEFFICIENTS (c_k)")
    print("From Wigner 3j Symbols")
    print("=" * 60)
    
    # Define orbital types
    interactions = [
        ("s-s", 0, 0),
        ("s-p", 0, 1),
        ("p-s", 1, 0),
        ("p-p", 1, 1),
        ("s-d", 0, 2),
        ("d-s", 2, 0),
        ("p-d", 1, 2),
        ("d-p", 2, 1),
        ("d-d", 2, 2),
    ]
    
    print(f"\n{'Interaction':>10} {'k':>3} {'c_k (symbolic)':>25} {'c_k (decimal)':>15}")
    print("-" * 60)
    
    for name, l1, l2 in interactions:
        k_min = abs(l1 - l2)
        k_max = l1 + l2
        
        for k in range(k_min, k_max + 1):
            c_k = compute_ck_slater(l1, l2, k)
            if c_k != 0:
                try:
                    c_k_float = float(c_k)
                except:
                    c_k_float = float(c_k.evalf())
                print(f"{name:>10} {k:>3} {str(c_k):>25} {c_k_float:>15.6f}")
    
    print("\n" + "=" * 60)
    print("COMPARISON WITH MY CURRENT VALUES")
    print("=" * 60)
    
    my_values = {
        (0, 0): {0: 1.0},
        (0, 1): {1: 1/3},
        (1, 1): {0: 1.0, 2: 2/5},
        (0, 2): {2: 1/5},
        (1, 2): {1: 2/15, 3: 3/35},
        (2, 2): {0: 1.0, 2: 2/7, 4: 2/7},
    }
    
    print(f"\n{'Interaction':>10} {'k':>3} {'My value':>12} {'Correct':>12} {'Match':>8}")
    print("-" * 50)
    
    for (l1, l2), ks in my_values.items():
        name = f"l={l1}-l={l2}"
        for k, my_val in ks.items():
            c_k = compute_ck_slater(l1, l2, k)
            try:
                correct = float(c_k)
            except:
                correct = float(c_k.evalf())
            
            match = "✓" if abs(my_val - correct) < 0.001 else "✗"
            print(f"{name:>10} {k:>3} {my_val:>12.6f} {correct:>12.6f} {match:>8}")


if __name__ == "__main__":
    main()
