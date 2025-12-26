# -*- coding: utf-8 -*-
"""
HELIUM DEBUGGING SCRIPT
========================
Verify each component of the HF calculation for He.
"""

import math
from scipy import integrate

# He 1s orbital from Koga basis
# Consists of 5 STOs with coefficients
HE_1S_BASIS = [
    {'n': 2, 'zeta': 6.438865513242302, 'coeff': 0.0008103},
    {'n': 1, 'zeta': 3.385077039750975, 'coeff': 0.0798826},
    {'n': 1, 'zeta': 2.178370004614139, 'coeff': 0.180161},
    {'n': 1, 'zeta': 1.4553870053179185, 'coeff': 0.7407925},
    {'n': 2, 'zeta': 1.3552466748849417, 'coeff': 0.0272015},
]

def factorial(n):
    if n <= 1:
        return 1
    return n * factorial(n - 1)


def sto_radial(r, n, zeta):
    """STO radial function: R(r) = N * r^(n-1) * exp(-zeta*r)
    
    Note: The STO convention in Koga is:
    χ_n,ζ(r) = N * r^(n-1) * exp(-ζr)
    
    Normalization: ∫ |R|² r² dr = 1
    => N² ∫ r^(2n) exp(-2ζr) dr = 1
    => N² * (2n)! / (2ζ)^(2n+1) = 1
    => N = (2ζ)^(n+1/2) / sqrt((2n)!)
    
    Standard form: N = (2ζ)^n * sqrt(2ζ/(2n)!)
    """
    N = (2 * zeta) ** n * math.sqrt(2 * zeta / factorial(2 * n))
    return N * (r ** (n - 1)) * math.exp(-zeta * r)


def orbital_1s(r):
    """He 1s orbital as linear combination of STOs"""
    total = 0.0
    for term in HE_1S_BASIS:
        total += term['coeff'] * sto_radial(r, term['n'], term['zeta'])
    return total


def test_normalization():
    """Check that ∫ |ψ|² r² dr = 1"""
    def integrand(r):
        psi = orbital_1s(r)
        return psi * psi * r * r
    
    result, error = integrate.quad(integrand, 0, 50)
    print(f"Normalization: ∫|ψ|² r² dr = {result:.10f} (should be 1.0)")
    print(f"  Integration error: {error:.2e}")
    return result


def test_kinetic_numerical():
    """Compute <T> numerically"""
    def T_psi(r):
        # T = -0.5 * (d²/dr² + 2/r * d/dr - l(l+1)/r²)
        # For numerical differentiation
        h = 1e-5
        psi_m = orbital_1s(r - h) if r > h else 0
        psi_0 = orbital_1s(r)
        psi_p = orbital_1s(r + h)
        
        d2psi = (psi_p - 2*psi_0 + psi_m) / (h*h)
        dpsi = (psi_p - psi_m) / (2*h)
        
        T_local = -0.5 * (d2psi + 2/r * dpsi) if r > 0.001 else 0
        return T_local
    
    def integrand(r):
        psi = orbital_1s(r)
        return psi * T_psi(r) * r * r
    
    result, error = integrate.quad(integrand, 0.001, 50)
    print(f"Kinetic <T> (numerical): {result:.10f} Ha")
    return result


def test_nuclear_numerical():
    """Compute <V_nuc> = <-Z/r> numerically"""
    Z = 2
    
    def integrand(r):
        psi = orbital_1s(r)
        return psi * (-Z / r) * psi * r * r
    
    result, error = integrate.quad(integrand, 0.001, 50)
    print(f"Nuclear <V_nuc> (numerical): {result:.10f} Ha")
    return result


def test_coulomb_numerical():
    """Compute J = <ψψ|1/r12|ψψ> numerically (expensive!)"""
    
    def inner_integral(r1):
        """∫ |ψ(r2)|² / r_> r2² dr2"""
        def integrand(r2):
            psi2 = orbital_1s(r2)
            r_greater = max(r1, r2)
            return psi2 * psi2 * r2 * r2 / r_greater
        
        result, _ = integrate.quad(integrand, 0, 30, limit=100)
        return result
    
    def outer_integrand(r1):
        psi1 = orbital_1s(r1)
        return psi1 * psi1 * r1 * r1 * inner_integral(r1)
    
    # This is a 2D integral, will be slow
    print("Computing J numerically (this may take a moment)...")
    result, error = integrate.quad(outer_integrand, 0, 30, limit=100)
    print(f"Coulomb J (numerical): {result:.10f} Ha")
    return result


def compute_total_energy():
    """E_tot = 2*T + 2*V_nuc + J for He"""
    
    print("=" * 60)
    print("HELIUM HF ENERGY - NUMERICAL VERIFICATION")
    print("=" * 60)
    
    norm = test_normalization()
    print()
    
    T = test_kinetic_numerical()
    V = test_nuclear_numerical()
    print()
    
    J = test_coulomb_numerical()
    print()
    
    # He has 2 electrons in 1s
    E_1e = 2 * (T + V)
    E_2e = J  # Only one electron pair
    E_tot = E_1e + E_2e
    
    print("=" * 60)
    print("ENERGY BREAKDOWN:")
    print(f"  <T>_1s     = {T:.10f} Ha")
    print(f"  <V_nuc>_1s = {V:.10f} Ha")
    print(f"  J_1s,1s    = {J:.10f} Ha")
    print()
    print(f"  E_1e = 2*(T+V) = {E_1e:.10f} Ha")
    print(f"  E_2e = J       = {E_2e:.10f} Ha")
    print(f"  E_tot          = {E_tot:.10f} Ha")
    print()
    print(f"  Reference E    = -2.861679996 Ha")
    print(f"  Error          = {100 * (E_tot - (-2.861679996)) / 2.861679996:.4f}%")
    print("=" * 60)
    
    # Virial theorem: -<V>/<T> = 2
    V_tot = 2 * V + J
    print(f"\nVirial ratio: -<V_tot>/<T_tot> = {-V_tot/(2*T):.6f} (should be 2.0)")


if __name__ == "__main__":
    compute_total_energy()
