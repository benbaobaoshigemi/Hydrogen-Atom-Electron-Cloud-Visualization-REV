# -*- coding: utf-8 -*-
"""
CORRECTED Y_0 INTEGRAL FOR HELIUM
=================================

The correct formula for Y_0:

Y_0(r) = ∫_0^∞ (r_< / r_>) * ρ(r') * r'^2 dr'

where r_< = min(r, r'), r_> = max(r, r')

Split:
- For r' < r: r_< = r', r_> = r, so factor = r'/r
- For r' > r: r_< = r, r_> = r', so factor = r/r'

Y_0(r) = (1/r) ∫_0^r r' * ρ(r') * r'^2 dr' + r * ∫_r^∞ (1/r') * ρ(r') * r'^2 dr'
       = (1/r) ∫_0^r ρ(r') r'^3 dr' + r ∫_r^∞ ρ(r') r' dr'

For the Coulomb potential:
V(r) = Y_0(r) / r

And J integral:
J = ∫_0^∞ ρ(r) * V(r) * r^2 dr = ∫_0^∞ ρ(r) * Y_0(r) / r * r^2 dr
  = ∫_0^∞ ρ(r) * Y_0(r) * r dr
"""

import numpy as np
from scipy.integrate import quad

# KOGA He 1s basis
HE_1S_BASIS = [
    {'n': 2, 'zeta': 6.438865513242302, 'coeff': 0.0008103},
    {'n': 1, 'zeta': 3.385077039750975, 'coeff': 0.0798826},
    {'n': 1, 'zeta': 2.178370004614139, 'coeff': 0.180161},
    {'n': 1, 'zeta': 1.4553870053179185, 'coeff': 0.7407925},
    {'n': 2, 'zeta': 1.3552466748849417, 'coeff': 0.0272015},
]


def factorial(n):
    if n <= 1:
        return 1.0
    return n * factorial(n - 1)


def sto_norm(n, zeta):
    return (2 * zeta) ** n * np.sqrt(2 * zeta / factorial(2 * n))


def R_1s(r):
    result = 0.0
    for term in HE_1S_BASIS:
        n = term['n']
        zeta = term['zeta']
        c = term['coeff']
        N = sto_norm(n, zeta)
        result += c * N * (r ** (n - 1)) * np.exp(-zeta * r)
    return result


def rho(r):
    """Radial density ρ(r) = R(r)²"""
    return R_1s(r) ** 2


def Y0_correct(r_val, r_max=100):
    """
    CORRECT Y_0 formula:
    Y_0(r) = (1/r) ∫_0^r ρ(r') r'^3 dr' + r ∫_r^∞ ρ(r') r' dr'
    """
    if r_val < 1e-12:
        return 0.0
    
    # Inner: (1/r) ∫_0^r ρ(r') r'^3 dr'
    def inner_integrand(rp):
        return rho(rp) * rp ** 3
    
    I_inner, _ = quad(inner_integrand, 0, r_val, limit=200)
    
    # Outer: r ∫_r^∞ ρ(r') r' dr'
    def outer_integrand(rp):
        return rho(rp) * rp
    
    I_outer, _ = quad(outer_integrand, r_val, r_max, limit=200)
    
    return I_inner / r_val + r_val * I_outer


def compute_J_correct():
    """
    J = ∫_0^∞ ρ(r) * Y_0(r) * r dr
    
    (Note: not r², because the volume element r² is in ρ already? 
     No wait, ρ(r) = R²(r), so we need r² from volume element.
     
     Actually: J = ∫∫ ρ(r1) ρ(r2) / r_> r1² r2² dr1 dr2
     We're computing: ∫ ρ(r1) r1² [∫ ρ(r2) / r_> r2² dr2] dr1
                    = ∫ ρ(r1) r1² [Y_0(r1) / r1] dr1
                    = ∫ ρ(r1) r1 Y_0(r1) dr1
    
    Hmm, that's not right either. Let me re-derive.
    
    J = ∫∫ ρ(r1) ρ(r2) (1/r_>) r1² r2² dr1 dr2
    
    Define: Y_0(r) = ∫_0^∞ (r_< / r_>) ρ(r') r'^2 dr' 
                   (this is the standard Hartree Y function)
    
    Then the potential at r due to ρ is: V(r) = ∫ ρ(r') r'^2 / r_> dr'
    
    But that's not exactly Y_0... Let me check the textbook definition.
    
    Standard Hartree Y_k function:
    Y_k(r) = r ∫_0^∞ (r_<^k / r_>^(k+1)) ρ(r') r'^2 dr'
    
    For k=0:
    Y_0(r) = r ∫_0^∞ (1/r_>) ρ(r') r'^2 dr'
           = r [ (1/r) ∫_0^r ρ r'^2 dr' + ∫_r^∞ (1/r') ρ r'^2 dr' ]
           = ∫_0^r ρ r'^2 dr' + r ∫_r^∞ ρ r' dr'
    
    Then: V(r) = Y_0(r) / r
    
    And: J = ∫_0^∞ ρ(r) V(r) r² dr = ∫_0^∞ ρ(r) [Y_0(r)/r] r² dr
           = ∫_0^∞ ρ(r) Y_0(r) r dr
    )    """
    def integrand(r):
        if r < 1e-12:
            return 0.0
        return rho(r) * Y0_correct(r) * r
    
    print("Computing J with CORRECT Y_0...")
    J, error = quad(integrand, 0, 50, limit=300)
    print(f"J (correct) = {J:.10f} Ha")
    print(f"Error estimate: {error:.2e}")
    return J


def Y0_standard(r_val, r_max=100):
    """
    Standard Y_0 definition from Slater's book:
    Y_0(r) = ∫_0^r ρ(r') r'^2 dr' + r ∫_r^∞ ρ(r') r' dr'
    """
    if r_val < 1e-12:
        # At r=0, first term is 0, second term is r * ∫_0^∞ ρ r' dr' → 0
        return 0.0
    
    # First term: ∫_0^r ρ(r') r'^2 dr'  (cumulative enclosed charge * r)
    def term1_integrand(rp):
        return rho(rp) * rp ** 2
    
    I1, _ = quad(term1_integrand, 0, r_val, limit=200)
    
    # Second term: r ∫_r^∞ ρ(r') r' dr'
    def term2_integrand(rp):
        return rho(rp) * rp
    
    I2, _ = quad(term2_integrand, r_val, r_max, limit=200)
    
    return I1 + r_val * I2


def compute_J_standard():
    """
    J = ∫_0^∞ ρ(r) [Y_0(r)/r] r² dr
      = ∫_0^∞ ρ(r) Y_0(r) r dr
    """
    def integrand(r):
        if r < 1e-12:
            return 0.0
        return rho(r) * Y0_standard(r) * r
    
    print("\nComputing J with STANDARD Y_0...")
    J, error = quad(integrand, 0, 50, limit=300)
    print(f"J (standard) = {J:.10f} Ha")
    print(f"Error estimate: {error:.2e}")
    return J


def main():
    print("=" * 60)
    print("Y_0 INTEGRAL DEBUGGING")
    print("=" * 60)
    
    # Check normalization
    norm, _ = quad(lambda r: rho(r) * r**2, 0, 100, limit=200)
    print(f"Normalization: {norm:.10f}")
    
    # Test Y_0 at a few points
    print("\nY_0 values at sample points:")
    for r in [0.1, 0.5, 1.0, 2.0, 5.0]:
        y_correct = Y0_correct(r)
        y_standard = Y0_standard(r)
        print(f"  r={r:.1f}: Y0_correct={y_correct:.6f}, Y0_standard={y_standard:.6f}")
    
    print()
    J_correct = compute_J_correct()
    J_standard = compute_J_standard()
    
    # Expected
    eps_ref = -0.9179556
    E_tot_ref = -2.861679996
    J_expected = 2 * eps_ref - E_tot_ref
    print(f"\nExpected J: {J_expected:.8f} Ha")
    
    # Also compute T and V for verification
    print("\n=== ONE-ELECTRON INTEGRALS ===")
    
    # Kinetic energy <T> for 1s
    # T = <R|-1/2 ∇²|R> = ∫ R * (-1/2) * [R'' + 2/r R' - 0] * r² dr
    # For numerical: use known result that T = <ψ|T|ψ> = E_tot - <V>
    
    # Nuclear <V> = -Z ∫ R² / r * r² dr = -Z ∫ R² r dr
    Z = 2
    V_nuc, _ = quad(lambda r: -Z * rho(r) * r, 0, 100, limit=200)
    print(f"<V_nuc> = {V_nuc:.8f} Ha")
    
    # From virial: T = -E_tot (for exact HF)
    # or T = (E_tot - V) where V is total potential
    
    # Let's compute h = T + V directly
    # h = ε - J (for He 1s)
    h_expected = eps_ref - J_expected
    print(f"h = ε - J = {h_expected:.8f} Ha")
    print(f"This implies T = h - V_nuc = {h_expected - V_nuc:.8f} Ha")
    
    # Verify E_tot = 2h + J
    E_check = 2 * h_expected + J_expected
    print(f"\nE_tot check: 2h + J = {E_check:.8f} Ha")
    print(f"E_tot (ref) = {E_tot_ref:.8f} Ha")


if __name__ == "__main__":
    main()
