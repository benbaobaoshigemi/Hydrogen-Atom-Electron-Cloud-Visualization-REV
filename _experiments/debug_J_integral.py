# -*- coding: utf-8 -*-
"""
DEBUG: Verify Y_0 integral for Helium 1s
========================================

For a 1s orbital with Z_eff=1.6875 (Slater):
R(r) = 2 * Z^(3/2) * exp(-Z*r)

J_11 = ∫∫ R²(r1) (1/r_>) R²(r2) r1² r2² dr1 dr2
     = 5/8 * Z (for hydrogenic)

With Z=1.6875: J = 5/8 * 1.6875 = 1.055 Ha
"""

import numpy as np
from scipy.integrate import quad, dblquad

# ==============================================================================
# HELIUM 1s FROM KOGA BASIS
# ==============================================================================

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
    """STO normalization"""
    return (2 * zeta) ** n * np.sqrt(2 * zeta / factorial(2 * n))


def R_1s(r):
    """Koga 1s radial wavefunction"""
    if isinstance(r, np.ndarray):
        result = np.zeros_like(r)
    else:
        result = 0.0
    
    for term in HE_1S_BASIS:
        n = term['n']
        zeta = term['zeta']
        c = term['coeff']
        N = sto_norm(n, zeta)
        result += c * N * (r ** (n - 1)) * np.exp(-zeta * r)
    
    return result


def test_normalization():
    """Test ∫ R² r² dr = 1"""
    def integrand(r):
        return R_1s(r) ** 2 * r ** 2
    
    result, error = quad(integrand, 0, 100, limit=200)
    print(f"Normalization: {result:.10f} (should be 1.0)")
    return result


def compute_Y0_analytical_simple(r_val, R_func, r_max=100):
    """
    Compute Y_0(1s,1s; r) using direct numerical integration.
    
    Y_0(r) = ∫_0^∞ min(r, r') / max(r, r') * R(r')² r'² dr'
           = (1/r) ∫_0^r R(r')² r'² dr' + r ∫_r^∞ R(r')² r' dr'
    
    Simpler form:
    Y_0(r) = (1/r) * I_inner + r * I_outer
    
    where:
    I_inner = ∫_0^r ρ(r') r'² dr'
    I_outer = ∫_r^∞ ρ(r') r' dr'  (one less power of r')
    """
    def rho(r):
        return R_func(r) ** 2
    
    # Inner integral: ∫_0^r ρ(r') r'² dr'  [cumulative charge]
    def inner_integrand(rp):
        return rho(rp) * rp ** 2
    
    # Outer integral: ∫_r^∞ ρ(r') r' dr'
    def outer_integrand(rp):
        return rho(rp) * rp
    
    if r_val < 1e-10:
        # At r=0, Y_0 = 0
        return 0.0
    
    I_inner, _ = quad(inner_integrand, 0, r_val, limit=200)
    I_outer, _ = quad(outer_integrand, r_val, r_max, limit=200)
    
    Y0 = I_inner / r_val + r_val * I_outer
    
    return Y0


def compute_J_integral():
    """
    Compute J = ∫ R(r)² * [Y_0(r)/r] * r² dr
             = ∫ ρ(r) * [Y_0(r)/r] * r² dr
             = ∫ ρ(r) * J(r) * r² dr
    
    where J(r) = Y_0(r) / r is the electrostatic potential.
    """
    def integrand(r):
        if r < 1e-10:
            return 0.0
        rho = R_1s(r) ** 2
        Y0 = compute_Y0_analytical_simple(r, R_1s)
        J_potential = Y0 / r
        return rho * J_potential * r ** 2
    
    print("Computing J integral (this may take a minute)...")
    result, error = quad(integrand, 0, 50, limit=300)
    print(f"J_11 = {result:.10f} Ha")
    print(f"Integration error estimate: {error:.2e}")
    return result


def compute_J_double_integral():
    """
    Alternative: Direct 2D integral.
    
    J = ∫∫ R(r1)² R(r2)² (1/r_>) r1² r2² dr1 dr2
    """
    def integrand(r2, r1):
        # Note: dblquad integrates y first, then x
        rho1 = R_1s(r1) ** 2
        rho2 = R_1s(r2) ** 2
        r_greater = max(r1, r2)
        return rho1 * rho2 * r1 ** 2 * r2 ** 2 / r_greater
    
    print("Computing J by 2D integral (slow!)...")
    result, error = dblquad(integrand, 0, 30, 0, 30, epsabs=1e-6, epsrel=1e-6)
    print(f"J_11 (2D) = {result:.10f} Ha")
    return result


def main():
    print("=" * 60)
    print("HELIUM J INTEGRAL VERIFICATION")
    print("=" * 60)
    
    test_normalization()
    print()
    
    J = compute_J_integral()
    print()
    
    # Expected from orbital energy:
    # ε = h + J (since K=J for same orbital)
    # E_tot = 2h + J = 2(ε - J) + J = 2ε - J
    # → J = 2ε - E_tot
    eps_ref = -0.9179556
    E_tot_ref = -2.861679996
    J_expected = 2 * eps_ref - E_tot_ref
    print(f"Expected J from ε and E_tot: {J_expected:.8f} Ha")
    
    # For comparison: hydrogenic with Z=1.6875
    Z_eff = 1.6875
    J_hydro = 5/8 * Z_eff
    print(f"Hydrogenic J (5/8 * Z_eff): {J_hydro:.8f} Ha")


if __name__ == "__main__":
    main()
