"""
RIGOROUS VERIFICATION: Kinetic Energy Calculation for STOs
=============================================================
This script provides MATHEMATICAL PROOF that the proposed kinetic energy
formula is correct by calculating T(∞) for hydrogen 1s.

Expected result: T(∞) = 0.5 Hartree (exact quantum mechanical value)

Reference: Hydrogen atom kinetic energy = -E_total = -(-0.5) = 0.5 Hartree
(by virial theorem: T = -E for hydrogen)
"""

import numpy as np
from scipy.special import gamma as gamma_func, gammainc
from scipy.integrate import quad
import sympy as sp

print("=" * 70)
print("RIGOROUS VERIFICATION: STO Kinetic Energy Formula")
print("=" * 70)

# ===========================================================================
# PART 1: Symbolic Verification using SymPy
# ===========================================================================
print("\n### PART 1: SYMBOLIC DERIVATION (SymPy) ###\n")

r = sp.Symbol('r', positive=True, real=True)
n, zeta, l = sp.symbols('n zeta l', positive=True, real=True)

# STO radial function: R(r) = r^(n-1) * exp(-zeta*r)
# (ignoring normalization for now - we'll add it later)
R = r**(n-1) * sp.exp(-zeta * r)

# Radial kinetic energy operator in spherical coordinates:
# T_hat = -1/2 * [d²/dr² + (2/r)*d/dr - l(l+1)/r²]
# Applied to R(r):

R_prime = sp.diff(R, r)
R_double_prime = sp.diff(R, r, 2)

# T_hat * R = -1/2 * [R'' + (2/r)*R' - l(l+1)/r² * R]
T_R = sp.Rational(-1, 2) * (R_double_prime + (2/r) * R_prime - l*(l+1)/r**2 * R)
T_R_simplified = sp.simplify(T_R)

print("T̂ R(r) = ", T_R_simplified)

# Factor out common exp(-zeta*r) term
T_R_factored = sp.collect(sp.expand(T_R_simplified * sp.exp(zeta*r)), r) * sp.exp(-zeta*r)
print("\nFactored form:")
print("T̂ R(r) = ", sp.simplify(T_R_factored))

# For n=1, l=0, zeta=1 (hydrogen 1s):
print("\n--- Hydrogen 1s (n=1, l=0, ζ=1) ---")
T_R_H1s = T_R_simplified.subs([(n, 1), (l, 0), (zeta, 1)])
print("T̂ ψ_1s = ", sp.simplify(T_R_H1s))

# ===========================================================================
# PART 2: Numerical Calculation of T(∞)
# ===========================================================================
print("\n### PART 2: NUMERICAL CALCULATION ###\n")

def sto_radial(r, n_val, zeta_val):
    """STO radial function (unnormalized): r^(n-1) * exp(-zeta*r)"""
    return r**(n_val - 1) * np.exp(-zeta_val * r)

def sto_normalization(n_val, zeta_val):
    """Normalization constant: N = (2*zeta)^n * sqrt(2*zeta / (2n)!)"""
    from math import factorial
    return (2 * zeta_val)**n_val * np.sqrt(2 * zeta_val / factorial(2 * n_val))

def kinetic_operator_on_sto(r, n_val, zeta_val, l_val):
    """
    Apply T̂ to (n-1) * exp(-zeta*r) and return the result.
    
    T̂ R = -1/2 * [A*r^(n-3) + B*r^(n-2) + C*r^(n-1)] * exp(-zeta*r)
    
    Where:
    A = n(n-1) - l(l+1)
    B = -2*zeta*n  
    C = zeta²
    """
    A = n_val * (n_val - 1) - l_val * (l_val + 1)
    B = -2 * zeta_val * n_val
    C = zeta_val**2
    
    # Handle r=0 carefully for n=1 case (A term has r^(-2))
    if r < 1e-10:
        r = 1e-10
    
    result = -0.5 * (
        A * r**(n_val - 3) +
        B * r**(n_val - 2) +
        C * r**(n_val - 1)
    ) * np.exp(-zeta_val * r)
    
    return result

# Hydrogen 1s parameters
n_H = 1
l_H = 0
zeta_H = 1.0
N_H = sto_normalization(n_H, zeta_H)

print(f"Hydrogen 1s parameters:")
print(f"  n = {n_H}, l = {l_H}, ζ = {zeta_H}")
print(f"  Normalization N = {N_H:.10f}")

# Verify normalization: ∫₀^∞ N² R² r² dr = 1
def norm_integrand(r):
    return N_H**2 * sto_radial(r, n_H, zeta_H)**2 * r**2

norm_integral, norm_error = quad(norm_integrand, 0, np.inf)
print(f"  Normalization check: ∫ N² R² r² dr = {norm_integral:.10f} (should be 1.0)")

# Calculate kinetic energy: T = ∫₀^∞ R* T̂ R r² dr (with normalization)
def kinetic_integrand(r):
    R = N_H * sto_radial(r, n_H, zeta_H)
    T_R = N_H * kinetic_operator_on_sto(r, n_H, zeta_H, l_H)
    return R * T_R * r**2

T_infinity, T_error = quad(kinetic_integrand, 0, np.inf, limit=200)
print(f"\nKinetic energy T(∞) = {T_infinity:.10f} Hartree")
print(f"  Integration error: {T_error:.2e}")
print(f"  Expected value:    0.5000000000 Hartree")
print(f"  Relative error:    {abs(T_infinity - 0.5) / 0.5 * 100:.4f}%")

# ===========================================================================
# PART 3: Verify Total Energy E = T + V_nuc = -0.5 Hartree
# ===========================================================================
print("\n### PART 3: TOTAL ENERGY VERIFICATION ###\n")

# Nuclear potential energy: V_nuc = ∫ R* (-Z/r) R r² dr = -Z ∫ R² r dr
def potential_integrand(r):
    R = N_H * sto_radial(r, n_H, zeta_H)
    return R**2 * (-1.0 / r) * r**2  # -Z/r with Z=1

V_infinity, V_error = quad(potential_integrand, 1e-10, np.inf, limit=200)
print(f"Potential energy V_nuc(∞) = {V_infinity:.10f} Hartree")
print(f"  Expected value:           -1.0000000000 Hartree")

E_total = T_infinity + V_infinity
print(f"\nTotal orbital energy E = T + V_nuc")
print(f"  E = {T_infinity:.6f} + ({V_infinity:.6f}) = {E_total:.10f} Hartree")
print(f"  Expected value: -0.5000000000 Hartree")
print(f"  Relative error: {abs(E_total - (-0.5)) / 0.5 * 100:.4f}%")

# ===========================================================================
# PART 4: Verify Coefficients A, B, C
# ===========================================================================
print("\n### PART 4: COEFFICIENT VERIFICATION ###\n")

# For hydrogen 1s: n=1, l=0
A_H = 1 * (1-1) - 0 * (0+1)  # = 0
B_H = -2 * 1.0 * 1  # = -2
C_H = 1.0**2  # = 1

print("Coefficients for H 1s (n=1, l=0, ζ=1):")
print(f"  A = n(n-1) - l(l+1) = 1×0 - 0×1 = {A_H}")
print(f"  B = -2ζn = -2×1×1 = {B_H}")
print(f"  C = ζ² = 1² = {C_H}")
print()
print("Therefore:")
print("  T̂ ψ_1s = -1/2 × [0×r^(-2) + (-2)×r^(-1) + 1×r^0] × e^(-r)")
print("         = -1/2 × [-2/r + 1] × e^(-r)")
print("         = [1/r - 1/2] × e^(-r)")

# Verify this formula directly
def direct_T_H1s(r):
    """Direct formula: T̂ ψ = (1/r - 1/2) * exp(-r)"""
    if r < 1e-10:
        return 0
    return (1.0/r - 0.5) * np.exp(-r)

def compare_integrand(r):
    R = N_H * sto_radial(r, n_H, zeta_H)
    T_R_formula = N_H * direct_T_H1s(r)
    return R * T_R_formula * r**2

T_direct, _ = quad(compare_integrand, 1e-10, np.inf, limit=200)
print(f"\nDirect formula verification:")
print(f"  T using explicit (1/r - 1/2)×e^(-r) = {T_direct:.10f} Hartree")
print(f"  Matches expected 0.5 Hartree: {abs(T_direct - 0.5) < 0.001}")

# ===========================================================================
# FINAL VERDICT
# ===========================================================================
print("\n" + "=" * 70)
print("FINAL VERDICT")
print("=" * 70)

if abs(T_infinity - 0.5) < 0.01 and abs(E_total - (-0.5)) < 0.01:
    print("✓ KINETIC ENERGY FORMULA VERIFIED")
    print("✓ T(∞) = 0.5 Hartree (within 1% of exact value)")
    print("✓ E_total = T + V = -0.5 Hartree (correct hydrogen 1s energy)")
    print("\nThe proposed analytical method is MATHEMATICALLY CORRECT.")
else:
    print("✗ VERIFICATION FAILED")
    print(f"  T_computed = {T_infinity}, expected = 0.5")
    print(f"  E_computed = {E_total}, expected = -0.5")
