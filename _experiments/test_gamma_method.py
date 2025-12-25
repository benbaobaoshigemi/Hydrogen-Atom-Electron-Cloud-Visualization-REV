"""
Testing Analytical Cumulative Energy Calculation
=================================================
Uses incomplete gamma functions to compute E(R) - the cumulative orbital energy.

For STO density ρ(r) = Σ wₖ r^Nₖ exp(-αₖr), the nuclear attraction integral is:
E_nuc(R) = -4πZ ∫₀^R r·ρ(r) dr = -4πZ Σ wₖ γ(Nₖ+2, αₖR) / αₖ^(Nₖ+2)

Test case: Hydrogen 1s should give E(∞) = -0.5 Hartree exactly.
"""

import numpy as np
from scipy.special import gammainc, gamma as gamma_func
import json
import matplotlib.pyplot as plt

def lower_incomplete_gamma(a, x):
    """
    Lower incomplete gamma function: γ(a, x) = ∫₀^x t^(a-1) e^(-t) dt
    
    scipy.special.gammainc returns the REGULARIZED incomplete gamma: P(a,x) = γ(a,x)/Γ(a)
    So: γ(a, x) = gammainc(a, x) * Γ(a)
    """
    return gammainc(a, x) * gamma_func(a)

def expand_sto_density(basis):
    """
    Expand |R(r)|² = (Σᵢ cᵢNᵢr^(nᵢ-1)e^(-ζᵢr))² into sum of terms.
    
    Each term in the expansion has form: w * r^N * exp(-α*r)
    where:
    - w = 2 * cᵢ * cⱼ * Nᵢ * Nⱼ  (cross terms)
    - N = (nᵢ-1) + (nⱼ-1) = nᵢ + nⱼ - 2
    - α = ζᵢ + ζⱼ
    
    Returns: list of (w, N, alpha) tuples
    """
    from math import factorial
    
    def sto_norm(n, zeta):
        """Normalization factor for STO: sqrt((2ζ)^(2n+1) / (2n)!)"""
        return np.sqrt((2 * zeta) ** (2 * n + 1) / factorial(2 * n))
    
    terms = []
    n_basis = len(basis)
    
    for i in range(n_basis):
        ni = basis[i]['nStar']
        zi = basis[i]['zeta']
        ci = basis[i]['coeff']
        Ni = sto_norm(ni, zi)
        
        for j in range(n_basis):
            nj = basis[j]['nStar']
            zj = basis[j]['zeta']
            cj = basis[j]['coeff']
            Nj = sto_norm(nj, zj)
            
            # Product term: (ci*Ni*r^(ni-1)*e^(-zi*r)) * (cj*Nj*r^(nj-1)*e^(-zj*r))
            # = ci*cj*Ni*Nj * r^(ni+nj-2) * e^(-(zi+zj)*r)
            w = ci * cj * Ni * Nj
            N = (ni - 1) + (nj - 1)  # = ni + nj - 2
            alpha = zi + zj
            
            terms.append((w, N, alpha))
    
    return terms

def calc_nuclear_attraction(R, Z, density_terms):
    """
    Calculate cumulative nuclear attraction energy:
    E_nuc(R) = ∫₀^R ρ(r) * (-Z/r) * r² dr   [Note: NO 4π because STO already normalized]
             = -Z ∫₀^R r * ρ(r) dr
             = -Z Σₖ wₖ ∫₀^R r^(Nₖ+1) exp(-αₖr) dr
             = -Z Σₖ wₖ * γ(Nₖ+2, αₖR) / αₖ^(Nₖ+2)
    """
    E = 0.0
    for w, N, alpha in density_terms:
        a = N + 2
        x = alpha * R
        integral = lower_incomplete_gamma(a, x) / (alpha ** a)
        E += -Z * w * integral
    return E

def verify_normalization(density_terms):
    """
    Verify that ∫₀^∞ ρ(r) * r² dr = 1 (normalized to 1 electron).
    This equals Σₖ wₖ * Γ(Nₖ+3) / αₖ^(Nₖ+3)  [Note: NO 4π]
    """
    total = 0.0
    for w, N, alpha in density_terms:
        a = N + 3
        integral = gamma_func(a) / (alpha ** a)
        total += w * integral
    return total

# =============================================================================
# TEST 1: Hydrogen 1s
# =============================================================================
print("=" * 70)
print("TEST 1: Hydrogen 1s")
print("Expected: E(∞) = -0.5 Hartree")
print("=" * 70)

# H 1s: single STO with n*=1, ζ=1.0, c=1.0
h_basis = [{"nStar": 1, "zeta": 1.0, "coeff": 1.0}]
h_terms = expand_sto_density(h_basis)
print(f"Density expansion: {len(h_terms)} terms")
print(f"Terms: {h_terms}")

# Verify normalization
h_norm = verify_normalization(h_terms)
print(f"Normalization check: {h_norm:.10f} (should be 1.0)")

# Calculate E(R) curve
R_grid = np.linspace(0.01, 20, 500)
E_H = np.array([calc_nuclear_attraction(R, Z=1, density_terms=h_terms) for R in R_grid])

# For H, the asymptotic value should be -0.5 Hartree
E_H_inf = calc_nuclear_attraction(1000, Z=1, density_terms=h_terms)
print(f"E(R=1000) = {E_H_inf:.10f} Hartree")
print(f"Expected:  -0.5000000000 Hartree")
print(f"Error:     {abs(E_H_inf - (-0.5)):.2e}")

# =============================================================================
# TEST 2: Helium 1s (using Koga basis)
# =============================================================================
print("\n" + "=" * 70)
print("TEST 2: Helium 1s (Koga basis)")
print("Expected: ε_1s = -0.9179556 Hartree (from Koga)")
print("=" * 70)

# Load Koga data
with open('koga_basis_database.json', 'r', encoding='utf-8') as f:
    koga_db = json.load(f)

he_basis = koga_db['He']['orbitals']['1s']['basis']
he_epsilon = koga_db['He']['orbitals']['1s']['energy']
print(f"Koga He 1s basis: {len(he_basis)} terms")

he_terms = expand_sto_density(he_basis)
print(f"Density expansion: {len(he_terms)} terms")

he_norm = verify_normalization(he_terms)
print(f"Normalization check: {he_norm:.10f} (should be 1.0)")

# Calculate nuclear attraction only (ignoring e-e repulsion for now)
E_He_nuc = np.array([calc_nuclear_attraction(R, Z=2, density_terms=he_terms) for R in R_grid])
E_He_nuc_inf = calc_nuclear_attraction(1000, Z=2, density_terms=he_terms)
print(f"\nNuclear attraction only:")
print(f"E_nuc(R=1000) = {E_He_nuc_inf:.10f} Hartree")
print(f"This should be MORE NEGATIVE than ε because e-e repulsion is POSITIVE")
print(f"Koga ε_1s = {he_epsilon} Hartree")

# =============================================================================
# TEST 3: Plot comparison
# =============================================================================
print("\n" + "=" * 70)
print("Generating verification plot...")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# H plot
axes[0].plot(R_grid, E_H, 'b-', linewidth=2, label='E(R) by γ-function')
axes[0].axhline(y=-0.5, color='r', linestyle='--', label='Exact: -0.5 Hartree')
axes[0].set_xlim(0, 15)
axes[0].set_ylim(-0.6, 0)
axes[0].set_xlabel('R (a₀)')
axes[0].set_ylabel('Cumulative E (Hartree)')
axes[0].set_title('Hydrogen 1s: E(R) Test')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# He plot
axes[1].plot(R_grid, E_He_nuc, 'b-', linewidth=2, label='E_nuc(R) nuclear only')
axes[1].axhline(y=he_epsilon, color='r', linestyle='--', label=f'Koga ε = {he_epsilon:.4f}')
axes[1].axhline(y=E_He_nuc_inf, color='g', linestyle=':', label=f'E_nuc(∞) = {E_He_nuc_inf:.4f}')
axes[1].set_xlim(0, 10)
axes[1].set_ylim(-5, 0)
axes[1].set_xlabel('R (a₀)')
axes[1].set_ylabel('Cumulative E (Hartree)')
axes[1].set_title('Helium 1s: E_nuc(R) [missing e-e repulsion]')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('test_gamma_method.png', dpi=150)
print("Plot saved to test_gamma_method.png")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("If H test passes (E→-0.5), the γ-function method is VALIDATED.")
print("He shows E_nuc < ε because we haven't added e-e repulsion yet.")
