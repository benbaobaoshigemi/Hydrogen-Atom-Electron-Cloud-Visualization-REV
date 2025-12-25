import numpy as np
import matplotlib.pyplot as plt

# Koga Helium 1s Data
he_basis = [
    {"nStar": 2, "zeta": 6.437494, "coeff": 0.0008103},
    {"nStar": 1, "zeta": 3.384356, "coeff": 0.0798826},
    {"nStar": 1, "zeta": 2.177906, "coeff": 0.180161},
    {"nStar": 1, "zeta": 1.455077, "coeff": 0.7407925},
    {"nStar": 2, "zeta": 1.354958, "coeff": 0.0272015}
]

def factorial(n):
    if n <= 1: return 1
    return n * factorial(n - 1)

def slater_radial(r, basis):
    val = 0.0
    for term in basis:
        n = term['nStar']
        z = term['zeta']
        c = term['coeff']
        # Normalization for R(r) (not wavefunction, just radial part)
        # R_nl(r) = N * r^(n-1) * exp(-z*r)
        # Normalization N: \int R^2 r^2 dr = 1
        # \int r^(2n-2) * exp(-2zr) * r^2 dr = \int r^(2n) * exp(-2zr) dr
        # = (2n)! / (2z)^(2n+1)
        # So N = sqrt( (2z)^(2n+1) / (2n)! )
        
        N = np.sqrt( (2*z)**(2*n + 1) / factorial(2*n) )
        val += c * N * (r**(n-1)) * np.exp(-z*r)
    return val

def compute_density(r_grid, basis):
    # Electron density of ONE electron in this orbital
    # rho(r) = |R(r)|^2  (spherically averaged, integrating out angles gives 4pi factor? 
    # Usually radial probability density P(r) = r^2 R^2.
    # Charge density rho(r) such that \int rho(r) 4pi r^2 dr = 1.
    # So rho(r) = R(r)^2 / (4pi) ? 
    # Wait, Potential V(r) calculation usually sums shells.
    # Let's use Q(r) = enclosed charge.
    
    R_vals = np.array([slater_radial(r, basis) for r in r_grid])
    # P(r) dr is probability in shell. P(r) = R^2 * r^2.
    # Enclosed charge Q(r) = \int_0^r R(x)^2 x^2 dx.
    return R_vals

def compute_hartree_potential(r_grid, basis):
    # V_H(r) = integral ( rho(r') / |r-r'| ) d3r'
    # For spherical rho(r), this is equivalent to classical shell theorem.
    # V(r) = (1/r) * integral_0^r rho(x) 4pi x^2 dx  +  integral_r^inf (rho(x) 4pi x^2 / x) dx
    
    dr = r_grid[1] - r_grid[0]
    R_vals = np.array([slater_radial(r, basis) for r in r_grid])
    P_vals = (R_vals**2) * (r_grid**2) # Radial probability density
    
    # Total charge check
    total_prob = np.sum(P_vals * dr)
    print(f"Total Probability (Normalization Check): {total_prob:.6f}")
    
    V_hartree = np.zeros_like(r_grid)
    
    for i, r in enumerate(r_grid):
        # Inner part: charge inside radius r behaves like point charge at origin
        # Q_in = \int_0^r P(x) dx
        Q_in = np.sum(P_vals[:i] * dr)
        term1 = Q_in / r if r > 1e-6 else 0
        
        # Outer part: charge outside radius r does not contribute to force, but contributes to POTENTIAL constant
        # V_out = \int_r^inf (P(x)/x) dx
        term2 = np.sum((P_vals[i:] / r_grid[i:]) * dr)
        
        V_hartree[i] = term1 + term2
        
    return V_hartree

# Simulation
r = np.linspace(0.01, 10, 500)
V_h = compute_hartree_potential(r, he_basis)
Z = 2
V_nuc = -Z / r

# Effective Potential for ONE electron
# The other electron creates the shielding.
# In Helium ground state, both electrons are in 1s.
# So the shielding potential is V_hartree calculated from the 1s orbital density.
# V_eff = V_nuc + V_hartree
V_eff = V_nuc + V_h 

# Plot
plt.figure(figsize=(10, 6))
plt.plot(r, V_nuc, '--', label='Bare Nucleus (-2/r)')
plt.plot(r, V_h, label='Hartree Shielding (from 1s)')
plt.plot(r, V_eff, 'r-', linewidth=2, label='Effective Potential')
plt.ylim(-5, 5)
plt.xlabel('r (a0)')
plt.ylabel('Potential (Hartree)')
plt.title('Helium Effective Potential Verification')
plt.legend()
plt.grid(True)
plt.savefig('verify_hartree.png')
print("Verification plot saved to verify_hartree.png")
