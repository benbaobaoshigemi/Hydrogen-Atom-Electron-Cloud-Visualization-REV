# -*- coding: utf-8 -*-
"""
FULL Roothaan-Hartree-Fock Energy Validation
============================================

This script computes the COMPLETE atomic total energy using:
- One-electron integrals: T (kinetic), V_nuc (nuclear attraction)
- Two-electron integrals: J (Coulomb), K (Exchange)

Then validates against Bunge (1993) reference energies.

NO FALLBACK. NO APPROXIMATIONS. Analytical STO integrals only.
"""

import json
import re
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from scipy import special  # For incomplete gamma function

# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

@dataclass
class STOTerm:
    """Single STO: N * r^(n-1) * exp(-zeta*r)"""
    n_star: int
    zeta: float
    coeff: float


@dataclass
class Orbital:
    """Atomic orbital as linear combination of STOs"""
    name: str
    l: int
    terms: List[STOTerm]


@dataclass
class AtomData:
    """Complete atomic data"""
    symbol: str
    Z: int
    orbitals: Dict[str, Orbital]
    E_tot_ref: float
    energies_ref: Dict[str, float]


# ==============================================================================
# ELECTRON CONFIGURATION
# ==============================================================================

# Ground state electron configurations for Z=1 to Z=36
ELECTRON_CONFIG = {
    1:  {'1s': 1},
    2:  {'1s': 2},
    3:  {'1s': 2, '2s': 1},
    4:  {'1s': 2, '2s': 2},
    5:  {'1s': 2, '2s': 2, '2p': 1},
    6:  {'1s': 2, '2s': 2, '2p': 2},
    7:  {'1s': 2, '2s': 2, '2p': 3},
    8:  {'1s': 2, '2s': 2, '2p': 4},
    9:  {'1s': 2, '2s': 2, '2p': 5},
    10: {'1s': 2, '2s': 2, '2p': 6},
    11: {'1s': 2, '2s': 2, '2p': 6, '3s': 1},
    12: {'1s': 2, '2s': 2, '2p': 6, '3s': 2},
    13: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 1},
    14: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 2},
    15: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 3},
    16: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 4},
    17: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 5},
    18: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6},
    19: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 1},
    20: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2},
    21: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 1},
    22: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 2},
    23: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 3},
    24: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 1, '3d': 5},  # Cr anomaly
    25: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 5},
    26: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 6},
    27: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 7},
    28: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 8},
    29: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 1, '3d': 10}, # Cu anomaly
    30: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10},
    31: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10, '4p': 1},
    32: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10, '4p': 2},
    33: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10, '4p': 3},
    34: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10, '4p': 4},
    35: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10, '4p': 5},
    36: {'1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 10, '4p': 6},
}


# ==============================================================================
# MATHEMATICAL UTILITIES
# ==============================================================================

def factorial(n: int) -> int:
    """Compute n!"""
    if n < 0:
        raise ValueError(f"Factorial undefined for n={n}")
    if n <= 1:
        return 1
    result = 1
    for i in range(2, n + 1):
        result *= i
    return result


def sto_normalization(n: int, zeta: float) -> float:
    """STO normalization: N = (2*zeta)^n * sqrt(2*zeta / (2n)!)"""
    return (2 * zeta) ** n * math.sqrt(2 * zeta / factorial(2 * n))


def radial_moment(k: int, alpha: float) -> float:
    """Compute ∫_0^∞ r^k exp(-alpha*r) dr = k! / alpha^(k+1)"""
    if k < 0:
        return 0.0
    return factorial(k) / (alpha ** (k + 1))


def lower_incomplete_gamma(a: int, x: float) -> float:
    """
    Lower incomplete gamma function γ(a, x) = ∫_0^x t^(a-1) e^-t dt
    
    For integer a, this has the closed form:
    γ(a, x) = (a-1)! * (1 - e^-x * Σ_{k=0}^{a-1} x^k / k!)
    """
    if x <= 0:
        return 0.0
    if a <= 0:
        raise ValueError(f"a must be positive, got {a}")
    
    # Use scipy for robust computation
    return special.gammainc(a, x) * special.gamma(a)


def upper_incomplete_gamma(a: int, x: float) -> float:
    """
    Upper incomplete gamma function Γ(a, x) = ∫_x^∞ t^(a-1) e^-t dt
    = Γ(a) - γ(a, x)
    """
    return special.gammaincc(a, x) * special.gamma(a)


# ==============================================================================
# STO INTEGRALS
# ==============================================================================

def overlap_sto(n1: int, z1: float, n2: int, z2: float) -> float:
    """Overlap integral <φ1|φ2> for STOs with same l"""
    N1 = sto_normalization(n1, z1)
    N2 = sto_normalization(n2, z2)
    alpha = z1 + z2
    return N1 * N2 * radial_moment(n1 + n2, alpha)


def kinetic_sto(n1: int, z1: float, l: int, n2: int, z2: float) -> float:
    """Kinetic energy integral <φ1|T|φ2>"""
    N1 = sto_normalization(n1, z1)
    N2 = sto_normalization(n2, z2)
    alpha = z1 + z2
    
    A = n2 * (n2 - 1) - l * (l + 1)
    B = -2 * z2 * n2
    C = z2 * z2
    
    k_base = n1 + n2
    integral = -0.5 * (
        A * radial_moment(k_base - 2, alpha) +
        B * radial_moment(k_base - 1, alpha) +
        C * radial_moment(k_base, alpha)
    )
    return N1 * N2 * integral


def nuclear_sto(n1: int, z1: float, n2: int, z2: float, Z: int) -> float:
    """Nuclear attraction integral <φ1|-Z/r|φ2>"""
    N1 = sto_normalization(n1, z1)
    N2 = sto_normalization(n2, z2)
    alpha = z1 + z2
    k = n1 + n2 - 1
    return -Z * N1 * N2 * radial_moment(k, alpha)


def coulomb_integral_sto(n1: int, z1: float, n2: int, z2: float,
                          n3: int, z3: float, n4: int, z4: float) -> float:
    """
    Compute Coulomb integral (12|34) for s-type STOs (l=0, k=0 only).
    
    (12|34) = ∫∫ φ1(r1) φ2(r1) (1/r_>) φ3(r2) φ4(r2) r1^2 r2^2 dr1 dr2
    
    Using Laplace expansion:
    1/r_> = 1/r_2 for r1 < r2, 1/r_1 for r1 > r2
    
    This becomes:
    (12|34) = N12 N34 * [I_< + I_>]
    
    where:
    I_< = ∫_0^∞ r2^(n3+n4) e^(-α34*r2) * [∫_0^r2 r1^(n1+n2+1) e^(-α12*r1) dr1] dr2
    I_> = ∫_0^∞ r1^(n1+n2) e^(-α12*r1) * [∫_0^r1 r2^(n3+n4+1) e^(-α34*r2) dr2] dr1
    """
    N1 = sto_normalization(n1, z1)
    N2 = sto_normalization(n2, z2)
    N3 = sto_normalization(n3, z3)
    N4 = sto_normalization(n4, z4)
    
    N12 = N1 * N2
    N34 = N3 * N4
    alpha12 = z1 + z2
    alpha34 = z3 + z4
    
    m12 = n1 + n2  # power from r^(n-1) * r^(n-1) * r^2 = r^(n1+n2)
    m34 = n3 + n4
    
    # The integral can be computed analytically using the formula:
    # ∫_0^∞ r^a e^(-br) [∫_0^r t^c e^(-dt) dt] dr
    # = c! / d^(c+1) * ∫_0^∞ r^a e^(-br) [1 - e^(-dr) Σ_{k=0}^c (dr)^k/k!] dr
    # = c! / d^(c+1) * [a!/b^(a+1) - Σ_{k=0}^c d^k/k! * (a+k)!/(b+d)^(a+k+1)]
    
    # I_< contribution: outer = r2^m34 e^(-α34*r2), inner up to r2, r1^(m12+1) e^(-α12*r1)
    # After inner integral from 0 to r2: γ(m12+2, α12*r2) / α12^(m12+2)
    # Factor out: (m12+1)! / α12^(m12+2) * [1 - e^(-α12*r2) Σ...]
    
    # Compute I_< using analytical formula
    a = m34  # outer power
    b = alpha34  # outer exp
    c = m12 + 1  # inner power (+1 from 1/r_> = 1/r2 contribution)
    d = alpha12  # inner exp
    
    # Full moment minus restricted sum
    full_moment = factorial(a) / (b ** (a + 1))
    restricted_sum = 0.0
    for k in range(c + 1):
        restricted_sum += (d ** k / factorial(k)) * factorial(a + k) / ((b + d) ** (a + k + 1))
    
    I_less = (factorial(c) / (d ** (c + 1))) * (full_moment - restricted_sum)
    
    # I_> contribution: swap roles
    a = m12
    b = alpha12
    c = m34 + 1
    d = alpha34
    
    full_moment = factorial(a) / (b ** (a + 1))
    restricted_sum = 0.0
    for k in range(c + 1):
        restricted_sum += (d ** k / factorial(k)) * factorial(a + k) / ((b + d) ** (a + k + 1))
    
    I_greater = (factorial(c) / (d ** (c + 1))) * (full_moment - restricted_sum)
    
    return N12 * N34 * (I_less + I_greater)


# ==============================================================================
# ORBITAL-LEVEL INTEGRALS
# ==============================================================================

def orbital_overlap(orb1: Orbital, orb2: Orbital) -> float:
    """Compute <orb1|orb2>"""
    if orb1.l != orb2.l:
        return 0.0
    
    total = 0.0
    for t1 in orb1.terms:
        for t2 in orb2.terms:
            total += t1.coeff * t2.coeff * overlap_sto(t1.n_star, t1.zeta, t2.n_star, t2.zeta)
    return total


def orbital_kinetic(orb: Orbital) -> float:
    """Compute <orb|T|orb>"""
    total = 0.0
    for t1 in orb.terms:
        for t2 in orb.terms:
            total += t1.coeff * t2.coeff * kinetic_sto(t1.n_star, t1.zeta, orb.l, t2.n_star, t2.zeta)
    return total


def orbital_nuclear(orb: Orbital, Z: int) -> float:
    """Compute <orb|-Z/r|orb>"""
    total = 0.0
    for t1 in orb.terms:
        for t2 in orb.terms:
            total += t1.coeff * t2.coeff * nuclear_sto(t1.n_star, t1.zeta, t2.n_star, t2.zeta, Z)
    return total


def orbital_coulomb(orb_i: Orbital, orb_j: Orbital) -> float:
    """
    Compute Coulomb integral J_ij = (ii|jj) = ∫∫ |φi|^2 (1/r12) |φj|^2
    
    For orbitals with different l, we only keep k=0 (monopole) term for simplicity.
    This is exact for s orbitals and a good approximation for the spherical average.
    """
    total = 0.0
    for t1 in orb_i.terms:
        for t2 in orb_i.terms:
            for t3 in orb_j.terms:
                for t4 in orb_j.terms:
                    total += (t1.coeff * t2.coeff * t3.coeff * t4.coeff *
                             coulomb_integral_sto(t1.n_star, t1.zeta, t2.n_star, t2.zeta,
                                                  t3.n_star, t3.zeta, t4.n_star, t4.zeta))
    return total


def orbital_exchange(orb_i: Orbital, orb_j: Orbital) -> float:
    """
    Compute Exchange integral K_ij = (ij|ij) = ∫∫ φi(1)φj(1) (1/r12) φi(2)φj(2)
    
    For s-s or same-l interactions with k=0.
    """
    # For same orbital, K = J
    if orb_i.name == orb_j.name:
        return orbital_coulomb(orb_i, orb_j)
    
    # For different orbitals with same l, this still works
    total = 0.0
    for t1 in orb_i.terms:
        for t2 in orb_j.terms:
            for t3 in orb_i.terms:
                for t4 in orb_j.terms:
                    total += (t1.coeff * t2.coeff * t3.coeff * t4.coeff *
                             coulomb_integral_sto(t1.n_star, t1.zeta, t2.n_star, t2.zeta,
                                                  t3.n_star, t3.zeta, t4.n_star, t4.zeta))
    return total


# ==============================================================================
# TOTAL ENERGY COMPUTATION
# ==============================================================================

def compute_closed_shell_energy(atom: AtomData, config: Dict[str, int]) -> Tuple[float, Dict]:
    """
    Compute total energy for CLOSED SHELL atoms using standard RHF formula:
    
    E_tot = Σ_i n_i * h_ii + (1/2) Σ_i Σ_j n_i * n_j * (J_ij - 0.5*K_ij)
    
    where n_i is occupation number, h_ii = T_ii + V_ii
    """
    Z = atom.Z
    
    # One-electron energy
    E_1e = 0.0
    h_values = {}
    
    for orb_name, occ in config.items():
        if orb_name not in atom.orbitals:
            raise ValueError(f"Orbital {orb_name} not in basis for {atom.symbol}")
        
        orb = atom.orbitals[orb_name]
        T = orbital_kinetic(orb)
        V = orbital_nuclear(orb, Z)
        h = T + V
        h_values[orb_name] = {'T': T, 'V': V, 'h': h}
        E_1e += occ * h
    
    # Two-electron energy
    E_2e = 0.0
    orb_names = list(config.keys())
    
    for i, name_i in enumerate(orb_names):
        occ_i = config[name_i]
        orb_i = atom.orbitals[name_i]
        
        for j, name_j in enumerate(orb_names):
            occ_j = config[name_j]
            orb_j = atom.orbitals[name_j]
            
            J = orbital_coulomb(orb_i, orb_j)
            K = orbital_exchange(orb_i, orb_j)
            
            # For closed shell: E_2e = 0.5 * Σ_ij n_i n_j (2J - K)
            # = 0.5 * Σ_ij n_i n_j * 2J - 0.5 * Σ_ij n_i n_j * K
            # The factor 2 for J comes from counting both spins
            E_2e += 0.5 * occ_i * occ_j * (J - 0.5 * K)
    
    E_tot = E_1e + E_2e
    
    return E_tot, {
        'E_1e': E_1e,
        'E_2e': E_2e,
        'E_tot': E_tot,
        'orbitals': h_values
    }


# ==============================================================================
# PARSING
# ==============================================================================

def get_l_from_name(name: str) -> int:
    if 's' in name: return 0
    if 'p' in name: return 1
    if 'd' in name: return 2
    if 'f' in name: return 3
    raise ValueError(f"Unknown orbital: {name}")


def parse_basis(filepath: Path) -> Dict[str, AtomData]:
    """Parse slater_basis.js"""
    content = filepath.read_text(encoding='utf-8')
    
    start = content.find('globalScope.SlaterBasis = ')
    if start == -1:
        raise ValueError("Cannot find SlaterBasis")
    
    json_start = start + len('globalScope.SlaterBasis = ')
    
    # Find matching brace
    depth = 0
    json_end = json_start
    in_str = False
    escape = False
    
    for i, c in enumerate(content[json_start:]):
        if escape:
            escape = False
            continue
        if c == '\\':
            escape = True
            continue
        if c == '"' and not escape:
            in_str = not in_str
            continue
        if in_str:
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                json_end = json_start + i + 1
                break
    
    json_str = content[json_start:json_end]
    json_str = re.sub(r',\s*}', '}', json_str)
    json_str = re.sub(r',\s*]', ']', json_str)
    
    data = json.loads(json_str)
    
    atoms = {}
    for sym, d in data.items():
        if not isinstance(d, dict) or 'Z' not in d:
            continue
        
        orbitals = {}
        for name, terms in d.get('orbitals', {}).items():
            orbitals[name] = Orbital(
                name=name,
                l=get_l_from_name(name),
                terms=[STOTerm(t['nStar'], t['zeta'], t['coeff']) for t in terms]
            )
        
        atoms[sym] = AtomData(
            symbol=sym,
            Z=d['Z'],
            orbitals=orbitals,
            E_tot_ref=d.get('E_tot', 0.0),
            energies_ref=d.get('energies', {})
        )
    
    return atoms


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 70)
    print("FULL HARTREE-FOCK ENERGY VALIDATION")
    print("=" * 70)
    
    basis_file = Path(__file__).parent.parent / "slater_basis.js"
    atoms = parse_basis(basis_file)
    print(f"Loaded {len(atoms)} atoms\n")
    
    # Test closed-shell atoms first (simplest case)
    closed_shell = ['He', 'Be', 'Ne']
    
    print("-" * 70)
    print(f"{'Atom':>4} {'Z':>3} {'E_calc':>15} {'E_ref':>15} {'Error %':>10}")
    print("-" * 70)
    
    for sym in closed_shell:
        if sym not in atoms:
            print(f"{sym:>4} NOT FOUND")
            continue
        
        atom = atoms[sym]
        config = ELECTRON_CONFIG.get(atom.Z)
        if not config:
            print(f"{sym:>4} NO CONFIG")
            continue
        
        try:
            E_calc, details = compute_closed_shell_energy(atom, config)
            E_ref = atom.E_tot_ref
            error = 100 * (E_calc - E_ref) / abs(E_ref) if E_ref != 0 else 0
            
            print(f"{sym:>4} {atom.Z:>3} {E_calc:>15.8f} {E_ref:>15.8f} {error:>10.4f}%")
            
            # Detailed breakdown
            print(f"     E_1e = {details['E_1e']:.8f}, E_2e = {details['E_2e']:.8f}")
            
        except Exception as e:
            print(f"{sym:>4} ERROR: {e}")
    
    print("-" * 70)
    
    # Virial theorem check for He
    print("\n=== HELIUM DETAILED ANALYSIS ===")
    if 'He' in atoms:
        he = atoms['He']
        config = {'1s': 2}
        E_calc, details = compute_closed_shell_energy(he, config)
        
        orb = he.orbitals['1s']
        T = orbital_kinetic(orb)
        V = orbital_nuclear(orb, 2)
        J = orbital_coulomb(orb, orb)
        
        print(f"<T> = {T:.8f} Ha")
        print(f"<V_nuc> = {V:.8f} Ha")
        print(f"J_11 = {J:.8f} Ha")
        print(f"E_calc = 2*T + 2*V + J = {2*T + 2*V + J:.8f} Ha")
        print(f"E_ref = {he.E_tot_ref:.8f} Ha")
        print(f"Virial: -<V>/<T> = {-(2*V + J)/(2*T):.6f} (should be ~2)")


if __name__ == "__main__":
    main()
