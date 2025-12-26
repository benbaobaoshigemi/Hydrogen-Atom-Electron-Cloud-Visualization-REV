# -*- coding: utf-8 -*-
"""
Roothaan Coupling Coefficients Validation Script
=================================================

PURPOSE:
This script validates the Roothaan coupling coefficients provided by an external advisor
by computing atomic total energies and comparing them against Bunge (1993) reference values.

METHODOLOGY:
1. Parse STO basis from slater_basis.js
2. Compute all required integrals (T, V_nuc, J, K) using analytical STO formulas
3. Assemble total energy using Roothaan average-energy formalism
4. Compare with E_tot from Bunge (1993)

CRITICAL CONSTRAINTS:
- NO fallback logic. If computation fails, raise an error.
- NO approximations. All integrals must be computed analytically.
- NO invented formulas. Only use textbook-verified equations.

Author: Antigravity Agent (validation mode)
Date: 2025-12-26
"""

import json
import re
import math
from pathlib import Path
from typing import Dict, List, Tuple, Any
from dataclasses import dataclass

# ==============================================================================
# ADVISOR'S ROOTHAAN COUPLING COEFFICIENTS (UNVERIFIED - TO BE TESTED)
# ==============================================================================
# Format: {Z: {'open': [orbital_keys], 'c': {(orb_i, orb_j): [2*alpha, beta]}}}
# Energy formula: E_open = f_i * f_j * (2*alpha * J - beta * K)

ATOM_ROOTHAAN_COEFFS = {
    # --- Closed shells (trivial) ---
    2:  {'term': '1S', 'open': [], 'c': {}},  # He
    4:  {'term': '1S', 'open': [], 'c': {}},  # Be
    10: {'term': '1S', 'open': [], 'c': {}},  # Ne
    12: {'term': '1S', 'open': [], 'c': {}},  # Mg
    18: {'term': '1S', 'open': [], 'c': {}},  # Ar
    20: {'term': '1S', 'open': [], 'c': {}},  # Ca
    30: {'term': '1S', 'open': [], 'c': {}},  # Zn
    36: {'term': '1S', 'open': [], 'c': {}},  # Kr
    
    # --- s-block open shells (Single s electron: 2S term) ---
    1:  {'term': '2S', 'open': ['1s'], 'c': {('1s','1s'): [0.0, 0.0]}},  # H
    3:  {'term': '2S', 'open': ['2s'], 'c': {('2s','2s'): [0.0, 0.0]}},  # Li
    11: {'term': '2S', 'open': ['3s'], 'c': {('3s','3s'): [0.0, 0.0]}},  # Na
    19: {'term': '2S', 'open': ['4s'], 'c': {('4s','4s'): [0.0, 0.0]}},  # K
    
    # --- p-block ---
    # p^1 (2P): B, Al, Ga
    5:  {'term': '2P', 'open': ['2p'], 'c': {('2p','2p'): [0.0, 0.0]}},
    13: {'term': '2P', 'open': ['3p'], 'c': {('3p','3p'): [0.0, 0.0]}},
    31: {'term': '2P', 'open': ['4p'], 'c': {('4p','4p'): [0.0, 0.0]}},
    
    # p^2 (3P): C, Si, Ge
    6:  {'term': '3P', 'open': ['2p'], 'c': {('2p','2p'): [0.5, 1.5]}},
    14: {'term': '3P', 'open': ['3p'], 'c': {('3p','3p'): [0.5, 1.5]}},
    32: {'term': '3P', 'open': ['4p'], 'c': {('4p','4p'): [0.5, 1.5]}},
    
    # p^3 (4S): N, P, As
    7:  {'term': '4S', 'open': ['2p'], 'c': {('2p','2p'): [1.0, 2.0]}},
    15: {'term': '4S', 'open': ['3p'], 'c': {('3p','3p'): [1.0, 2.0]}},
    33: {'term': '4S', 'open': ['4p'], 'c': {('4p','4p'): [1.0, 2.0]}},
    
    # p^4 (3P): O, S, Se
    8:  {'term': '3P', 'open': ['2p'], 'c': {('2p','2p'): [0.5, 1.5]}},
    16: {'term': '3P', 'open': ['3p'], 'c': {('3p','3p'): [0.5, 1.5]}},
    34: {'term': '3P', 'open': ['4p'], 'c': {('4p','4p'): [0.5, 1.5]}},
    
    # p^5 (2P): F, Cl, Br
    9:  {'term': '2P', 'open': ['2p'], 'c': {('2p','2p'): [0.0, 0.0]}},
    17: {'term': '2P', 'open': ['3p'], 'c': {('3p','3p'): [0.0, 0.0]}},
    35: {'term': '2P', 'open': ['4p'], 'c': {('4p','4p'): [0.0, 0.0]}},
    
    # --- d-block ---
    # d^1 (2D): Sc
    21: {'term': '2D', 'open': ['3d'], 'c': {('3d','3d'): [0.0, 0.0]}},
    # d^2 (3F): Ti
    22: {'term': '3F', 'open': ['3d'], 'c': {('3d','3d'): [0.5, 1.5]}},
    # d^3 (4F): V
    23: {'term': '4F', 'open': ['3d'], 'c': {('3d','3d'): [2/3, 5/3]}},
    
    # Cr: 4s^1 3d^5 (7S) - DOUBLE OPEN SHELL
    24: {'term': '7S', 'open': ['3d', '4s'], 'c': {
        ('3d','3d'): [1.0, 2.0],  # d^5 self
        ('4s','4s'): [0.0, 0.0],  # s^1 self
        ('3d','4s'): [1.0, 1.0],  # Cross-term
        ('4s','3d'): [1.0, 1.0],
    }},
    
    # Mn: 4s^2 3d^5 (6S) - 4s closed
    25: {'term': '6S', 'open': ['3d'], 'c': {('3d','3d'): [1.0, 2.0]}},
    
    # Fe: d^6 (5D)
    26: {'term': '5D', 'open': ['3d'], 'c': {('3d','3d'): [0.75, 1.75]}},
    
    # Co: d^7 (4F)
    27: {'term': '4F', 'open': ['3d'], 'c': {('3d','3d'): [2/3, 5/3]}},
    
    # Ni: d^8 (3F)
    28: {'term': '3F', 'open': ['3d'], 'c': {('3d','3d'): [0.5, 1.5]}},
    
    # Cu: 4s^1 3d^10 (2S) - 3d closed
    29: {'term': '2S', 'open': ['4s'], 'c': {('4s','4s'): [0.0, 0.0]}},
}


# ==============================================================================
# STO BASIS DATA STRUCTURES
# ==============================================================================

@dataclass
class STOTerm:
    """Single STO primitive: N * r^(n-1) * exp(-zeta*r)"""
    n_star: int
    zeta: float
    coeff: float


@dataclass
class Orbital:
    """Atomic orbital as linear combination of STOs"""
    name: str  # e.g., '1s', '2p', '3d'
    l: int     # Angular momentum quantum number
    terms: List[STOTerm]


@dataclass
class AtomData:
    """Complete atomic data"""
    symbol: str
    Z: int
    orbitals: Dict[str, Orbital]
    E_tot_ref: float  # Bunge reference energy
    energies_ref: Dict[str, float]  # Orbital energies


# ==============================================================================
# PARSING FUNCTIONS
# ==============================================================================

def get_l_from_orbital_name(name: str) -> int:
    """Extract angular momentum from orbital name"""
    if 's' in name:
        return 0
    elif 'p' in name:
        return 1
    elif 'd' in name:
        return 2
    elif 'f' in name:
        return 3
    else:
        raise ValueError(f"Unknown orbital type: {name}")


def parse_slater_basis_js(filepath: Path) -> Dict[str, AtomData]:
    """Parse slater_basis.js and extract atomic data using JSON parsing"""
    content = filepath.read_text(encoding='utf-8')
    
    # Find the JSON object within the JavaScript file
    # Look for: globalScope.SlaterBasis = { ... };
    start_marker = 'globalScope.SlaterBasis = '
    start_idx = content.find(start_marker)
    if start_idx == -1:
        raise ValueError("Could not find SlaterBasis object in file")
    
    json_start = start_idx + len(start_marker)
    
    # Find the matching closing brace
    brace_count = 0
    json_end = json_start
    in_string = False
    escape_next = False
    
    for i, c in enumerate(content[json_start:]):
        if escape_next:
            escape_next = False
            continue
        if c == '\\':
            escape_next = True
            continue
        if c == '"' and not escape_next:
            in_string = not in_string
            continue
        if in_string:
            continue
        if c == '{':
            brace_count += 1
        elif c == '}':
            brace_count -= 1
            if brace_count == 0:
                json_end = json_start + i + 1
                break
    
    json_str = content[json_start:json_end]
    
    # Clean up JavaScript-specific syntax that's not valid JSON
    # 1. Remove trailing commas before closing braces/brackets
    import re as re_module
    json_str = re_module.sub(r',\s*}', '}', json_str)
    json_str = re_module.sub(r',\s*]', ']', json_str)
    
    # Parse the JSON
    try:
        data = json.loads(json_str)
    except json.JSONDecodeError as e:
        # Save failed JSON for debugging
        debug_path = filepath.parent / "_debug_json_parse_fail.txt"
        with open(debug_path, 'w', encoding='utf-8') as f:
            f.write(f"Error at position {e.pos}: {e.msg}\n")
            f.write(f"Context: ...{json_str[max(0,e.pos-50):e.pos+50]}...\n")
        raise ValueError(f"Failed to parse JSON: {e}")
    
    # Convert to AtomData objects
    atoms = {}
    for symbol, atom_dict in data.items():
        # Skip non-element entries
        if not isinstance(atom_dict, dict) or 'Z' not in atom_dict:
            continue
        
        orbitals = {}
        for orb_name, terms_list in atom_dict.get('orbitals', {}).items():
            terms = [
                STOTerm(
                    n_star=term['nStar'],
                    zeta=term['zeta'],
                    coeff=term['coeff']
                )
                for term in terms_list
            ]
            orbitals[orb_name] = Orbital(
                name=orb_name,
                l=get_l_from_orbital_name(orb_name),
                terms=terms
            )
        
        atoms[symbol] = AtomData(
            symbol=symbol,
            Z=atom_dict['Z'],
            orbitals=orbitals,
            E_tot_ref=atom_dict.get('E_tot', 0.0),
            energies_ref=atom_dict.get('energies', {})
        )
    
    return atoms


# ==============================================================================
# ANALYTICAL STO INTEGRALS
# ==============================================================================

def factorial(n: int) -> int:
    """Compute n!"""
    if n < 0:
        raise ValueError(f"Factorial undefined for negative n: {n}")
    if n <= 1:
        return 1
    result = 1
    for i in range(2, n + 1):
        result *= i
    return result


def double_factorial(n: int) -> int:
    """Compute n!! = n * (n-2) * (n-4) * ..."""
    if n <= 0:
        return 1
    result = 1
    while n > 0:
        result *= n
        n -= 2
    return result


def sto_normalization(n: int, zeta: float) -> float:
    """
    Normalization constant for STO: N = (2*zeta)^n * sqrt(2*zeta / (2n)!)
    
    The radial STO is: R(r) = N * r^(n-1) * exp(-zeta * r)
    """
    power_term = (2 * zeta) ** n
    sqrt_term = math.sqrt(2 * zeta / factorial(2 * n))
    return power_term * sqrt_term


def overlap_integral_sto(n1: int, z1: float, n2: int, z2: float) -> float:
    """
    Compute overlap integral <phi_1|phi_2> for two STOs (same l).
    
    S = ∫ R1(r) R2(r) r^2 dr = N1 N2 ∫ r^(n1+n2) exp(-(z1+z2)*r) dr
      = N1 N2 * (n1+n2)! / (z1+z2)^(n1+n2+1)
    """
    N1 = sto_normalization(n1, z1)
    N2 = sto_normalization(n2, z2)
    alpha = z1 + z2
    power = n1 + n2  # exponent of r in integrand (after r^2 from measure)
    
    integral = factorial(power) / (alpha ** (power + 1))
    return N1 * N2 * integral


def kinetic_integral_sto(n1: int, z1: float, l1: int, n2: int, z2: float, l2: int) -> float:
    """
    Compute kinetic energy integral <phi_1|-0.5*nabla^2|phi_2> for STOs.
    
    Using the radial kinetic energy operator:
    T = -0.5 * (d^2/dr^2 + 2/r * d/dr - l(l+1)/r^2)
    
    For STO: chi = r^(n-1) exp(-z*r)
    T chi = -0.5 * [A*r^(n-3) + B*r^(n-2) + C*r^(n-1)] exp(-z*r)
    
    where:
    A = n(n-1) - l(l+1)
    B = -2*z*n
    C = z^2
    """
    if l1 != l2:
        return 0.0  # Orthogonal by angular momentum
    
    l = l1
    N1 = sto_normalization(n1, z1)
    N2 = sto_normalization(n2, z2)
    alpha = z1 + z2
    
    # Coefficients for T|phi_2>
    A = n2 * (n2 - 1) - l * (l + 1)
    B = -2 * z2 * n2
    C = z2 * z2
    
    # Integral: <phi_1|T|phi_2> = -0.5 * N1 * N2 * ∫ r^(n1-1) * [A*r^(n2-3) + B*r^(n2-2) + C*r^(n2-1)] * exp(-alpha*r) * r^2 dr
    # = -0.5 * N1 * N2 * [A * I(n1+n2-2) + B * I(n1+n2-1) + C * I(n1+n2)]
    # where I(k) = ∫ r^k exp(-alpha*r) dr = k! / alpha^(k+1)
    
    def radial_moment(k):
        if k < 0:
            return 0.0  # These terms don't contribute for physical STOs
        return factorial(k) / (alpha ** (k + 1))
    
    k_base = n1 + n2
    integral = -0.5 * (
        A * radial_moment(k_base - 2) +
        B * radial_moment(k_base - 1) +
        C * radial_moment(k_base)
    )
    
    return N1 * N2 * integral


def nuclear_integral_sto(n1: int, z1: float, n2: int, z2: float, Z_nuc: int) -> float:
    """
    Compute nuclear attraction integral <phi_1|-Z/r|phi_2> for STOs.
    
    V = ∫ R1(r) (-Z/r) R2(r) r^2 dr = -Z * N1 N2 ∫ r^(n1+n2-1) exp(-alpha*r) dr
      = -Z * N1 N2 * (n1+n2-1)! / alpha^(n1+n2)
    """
    N1 = sto_normalization(n1, z1)
    N2 = sto_normalization(n2, z2)
    alpha = z1 + z2
    k = n1 + n2 - 1  # power of r after -1 from 1/r
    
    if k < 0:
        raise ValueError(f"Invalid radial power {k} in nuclear integral")
    
    integral = factorial(k) / (alpha ** (k + 1))
    return -Z_nuc * N1 * N2 * integral


def coulomb_integral_sto_k0(n1: int, z1: float, n2: int, z2: float,
                             n3: int, z3: float, n4: int, z4: float) -> float:
    """
    Compute k=0 (monopole) Coulomb integral: (12|34) = <phi_1 phi_2|1/r_12|phi_3 phi_4>
    
    For spherically symmetric (k=0) case with same-l orbitals:
    J = ∫∫ rho_12(r1) * (1/r_>) * rho_34(r2) * r1^2 r2^2 dr1 dr2
    
    Using the Laplace expansion for 1/r_>:
    1/r_> = 1/r_< for r1 < r2, 1/r_1 for r1 > r2
    
    This is the most complex integral. For STO, it involves 
    incomplete gamma functions or can be computed via recursion.
    
    Simplified formula for J_0 (k=0 Slater integral):
    """
    # This is a placeholder - the actual implementation requires
    # proper handling of the r_< and r_> terms
    # For now, we'll use numerical integration as a fallback check
    
    # For a proper implementation, see:
    # Harris, F.E., "Analytic Evaluation of Two-Center STO Integrals"
    # or use the Yk(r) potential method from the advisor's suggestion
    
    raise NotImplementedError(
        "Full analytical two-electron integrals not implemented. "
        "This validation requires proper Yk(r) potential method."
    )


# ==============================================================================
# SIMPLIFIED VALIDATION: ONE-ELECTRON INTEGRALS ONLY
# ==============================================================================

def compute_one_electron_energy(atom: AtomData) -> Tuple[float, Dict[str, float]]:
    """
    Compute one-electron energy: E_1e = sum_i n_i * (T_i + V_nuc_i)
    
    This is a PARTIAL validation - it only checks the one-electron part.
    The full energy requires two-electron integrals.
    """
    Z = atom.Z
    
    # Determine electron configuration
    # This is a simplification - should come from data
    electron_counts = {
        '1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '3d': 10, '4s': 2, '4p': 6, '4d': 10, '4f': 14
    }
    
    total_1e = 0.0
    orbital_1e = {}
    
    for orb_name, orbital in atom.orbitals.items():
        # Get occupation (simplified)
        max_occ = 2 * (2 * orbital.l + 1)
        # Determine actual occupation from Z
        # This is a rough approximation
        
        # Compute <T> + <V_nuc> for this orbital
        T_orb = 0.0
        V_orb = 0.0
        
        for i, term_i in enumerate(orbital.terms):
            for j, term_j in enumerate(orbital.terms):
                # Kinetic
                T_ij = kinetic_integral_sto(
                    term_i.n_star, term_i.zeta, orbital.l,
                    term_j.n_star, term_j.zeta, orbital.l
                )
                T_orb += term_i.coeff * term_j.coeff * T_ij
                
                # Nuclear
                V_ij = nuclear_integral_sto(
                    term_i.n_star, term_i.zeta,
                    term_j.n_star, term_j.zeta,
                    Z
                )
                V_orb += term_i.coeff * term_j.coeff * V_ij
        
        # Also check normalization
        S_orb = 0.0
        for i, term_i in enumerate(orbital.terms):
            for j, term_j in enumerate(orbital.terms):
                S_ij = overlap_integral_sto(
                    term_i.n_star, term_i.zeta,
                    term_j.n_star, term_j.zeta
                )
                S_orb += term_i.coeff * term_j.coeff * S_ij
        
        orbital_1e[orb_name] = {
            'T': T_orb,
            'V': V_orb,
            'E_1e': T_orb + V_orb,
            'S': S_orb  # Should be ~1.0 for normalized orbital
        }
    
    return total_1e, orbital_1e


# ==============================================================================
# MAIN VALIDATION
# ==============================================================================

def main():
    print("=" * 70)
    print("ROOTHAAN COUPLING COEFFICIENTS VALIDATION")
    print("=" * 70)
    print()
    
    # Find the slater_basis.js file
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    basis_file = project_root / "slater_basis.js"
    
    if not basis_file.exists():
        raise FileNotFoundError(f"Cannot find basis file: {basis_file}")
    
    print(f"Parsing basis file: {basis_file}")
    atoms = parse_slater_basis_js(basis_file)
    print(f"Loaded {len(atoms)} atoms\n")
    
    # Test atoms with known Roothaan coefficients
    test_atoms = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
    
    print("-" * 70)
    print(f"{'Atom':>4} {'Z':>3} {'E_tot(ref)':>15} {'<T>+<V>':>15} {'S(1s)':>10} {'S(2s)':>10}")
    print("-" * 70)
    
    for symbol in test_atoms:
        if symbol not in atoms:
            print(f"{symbol:>4} -- NOT FOUND --")
            continue
        
        atom = atoms[symbol]
        total_1e, orb_1e = compute_one_electron_energy(atom)
        
        # Sum orbital contributions
        sum_TV = sum(d['E_1e'] for d in orb_1e.values())
        
        # Normalization check
        S_1s = orb_1e.get('1s', {}).get('S', 0)
        S_2s = orb_1e.get('2s', {}).get('S', 0)
        
        print(f"{symbol:>4} {atom.Z:>3} {atom.E_tot_ref:>15.8f} {sum_TV:>15.8f} {S_1s:>10.6f} {S_2s:>10.6f}")
    
    print("-" * 70)
    print()
    print("NOTES:")
    print("- <T>+<V> is the one-electron energy sum (NOT total energy)")
    print("- S values should be ~1.0 for normalized orbitals")
    print("- Full validation requires two-electron integrals (J, K)")
    print()
    
    # Detailed breakdown for Carbon (the user's target)
    print("=" * 70)
    print("DETAILED BREAKDOWN: CARBON (Z=6)")
    print("=" * 70)
    
    if 'C' in atoms:
        carbon = atoms['C']
        _, orb_1e = compute_one_electron_energy(carbon)
        
        for orb_name, data in orb_1e.items():
            print(f"\n{orb_name} orbital:")
            print(f"  Kinetic energy <T>:     {data['T']:>15.8f} Ha")
            print(f"  Nuclear energy <V>:     {data['V']:>15.8f} Ha")
            print(f"  One-electron <T+V>:     {data['E_1e']:>15.8f} Ha")
            print(f"  Normalization <psi|psi>: {data['S']:>15.8f} (should be 1.0)")
            
            if orb_name in carbon.energies_ref:
                ref_e = carbon.energies_ref[orb_name]
                print(f"  Reference orbital ε:    {ref_e:>15.8f} Ha")
        
        print(f"\nReference total energy: {carbon.E_tot_ref:>15.8f} Ha")
        print(f"(Full calculation requires J and K integrals)")
    
    print()
    print("=" * 70)
    print("VALIDATION STATUS: PARTIAL (one-electron only)")
    print("=" * 70)
    print()
    print("NEXT STEPS:")
    print("1. Implement Yk(r) potential for two-electron integrals")
    print("2. Validate Roothaan coefficients against full E_tot")
    print("3. Compare with Bunge (1993) orbital energies")


if __name__ == "__main__":
    main()
