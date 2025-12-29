# -*- coding: utf-8 -*-
"""
PRODUCTION: Generate Vee Cache for All Atoms (H-Kr)
Exports to vee_cache.js for frontend use.
"""
import numpy as np
from scipy.special import gamma, gammainc, gammaincc, factorial
import json, re
from pathlib import Path

# ==================== MATH KERNELS ====================
def sto_normalization(n, zeta):
    return (2 * zeta) ** n * np.sqrt(2 * zeta / float(factorial(2 * n)))

def lower_incomplete_gamma(s, x):
    return gamma(s) * gammainc(s, x)

def upper_incomplete_gamma(s, x):
    return gamma(s) * gammaincc(s, x)

def compute_Yk_analytic(basis_a, basis_b, k, r_grid):
    result = np.zeros_like(r_grid)
    for term_a in basis_a:
        na, za, ca = term_a['nStar'], term_a['zeta'], term_a['coeff']
        Na = sto_normalization(na, za)
        for term_b in basis_b:
            nb, zb, cb = term_b['nStar'], term_b['zeta'], term_b['coeff']
            Nb = sto_normalization(nb, zb)
            pref = ca * cb * Na * Nb
            N_power = na + nb
            Z_exp = za + zb
            inner_power = k + N_power + 1
            term1 = (1.0 / np.power(r_grid, k+1)) * (1.0 / Z_exp**inner_power) * lower_incomplete_gamma(inner_power, Z_exp * r_grid)
            outer_power = N_power - k
            if outer_power > 0:
                term2 = np.power(r_grid, k) * (1.0 / Z_exp**outer_power) * upper_incomplete_gamma(outer_power, Z_exp * r_grid)
            else:
                term2 = 0.0
            result += pref * (term1 + term2)
    return result

EXCHANGE_COEFFS = {
    (0, 0): [(0, 1.0)], (0, 1): [(1, 1/3)], (1, 0): [(1, 1/3)],
    (1, 1): [(0, 1/3), (2, 2/15)], (0, 2): [(2, 1/5)], (2, 0): [(2, 1/5)],
    (1, 2): [(1, 2/15), (3, 3/35)], (2, 1): [(1, 2/15), (3, 3/35)],
    (2, 2): [(0, 1/5), (2, 2/35), (4, 2/35)],
}

def get_l(orb_name):
    if 's' in orb_name: return 0
    if 'p' in orb_name: return 1
    if 'd' in orb_name: return 2
    return 0

def parse_slater_basis(filepath):
    content = Path(filepath).read_text(encoding='utf-8')
    marker = 'globalScope.SlaterBasis = '
    start_idx = content.find(marker)
    json_start = start_idx + len(marker)
    depth = 0
    json_end = -1
    for i, c in enumerate(content[json_start:]):
        if c == '{': depth += 1
        elif c == '}': 
            depth -= 1
            if depth == 0: json_end = json_start + i + 1; break
    json_str = content[json_start:json_end]
    json_str = re.sub(r'//.*', '', json_str)
    json_str = re.sub(r',(\s*[}\]])', r'\1', json_str)
    return json.loads(json_str)

# Full H-Kr configs
ELECTRON_CONFIGS = {
    "H":  {"1s": 1}, "He": {"1s": 2},
    "Li": {"1s": 2, "2s": 1}, "Be": {"1s": 2, "2s": 2},
    "B":  {"1s": 2, "2s": 2, "2p": 1}, "C":  {"1s": 2, "2s": 2, "2p": 2},
    "N":  {"1s": 2, "2s": 2, "2p": 3}, "O":  {"1s": 2, "2s": 2, "2p": 4},
    "F":  {"1s": 2, "2s": 2, "2p": 5}, "Ne": {"1s": 2, "2s": 2, "2p": 6},
    "Na": {"1s": 2, "2s": 2, "2p": 6, "3s": 1}, "Mg": {"1s": 2, "2s": 2, "2p": 6, "3s": 2},
    "Al": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 1}, "Si": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 2},
    "P":  {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 3}, "S":  {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 4},
    "Cl": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 5}, "Ar": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6},
    "K":  {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "4s": 1},
    "Ca": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "4s": 2},
    "Sc": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 1, "4s": 2},
    "Ti": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 2, "4s": 2},
    "V":  {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 3, "4s": 2},
    "Cr": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 5, "4s": 1},
    "Mn": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 5, "4s": 2},
    "Fe": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 6, "4s": 2},
    "Co": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 7, "4s": 2},
    "Ni": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 8, "4s": 2},
    "Cu": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10, "4s": 1},
    "Zn": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10, "4s": 2},
    "Ga": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10, "4s": 2, "4p": 1},
    "Ge": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10, "4s": 2, "4p": 2},
    "As": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10, "4s": 2, "4p": 3},
    "Se": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10, "4s": 2, "4p": 4},
    "Br": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10, "4s": 2, "4p": 5},
    "Kr": {"1s": 2, "2s": 2, "2p": 6, "3s": 2, "3p": 6, "3d": 10, "4s": 2, "4p": 6},
}

def get_spin_factor(orb_name, n_occupancy):
    """
    Calculate Spin Polarization Factor chi = (N_up^2 + N_down^2) / N^2
    Based on Hund's Rule (High Spin Ground State assumption).
    
    This corrects the Exchange interaction:
    - Closed Shells (Singlet): chi = 0.5 (Corrects the 'Ferromagnetic Error')
    - Open Shells (High Spin): chi = (N_up^2 + N_down^2)/N^2 (e.g. O 2p4 -> 0.625)
    """
    l = get_l(orb_name)
    capacity = 2 * (2 * l + 1)
    
    # Hund's Rule allocation
    half_cap = capacity // 2
    n_up = min(n_occupancy, half_cap)
    n_down = max(0, n_occupancy - half_cap)
    
    # chi = N_eff / N = (N_up^2 + N_down^2) / N^2
    return (float(n_up)**2 + float(n_down)**2) / (float(n_occupancy)**2)

def generate_vee_cache(atom_symbol, basis_data, config, r_grid):
    orbitals = basis_data['orbitals']
    cache = {}
    n_total = int(sum(config.values()))
    # For a neutral atom (or a configuration with N electrons), an electron outside the atom
    # must feel a net +1 charge: V_ee(r) -> (N-1)/r as r -> ∞ (exchange is short-ranged).
    a_tail = max(0, n_total - 1)
    for target_orb, target_basis in orbitals.items():
        if target_orb not in config: continue
        l_target = get_l(target_orb)
        Vee = np.zeros_like(r_grid)
        for source_orb, source_basis in orbitals.items():
            if source_orb not in config: continue
            n_source = config[source_orb]
            l_source = get_l(source_orb)
            
            # Direct Coulomb (Hartree Term): N * Y0(aa)
            Y0 = compute_Yk_analytic(source_basis, source_basis, 0, r_grid)
            Vee += n_source * Y0
            
            # Exchange Interaction (Fock Term)
            # HF Potential: V_eff = V_Coulomb - V_Exchange
            # V_Exchange occurs only between parallel spins.
            # We apply the Spin Polarization Factor to scale the "All-Parallel" geometrical term.
            key = (l_target, l_source)
            if key in EXCHANGE_COEFFS:
                # Calculate the rigorous spin factor for the source shell
                spin_factor = get_spin_factor(source_orb, n_source)
                
                for k, coeff in EXCHANGE_COEFFS[key]:
                    Yk = compute_Yk_analytic(target_basis, source_basis, k, r_grid)
                    Vee -= coeff * n_source * Yk * spin_factor

                # Enforce the physically required far-field asymptote:
                # For a neutral atom, Vee(r) -> (N-1)/r as r -> ∞ (exchange is short-ranged).
                # This prevents open-shell p/d subshell exchange-averaging artifacts from polluting
                # Z_eff and other long-range derived diagnostics.
        if a_tail > 0:
            tail_count = min(32, len(r_grid))
            tail_start = len(r_grid) - tail_count
            Vee[tail_start:] = (a_tail / r_grid[tail_start:])
                    
        cache[target_orb] = Vee.tolist()
    return cache

def generate_full_cache():
    print("Generating Full Vee Cache for H-Kr...")
    
    basis_path = "../slater_basis.js"
    if not Path(basis_path).exists(): basis_path = "slater_basis.js"
    data = parse_slater_basis(basis_path)
    
    # 700 point log-spaced grid, extended to 50.0 a.u.
    r_grid = np.geomspace(1e-4, 50.0, 700).tolist()
    
    full_cache = {
        "r_grid": r_grid,
        "atoms": {}
    }
    
    elements = list(ELECTRON_CONFIGS.keys())
    
    for i, atom in enumerate(elements):
        if atom not in data:
            print(f"  [SKIP] {atom}: Not in basis")
            continue
            
        cache = generate_vee_cache(atom, data[atom], ELECTRON_CONFIGS[atom], np.array(r_grid))
        full_cache["atoms"][atom] = cache
        print(f"  [{i+1}/{len(elements)}] {atom}: {len(cache)} orbitals")
    
    # Export as JS file
    output_path = Path("../vee_cache.js")
    
    js_content = f"""// Auto-generated Vee(r) cache for Hartree-Fock potential
// Generated by generate_production_cache.py
// Contains pre-computed electron-electron repulsion potentials for H-Kr

(typeof self !== 'undefined' ? self : window).VeeCache = {json.dumps(full_cache, indent=2)};
"""
    
    output_path.write_text(js_content, encoding='utf-8')
    
    file_size = output_path.stat().st_size / 1024
    print(f"\nExported to: {output_path.absolute()}")
    print(f"File size: {file_size:.1f} KB")
    print("DONE!")

if __name__ == "__main__":
    generate_full_cache()
