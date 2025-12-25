import os
import json
import re

SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
]

def parse_koga_file(path, symbol, Z):
    if not os.path.exists(path):
        return None
        
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    data = {"Z": Z, "name": symbol, "orbitals": {}}
    
    # Extract total energy
    e_match = re.search(r"E\s*=\s*(-?\d+\.\d+)", content)
    if e_match:
        data["E_tot"] = float(e_match.group(1))
    
    # Split into blocks of S, P, D, F
    blocks = re.split(r"\s+([SPDF])\s+", content)
    
    for b_idx in range(1, len(blocks), 2):
        am_char = blocks[b_idx]
        am_content = blocks[b_idx+1]
        lines = am_content.strip().split("\n")
        if not lines: continue
        
        # Header line contains orbital names like "1S 2S 3S 4S"
        # The line starting with BASIS/ORB.ENERGY actually has the energies
        orb_names = []
        orb_energies = []
        
        # Find the line starting with "S" or "P" etc. at the very top of the block
        # Wait, the split might have removed it.
        # Let's re-parse simply.
        pass
    
    # Using a more robust line-by-line parsing
    lines = content.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        parts = line.split()
        if parts and parts[0] in ["S", "P", "D", "F"] and len(parts) > 1:
            am_char = parts[0]
            orb_names = [n.lower() for n in parts[1:]]
            i += 1
            # Next line should be BASIS/ORB.ENERGY
            if i < len(lines) and "BASIS/ORB.ENERGY" in lines[i]:
                orb_energies = [float(e) for e in lines[i].replace("BASIS/ORB.ENERGY", "").split()]
                i += 1
                if i < len(lines) and "CUSP" in lines[i]:
                    i += 1
                
                # Basis functions
                basis_funcs = []
                while i < len(lines):
                    row = lines[i].split()
                    if not row or row[0] in ["S", "P", "D", "F"]:
                        break
                    if row[0][-1] == am_char:
                        n_star = int(row[0][:-1])
                        zeta = float(row[1])
                        coeffs = [float(c) for c in row[2:]]
                        basis_funcs.append({"nStar": n_star, "zeta": zeta, "coeffs": coeffs})
                        i += 1
                    else:
                        break
                
                # Map
                for j, orb_name in enumerate(orb_names):
                    orb_basis = []
                    for bf in basis_funcs:
                        if j < len(bf["coeffs"]):
                            orb_basis.append({
                                "nStar": bf["nStar"],
                                "zeta": bf["zeta"],
                                "coeff": bf["coeffs"][j]
                            })
                    data["orbitals"][orb_name] = {
                        "energy": orb_energies[j] if j < len(orb_energies) else 0,
                        "basis": orb_basis
                    }
                continue
        i += 1
    return data

def main():
    root = os.getcwd()
    k99_dir = os.path.join(root, 'AtomicOrbitals-master', 'literature_data', 'k99l', 'neutral')
    k00_dir = os.path.join(root, 'AtomicOrbitals-master', 'literature_data', 'k00heavy')
    
    final_db = {}
    for i, symbol in enumerate(SYMBOLS):
        Z = i + 1
        path = os.path.join(k99_dir, symbol.lower())
        if not os.path.exists(path):
            path = os.path.join(k00_dir, symbol.lower())
            
        if os.path.exists(path):
            try:
                elem_data = parse_koga_file(path, symbol, Z)
                if elem_data:
                    final_db[symbol] = elem_data
            except Exception as e:
                print(f"Error parsing {symbol}: {e}")
                
    with open('koga_basis_database.json', 'w', encoding='utf-8') as f:
        json.dump(final_db, f, indent=2)
    print(f"Extracted {len(final_db)} elements.")

if __name__ == "__main__":
    main()
