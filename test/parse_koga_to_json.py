
import os
import json
import re

elements = {
    "He": "he",
    "Li": "li",
    "Be": "be",
    "B": "b",
    "C": "c",
    "N": "n",
    "O": "o",
    "F": "f",
    "Ne": "ne"
}

base_path = "AtomicOrbitals-master/literature_data/k99l/neutral"

def parse_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    orbitals_data = {}
    current_block_type = None # S, P, D, F
    orbital_labels = []
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
            
        if "BASIS/ORB.ENERGY" in line:
            # This line contains orbital energies. 
            # The line BEFORE this (or 2 lines before) contained the labels
            # Let's look back to find the header line "S   1S   2S"
            # It seems the header is a few lines up, after "ORBITAL ENERGIES..."
            # Let's search for the line starting with S, P, D, F indented
            
            # Actually, let's scan for the block headers
            # The format is:
            #         S                    1S             2S 
            #   BASIS/ORB.ENERGY       -2.4777413     -0.1963228
            
            # So we look for the line starting with S, P, D, F
            pass

        # Check for block start
        parts = line.split()
        if len(parts) > 0 and parts[0] in ['S', 'P', 'D', 'F']:
            current_block_type = parts[0]
            # The rest of the line are orbital names, e.g. "1S", "2S"
            orbital_labels = parts[1:]
            
            # Initialize lists for these orbitals
            for label in orbital_labels:
                orb_key = label.lower()
                if orb_key not in orbitals_data:
                    orbitals_data[orb_key] = []
            
            i += 1
            # Next line is BASIS/ORB.ENERGY
            while i < len(lines) and "BASIS/ORB.ENERGY" not in lines[i]:
                i += 1
            i += 1 # Skip BASIS line
            if i < len(lines) and "CUSP" in lines[i]:
                i += 1 # Skip CUSP line
            
            # Now read coefficients
            while i < len(lines):
                line = lines[i].strip()
                if not line or line.split()[0] in ['S', 'P', 'D', 'F', 'ORBITAL']:
                    i -= 1 # Step back so outer loop handles it
                    break
                
                parts = line.split()
                # Format: Label Zeta Coeff1 Coeff2 ...
                # e.g. 1S 10.335672 0.0014270 0.0002728
                
                # Check if valid line
                try:
                    # Sometimes the label has 'S', 'P' etc.
                    # Koga format: "1S", "2S" etc.
                    n_label = parts[0] # e.g. "1S"
                    # n is the first character
                    n = int(n_label[0])
                    zeta = float(parts[1])
                    coeffs = [float(x) for x in parts[2:]]
                    
                    for idx, coeff in enumerate(coeffs):
                        if idx < len(orbital_labels):
                            orb_name = orbital_labels[idx].lower()
                            orbitals_data[orb_name].append({
                                "n": n,
                                "zeta": zeta,
                                "coeff": coeff
                            })
                except (ValueError, IndexError):
                    pass
                
                i += 1
        i += 1

    return orbitals_data

result = {}
for el, filename in elements.items():
    path = os.path.join(base_path, filename)
    if os.path.exists(path):
        result[el] = {"orbitals": parse_file(path)}

print(json.dumps(result, indent=2))
