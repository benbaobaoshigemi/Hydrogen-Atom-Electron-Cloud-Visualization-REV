import re
import json

def parse_clementi_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    basis_data = {}
    current_atom = None
    current_block_type = None
    current_orbitals = []
    
    # regex patterns
    atom_header_re = re.compile(r'^\s*([A-Za-z]+(-[A-Za-z]+)?):Clementi-Roetti')
    orbital_def_re = re.compile(r'^\s*([SPD])\s*\{\s*(.*?)\s*\}') # Capture S { 1S 2S }
    
    atomic_numbers = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
        'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
        'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
        'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
        'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54
    }

    names = {
        'H': 'Hydrogen', 'He': 'Helium', 'Li': 'Lithium', 'Be': 'Beryllium', 'B': 'Boron', 'C': 'Carbon', 'N': 'Nitrogen', 'O': 'Oxygen', 'F': 'Fluorine', 'Ne': 'Neon',
        'Na': 'Sodium', 'Mg': 'Magnesium', 'Al': 'Aluminum', 'Si': 'Silicon', 'P': 'Phosphorus', 'S': 'Sulfur', 'Cl': 'Chlorine', 'Ar': 'Argon',
        'K': 'Potassium', 'Ca': 'Calcium', 'Sc': 'Scandium', 'Ti': 'Titanium', 'V': 'Vanadium', 'Cr': 'Chromium', 'Mn': 'Manganese', 'Fe': 'Iron', 'Co': 'Cobalt', 'Ni': 'Nickel', 'Cu': 'Copper', 'Zn': 'Zinc',
        'Ga': 'Gallium', 'Ge': 'Germanium', 'As': 'Arsenic', 'Se': 'Selenium', 'Br': 'Bromine', 'Kr': 'Krypton', 'Rb': 'Rubidium', 'Sr': 'Strontium', 'Y': 'Yttrium', 'Zr': 'Zirconium', 'Nb': 'Niobium', 'Mo': 'Molybdenum', 'Tc': 'Technetium', 'Ru': 'Ruthenium', 'Rh': 'Rhodium', 'Pd': 'Palladium', 'Ag': 'Silver', 'Cd': 'Cadmium', 'In': 'Indium', 'Sn': 'Tin', 'Sb': 'Antimony', 'Te': 'Tellurium', 'I': 'Iodine', 'Xe': 'Xenon'
    }

    # Temporary storage for current parsing context
    temp_orbitals_data = {} # Map "1S" -> list of primitives

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # New Atom Block
        atom_match = atom_header_re.match(line)
        if atom_match:
            # Save previous atom if exists
            if current_atom:
                # Add version if modified (e.g. H-normal) but let's stick to standard names
                # For basic mapping we filter for standard element symbols
                if current_atom in atomic_numbers:
                    basis_data[current_atom] = {
                        'Z': atomic_numbers[current_atom],
                        'name': names.get(current_atom, current_atom),
                        'groundState': current_config, # Extracted below
                        'orbitals': temp_orbitals_data
                    }
            
            current_atom = atom_match.group(1)
            temp_orbitals_data = {}
            # Next line usually contains config string like "1S(2)2S(1)"
            if i + 1 < len(lines):
                 # Try to grab config line (skipping potential blanks)
                 # Actually the file format is Atom:CR \n Config \n {
                 config_line = lines[i+1].strip()
                 # Remove comments like "! 2S ?"
                 if '!' in config_line:
                     config_part = config_line.split('!')[0].strip()
                 else:
                     config_part = config_line
                 current_config = config_part
            
            i += 1
            continue

        # Orbital Block Definition: S { 1S 2S }
        orb_match = orbital_def_re.match(line)
        if orb_match:
            orb_type = orb_match.group(1) # S, P, D
            orb_names = orb_match.group(2).split() # ['1S', '2S']
            current_orbitals = orb_names
            
            # Ensure orbital lists exist
            for orb in current_orbitals:
                if orb.lower() not in temp_orbitals_data: # Normalize keys to lowercase later? Or keep standard
                    temp_orbitals_data[orb.lower()] = [] 
            
            # Consume opening brace {
            i += 1
            # Now read data lines until closing brace }
            while i < len(lines):
                data_line = lines[i].strip()
                if data_line.startswith('}'):
                    break
                
                parts = data_line.split()
                if len(parts) >= 2: # At least n and zeta
                    try:
                        n_star = int(parts[0])
                        zeta = float(parts[1])
                        coeffs = [float(x) for x in parts[2:]]
                        
                        # Assign coeffs to respective orbitals
                        for idx, coeff in enumerate(coeffs):
                            if idx < len(current_orbitals):
                                orb_name = current_orbitals[idx].lower()
                                if abs(coeff) > 1e-8: # Filter negligible comparisons
                                    temp_orbitals_data[orb_name].append({
                                        'nStar': n_star,
                                        'zeta': zeta,
                                        'coeff': coeff
                                    })
                    except ValueError:
                        pass # Header lines or braces
                i += 1
            continue

        i += 1

    # Save last atom
    if current_atom and current_atom in atomic_numbers:
         basis_data[current_atom] = {
            'Z': atomic_numbers[current_atom],
            'name': names.get(current_atom, current_atom),
            'groundState': current_config,
            'orbitals': temp_orbitals_data
        }

    return basis_data

def generate_js_file(data, output_path):
    js_content = "// Clementi-Roetti 1974 Slater-Type Orbital Basis Sets\n"
    js_content += "// Generated from Clementi-Roetti.txt\n"
    js_content += "// 兼容 Worker 环境和主线程环境\n"
    js_content += "const globalScope = typeof self !== 'undefined' ? self : window;\n\n"
    js_content += "globalScope.SlaterBasis = {\n"
    js_content += "  version: '1.2',\n"
    js_content += "  reference: 'Clementi-Roetti (1974)',\n"
    
    # Sort elements by Z
    sorted_elements = sorted(data.items(), key=lambda x: x[1]['Z'])
    
    for symbol, info in sorted_elements:
        js_content += f"  '{symbol}': {{\n"
        js_content += f"    Z: {info['Z']},\n"
        js_content += f"    name: '{info['name']}',\n"
        js_content += f"    groundState: '{info['groundState']}',\n"
        js_content += "    orbitals: {\n"
        
        # Sort orbitals (1s, 2s, 2p...)
        # Custom sort key
        def orbital_sort_key(orb_name):
            n = int(orb_name[0])
            l_char = orb_name[1].lower()
            l = {'s':0, 'p':1, 'd':2, 'f':3}.get(l_char, 0)
            return n * 10 + l
            
        sorted_orbs = sorted(info['orbitals'].items(), key=lambda x: orbital_sort_key(x[0]))
        
        for orb_name, primitives in sorted_orbs:
            js_content += f"      '{orb_name}': [\n"
            for p in primitives:
               js_content += f"        {{ nStar: {p['nStar']}, zeta: {p['zeta']}, coeff: {p['coeff']} }},\n"
            js_content += "      ],\n"
        
        js_content += "    }\n"
        js_content += "  },\n"

    js_content += "};\n\n"
    js_content += "// Helper functions\n"
    js_content += "globalScope.getOrbitalData = function(symbol, orbitalName) {\n"
    js_content += "    const atom = globalScope.SlaterBasis[symbol];\n"
    js_content += "    if (!atom || !atom.orbitals[orbitalName]) return null;\n"
    js_content += "    return atom.orbitals[orbitalName];\n"
    js_content += "};\n"

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(js_content)

if __name__ == "__main__":
    filepath = "Clementi-Roetti.txt"
    output_path = "slater_basis.js"
    try:
        data = parse_clementi_file(filepath)
        print(f"Parsed {len(data)} atoms.")
        generate_js_file(data, output_path)
        print(f"Successfully generated {output_path}")
    except Exception as e:
        print(f"Error parsing file: {e}")
