import json
import os

def build_engine():
    json_path = 'koga_basis_database.json'
    target_js = 'slater_basis.js'
    
    with open(json_path, 'r', encoding='utf-8') as f:
        db = json.load(f)
    
    # Flatten the structure
    flattened_db = {}
    for symbol, data in db.items():
        # Keep all available atoms (Z=1 to 103)

        new_atom = {
            "Z": data["Z"],
            "name": data["name"],
            "orbitals": {},
            "E_tot": data.get("E_tot", 0) # Inject Total Energy
        }
        
        # Original groundStates from Koga are preferred
        if "groundState" in data:
            new_atom["groundState"] = data["groundState"]
        elif symbol == "H":
            new_atom["groundState"] = "1s¹"
        
        for orb_key, orb_data in data["orbitals"].items():
            # Inject Orbital Energy as a property of the orbital array
            # We attach it to the array object in JS, or better, change structure?
            # User wants minimal changes. Let's make "orbitals" map to an object if possible?
            # Wait, existing code expects an array of basis functions.
            # To be safe and compatible, we can't change the structure of `new_atom["orbitals"][orb_key]` from array to object without breaking density code.
            # BUT we can add a parallel structure or use a property on the array (less clean in JSON).
            # Solution: We will keep the array as the value of the key, but maybe add a separate "energies" map?
            # No, 'slater_basis.js' is the only source.
            # Let's check how `radialR` accesses it. `window.SlaterBasis[atomType].orbitals[orbitalKey]`.
            # If I make `orbitals[orbitalKey]` an object { basis: [], energy: -0.5 }, I break `basis.length` checks.
            # Let's look at `physics.js`: `const basis = atomData.orbitals[orbitalKey];` ... `for (let i = 0; i < basis.length; i++)`
            # YES, it breaks.
            # Alternative: Add `energies` map to the atom object.
            new_atom["orbitals"][orb_key] = orb_data["basis"]
        
        # Add energies map
        new_atom["energies"] = {}
        for orb_key, orb_data in data["orbitals"].items():
             new_atom["energies"][orb_key] = orb_data.get("energy", 0)
            
        flattened_db[symbol] = new_atom

    # Generate JS content
    js_header = """// Koga (1999)/Koga (2000) High-Precision Slater-Type Orbital Basis Sets
// Expanded for Z=1 to Z=103 + Relativistic Gold (Au_R)
// Compatible with the Electron Cloud Visualization physics engine

const globalScope = typeof self !== 'undefined' ? self : window;

globalScope.SlaterBasis = """

    with open(target_js, 'w', encoding='utf-8') as f:
        f.write(js_header)
        json.dump(flattened_db, f, indent=2)
        f.write(";\n")
        
    print(f"Successfully generated {target_js} with {len(flattened_db)} atoms.")

if __name__ == "__main__":
    build_engine()
