import json

def validate_db():
    with open('koga_basis_database.json', 'r', encoding='utf-8') as f:
        db = json.load(f)
    
    print(f"Total elements: {len(db)}")
    
    stats = {
        "missing_orbitals": [],
        "empty_basis": [],
        "orbital_counts": {}
    }
    
    for symbol, data in db.items():
        if symbol == "Au_R": continue
        
        orbitals = data.get("orbitals", {})
        if not orbitals:
            stats["missing_orbitals"].append(symbol)
            continue
            
        stats["orbital_counts"][symbol] = len(orbitals)
        
        for orb_name, orb_data in orbitals.items():
            basis = orb_data.get("basis", [])
            if not basis:
                stats["empty_basis"].append(f"{symbol}-{orb_name}")
                
    # Check Ti specifically
    if "Ti" in db:
        ti = db["Ti"]
        print(f"Ti orbitals: {list(ti['orbitals'].keys())}")
        print(f"Ti 4s basis count: {len(ti['orbitals']['4s']['basis'])}")
        
    # Check Xe specifically
    if "Xe" in db:
        xe = db["Xe"]
        print(f"Xe orbitals: {list(xe['orbitals'].keys())}")

    print(f"Missing orbitals count: {len(stats['missing_orbitals'])}")
    if stats["missing_orbitals"]:
        print(f"Missing: {stats['missing_orbitals']}")
        
    print(f"Empty basis count: {len(stats['empty_basis'])}")
    if stats["empty_basis"]:
        print(f"Empty: {stats['empty_basis'][:10]}...")

if __name__ == "__main__":
    validate_db()
