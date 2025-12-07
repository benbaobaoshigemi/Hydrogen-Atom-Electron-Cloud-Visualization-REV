import urllib.request
import json
import sys

def get_basis_list():
    url = "https://www.basissetexchange.org/api/basis/"
    try:
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            return data
    except Exception as e:
        print(f"Error fetching basis list: {e}")
        return None

def get_basis_data(basis_name, elements):
    # elements is a list of integers (Z)
    # elements_str = ",".join(map(str, elements))
    # API format: /api/basis/<basis_name>/format/json?elements=<elements>
    # Try getting all available elements if not specified? 
    # Usually we can get the whole basis.
    url = f"https://www.basissetexchange.org/api/basis/{basis_name}/format/json"
    try:
        print(f"Fetching {url}...")
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            return data
    except Exception as e:
        print(f"Error fetching basis data for {basis_name}: {e}")
        # Try raw format or other endpoints if needed
        return None

def main():
    print("Fetching basis set list...")
    basis_list = get_basis_list()
    
    if not basis_list:
        print("Failed to get basis list.")
        return

    # Case insensitive search
    candidates = []
    for key in basis_list:
        if "clementi" in key.lower() or "roetti" in key.lower() or "sto" in key.lower() or "bunge" in key.lower():
            candidates.append(key)
    
    print(f"Found {len(candidates)} candidate basis sets:")
    for c in candidates:
        print(f" - {c}")

    # If we find a good candidate, try to fetch it.
    # Prioritize 'clementi-roetti' matches
    target_basis = None
    for c in candidates:
        if "clementi" in c.lower():
            target_basis = c
            break
    
    if not target_basis and candidates:
        # Fallback to STO-6G or similar if specific CR not found? No, user wants CR.
        pass

    if target_basis:
        print(f"\nAttempting to download data for: {target_basis}")
        data = get_basis_data(target_basis, range(1, 55))
        if data:
            with open("clementi_data.json", "w") as f:
                json.dump(data, f, indent=2)
            print("Successfully saved clementi_data.json")
        else:
            print("Failed to download data.")
    else:
        print("No specific Clementi-Roetti basis found in the list.")

if __name__ == "__main__":
    main()
