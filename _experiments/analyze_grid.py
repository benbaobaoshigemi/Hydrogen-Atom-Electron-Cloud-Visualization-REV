# -*- coding: utf-8 -*-
"""Analyze grid requirements for heavy atoms"""
import numpy as np
from pathlib import Path
import json, re

def parse_slater_basis(filepath):
    content = filepath.read_text(encoding='utf-8')
    start = content.find('globalScope.SlaterBasis = ')
    json_start = start + len('globalScope.SlaterBasis = ')
    depth, json_end, in_str, escape = 0, json_start, False, False
    for i, c in enumerate(content[json_start:]):
        if escape: escape = False; continue
        if c == '\\': escape = True; continue
        if c == '"': in_str = not in_str; continue
        if in_str: continue
        if c == '{': depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0: json_end = json_start + i + 1; break
    json_str = content[json_start:json_end]
    json_str = re.sub(r',\s*}', '}', json_str)
    json_str = re.sub(r',\s*]', ']', json_str)
    return json.loads(json_str)

data = parse_slater_basis(Path('../slater_basis.js'))

print("Zeta ranges for closed-shell atoms:")
print("="*60)
for atom in ['He', 'Ne', 'Ar', 'Kr']:
    atom_data = data[atom]
    all_zetas = []
    for orb_name, terms in atom_data['orbitals'].items():
        zetas = [t['zeta'] for t in terms]
        all_zetas.extend(zetas)
        print(f'{atom} {orb_name}: zeta_min={min(zetas):.2f}, zeta_max={max(zetas):.2f}')
    print(f'  -> needs r_min < {1/max(all_zetas):.1e}, r_max > {5/min(all_zetas):.1f}')
    print()

print("Current grid: r_min=1e-6, r_max=80, N=3000")
print("This may be insufficient for Kr with zeta > 50")
