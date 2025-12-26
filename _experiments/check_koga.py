# -*- coding: utf-8 -*-
"""Check Koga basis parameters"""
import json, re
from pathlib import Path

content = Path('../slater_basis.js').read_text(encoding='utf-8')
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
data = json.loads(json_str)

# 检查 Zn 3d 的基组参数
print('Zn 3d orbital basis:')
for term in data['Zn']['orbitals']['3d']:
    print(f"  nStar={term['nStar']}, zeta={term['zeta']:.4f}, coeff={term['coeff']:.6f}")

print()
print('Zn 4s orbital basis:')
for term in data['Zn']['orbitals']['4s']:
    print(f"  nStar={term['nStar']}, zeta={term['zeta']:.4f}, coeff={term['coeff']:.6f}")

print()
print('Zn orbital energies (reference):')
for orb, e in data['Zn']['energies'].items():
    print(f'  {orb}: {e:.6f} Ha')
