# -*- coding: utf-8 -*-
"""
杂化轨道几何优化规模测试 v2
使用 Thomson 问题（球面电荷排斥）生成最优方向
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from scipy.optimize import minimize
from itertools import combinations
import time
import matplotlib.pyplot as plt
from collections import defaultdict

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# 所有可用轨道定义 (name, l, m)
ALL_ORBITALS = [
    ('2s', 0, 0), ('2px', 1, 1), ('2py', 1, -1), ('2pz', 1, 0),
    ('3s', 0, 0), ('3px', 1, 1), ('3py', 1, -1), ('3pz', 1, 0),
    ('3dz2', 2, 0), ('3dxz', 2, 1), ('3dyz', 2, -1), ('3dx2y2', 2, 2), ('3dxy', 2, -2),
    ('4s', 0, 0), ('4px', 1, 1), ('4py', 1, -1), ('4pz', 1, 0),
]

# 缓存 Thomson 优化结果
_thomson_cache = {}

def thomson_energy(flat_points, n):
    """计算 N 个球面点的静电势能"""
    points = flat_points.reshape(n, 3)
    # 投影到球面
    norms = np.linalg.norm(points, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    points = points / norms
    
    energy = 0.0
    for i in range(n):
        for j in range(i+1, n):
            dist = np.linalg.norm(points[i] - points[j])
            energy += 1.0 / (dist + 1e-10)
    return energy

def optimal_directions_thomson(n):
    """用 Thomson 问题找到 N 个最优分布的方向"""
    if n in _thomson_cache:
        return _thomson_cache[n].copy()
    
    if n == 1:
        return np.array([[0, 0, 1]])
    
    # 多次随机初始化，取最优
    best_energy = float('inf')
    best_dirs = None
    
    for _ in range(3):  # 3次尝试
        x0 = np.random.randn(n * 3)
        result = minimize(thomson_energy, x0, args=(n,), method='L-BFGS-B', 
                         options={'maxiter': 100})
        
        if result.fun < best_energy:
            best_energy = result.fun
            directions = result.x.reshape(n, 3)
            norms = np.linalg.norm(directions, axis=1, keepdims=True)
            best_dirs = directions / np.maximum(norms, 1e-10)
    
    _thomson_cache[n] = best_dirs
    return best_dirs.copy()

def angular_function(l, m, theta, phi):
    """实球谐函数在给定方向的值"""
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    
    if l == 0:  # s
        return 1.0 / np.sqrt(4 * np.pi)
    elif l == 1:
        norm = np.sqrt(3 / (4 * np.pi))
        if m == 1:   return norm * x
        if m == -1:  return norm * y
        if m == 0:   return norm * z
    elif l == 2:
        if m == 0:   return np.sqrt(5 / (16 * np.pi)) * (3 * z**2 - 1)
        if m == 1:   return np.sqrt(15 / (4 * np.pi)) * x * z
        if m == -1:  return np.sqrt(15 / (4 * np.pi)) * y * z
        if m == 2:   return np.sqrt(15 / (16 * np.pi)) * (x**2 - y**2)
        if m == -2:  return np.sqrt(15 / (4 * np.pi)) * x * y
    return 0

def direction_to_angles(d):
    d = d / np.linalg.norm(d)
    theta = np.arccos(np.clip(d[2], -1, 1))
    phi = np.arctan2(d[1], d[0])
    return theta, phi

def build_direction_matrix(directions, orbitals):
    n_dirs = len(directions)
    n_orbitals = len(orbitals)
    A = np.zeros((n_dirs, n_orbitals))
    
    for i, d in enumerate(directions):
        theta, phi = direction_to_angles(d)
        for j, (name, l, m) in enumerate(orbitals):
            A[i, j] = angular_function(l, m, theta, phi)
    return A

def optimize_hybrid_coefficients(directions, orbitals):
    A = build_direction_matrix(directions, orbitals)
    U, S, Vt = np.linalg.svd(A, full_matrices=True)
    n = min(len(directions), len(orbitals))
    C = Vt[:n, :].T @ U[:, :n].T
    return C.T

def verify_orthonormality(C, tol=1e-6):
    if C.shape[0] != C.shape[1]:
        return False
    result = C @ C.T
    identity = np.eye(C.shape[0])
    return np.max(np.abs(result - identity)) < tol

def analyze_s_content(C, orbitals):
    s_indices = [i for i, (name, l, m) in enumerate(orbitals) if l == 0]
    if not s_indices:
        return np.zeros(C.shape[0])
    return sum(C[:, i] ** 2 for i in s_indices)

def run_test():
    results = defaultdict(list)
    timing_stats = defaultdict(list)
    standard_hybrids = {}
    
    print("Thomson 方向优化 + 杂化系数计算")
    print("=" * 70)
    
    # 预热 Thomson 缓存
    print("\n预热 Thomson 缓存 (n=2..6)...")
    for n in range(2, 7):
        start = time.perf_counter()
        optimal_directions_thomson(n)
        elapsed = (time.perf_counter() - start) * 1000
        print(f"  n={n}: {elapsed:.1f} ms")
    
    print("\n" + "=" * 70)
    
    # 只测试部分组合（减少时间）
    for n in range(2, 7):
        combos = list(combinations(range(len(ALL_ORBITALS)), n))
        # 限制测试数量
        max_test = min(len(combos), 500)
        test_combos = combos[:max_test]
        
        print(f"\n{n}轨道组合: 测试 {max_test}/{len(combos)} 种")
        
        directions = optimal_directions_thomson(n)  # 使用缓存的方向
        
        for combo in test_combos:
            orbitals = [ALL_ORBITALS[i] for i in combo]
            orbital_names = [o[0] for o in orbitals]
            
            start = time.perf_counter()
            C = optimize_hybrid_coefficients(directions, orbitals)
            elapsed = (time.perf_counter() - start) * 1000
            
            timing_stats[n].append(elapsed)
            is_valid = verify_orthonormality(C)
            
            results[n].append({
                'orbitals': orbital_names,
                'coefficients': C,
                'valid': is_valid,
                'time_ms': elapsed
            })
            
            # 检测标准杂化
            orbital_set = set(orbital_names)
            
            if orbital_set == {'2s', '2pz'}:
                standard_hybrids['sp'] = (orbital_names, C, elapsed)
            if orbital_set == {'2s', '2px', '2py'}:
                standard_hybrids['sp2'] = (orbital_names, C, elapsed)
            if orbital_set == {'2s', '2px', '2py', '2pz'}:
                standard_hybrids['sp3'] = (orbital_names, C, elapsed)
            if orbital_set == {'2s', '2px', '2py', '2pz', '3dz2'}:
                standard_hybrids['sp3d'] = (orbital_names, C, elapsed)
            if orbital_set == {'2s', '2px', '2py', '2pz', '3dz2', '3dx2y2'}:
                standard_hybrids['sp3d2'] = (orbital_names, C, elapsed)
        
        times = timing_stats[n]
        valid_count = sum(1 for r in results[n] if r['valid'])
        print(f"  平均时间: {np.mean(times):.3f} ms")
        print(f"  最大时间: {np.max(times):.3f} ms")
        print(f"  有效率: {valid_count / len(results[n]) * 100:.1f}%")
    
    return results, timing_stats, standard_hybrids

def print_standard_hybrids(standard_hybrids):
    print("\n" + "=" * 70)
    print("标准杂化轨道系数 (Thomson 方向优化)")
    print("=" * 70)
    
    for name in ['sp', 'sp2', 'sp3', 'sp3d', 'sp3d2']:
        if name not in standard_hybrids:
            print(f"\n{name}: 未找到")
            continue
        
        orbitals, C, elapsed = standard_hybrids[name]
        display_name = name.replace('2', '\u00b2').replace('3', '\u00b3')
        print(f"\n[{display_name}] ({elapsed:.3f} ms)")
        print(f"轨道: {', '.join(orbitals)}")
        print("-" * 50)
        
        # 系数表
        print(f"{'hybrid':<10}", end="")
        for o in orbitals:
            print(f"{o:>8}", end="")
        print()
        
        for i, row in enumerate(C):
            print(f"phi_{i+1:<6}", end="")
            for val in row:
                print(f"{val:>8.4f}", end="")
            print()
        
        # s成分
        s_idx = [j for j, o in enumerate(orbitals) if o.endswith('s') and 'p' not in o and 'd' not in o]
        if s_idx:
            print(f"\ns 成分:")
            for i, row in enumerate(C):
                s_content = sum(row[j]**2 for j in s_idx)
                bar = "#" * int(s_content * 40)
                print(f"  phi_{i+1}: {s_content:.4f} ({s_content*100:.1f}%) {bar}")
            total_s = sum(sum(C[i, j]**2 for j in s_idx) for i in range(C.shape[0]))
            print(f"  总和: {total_s:.4f}")

if __name__ == "__main__":
    total_start = time.time()
    
    results, timing_stats, standard_hybrids = run_test()
    
    print(f"\n总时间: {time.time() - total_start:.2f} 秒")
    
    print_standard_hybrids(standard_hybrids)
    
    print("\n" + "=" * 70)
    print("完成!")
