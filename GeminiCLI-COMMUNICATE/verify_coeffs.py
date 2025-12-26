
import math

def factorial(n):
    return math.factorial(n)

def wigner_3j(j1, j2, j3, m1, m2, m3):
    if (m1 + m2 + m3 != 0): return 0
    if not (abs(j1 - j2) <= j3 <= j1 + j2): return 0
    
    # Precompute factorials
    # ... straightforward implementation of Racah formula
    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2
    
    t_min = max(0, -(j3 - j2 + m1), -(j3 - j1 - m2))
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    
    sum_t = 0
    for t in range(t_min, t_max + 1):
        num = (-1)**t
        den = factorial(t) * factorial(j1 + j2 - j3 - t) * \
              factorial(j1 - m1 - t) * factorial(j2 + m2 - t) * \
              factorial(j3 - j2 + m1 + t) * factorial(j3 - j1 - m2 + t)
        sum_t += num / den
        
    triangle = factorial(j1 + j2 - j3) * factorial(j1 + j3 - j2) * factorial(j2 + j3 - j1) / factorial(j1 + j2 + j3 + 1)
    prefactor = math.sqrt(triangle) * math.sqrt(factorial(j1 + m1) * factorial(j1 - m1) * factorial(j2 + m2) * factorial(j2 - m2) * factorial(j3 + m3) * factorial(j3 - m3))
    
    return prefactor * sum_t * (-1)**(j1 - j2 - m3)

def calc_exchange_coeff(l_target, l_source, k):
    # Sum over m_source of the angular part squared
    # Angular part for a pair lm, l'm' is (2l+1)(2l'+1) [ (3j_m) * (3j_0) ]^2 / (4pi/(2k+1))? 
    # Let's stick to the analytical result derived: Coeff = (2l_source + 1) * (3j_000)^2
    
    # Calculate 3j_000
    tj = wigner_3j(l_target, k, l_source, 0, 0, 0)
    bare = tj**2
    
    full_coeff = (2*l_source + 1) * bare
    return bare, full_coeff

print("Verification of Coefficients:")
print("-" * 40)
pairs = [
    ("s-s", 0, 0, [0]),
    ("s-p", 0, 1, [1]),
    ("p-s", 1, 0, [1]),
    ("p-p", 1, 1, [0, 2]),
    ("p-d", 1, 2, [1, 3]),
    ("d-p", 2, 1, [1, 3]),
    ("d-d", 2, 2, [0, 2, 4])
]

for label, l1, l2, ks in pairs:
    print(f"Interaction {label} (Target l={l1}, Source l={l2}):")
    for k in ks:
        bare, full = calc_exchange_coeff(l1, l2, k)
        print(f"  k={k}: Bare 3j^2 = {bare:.5f} ({bare.as_integer_ratio()}), Full (with 2l'+1) = {full:.5f} ({full.as_integer_ratio()})")
