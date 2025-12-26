
import math

def factorial(n):
    return math.factorial(n)

def wigner_3j(j1, j2, j3, m1, m2, m3):
    if (m1 + m2 + m3 != 0): return 0.0
    if not (abs(j1 - j2) <= j3 <= j1 + j2): return 0.0
    
    t_min = max(0, -(j3 - j2 + m1), -(j3 - j1 - m2))
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    
    sum_t = 0.0
    for t in range(t_min, t_max + 1):
        num = (-1)**t
        den = factorial(t) * factorial(j1 + j2 - j3 - t) * \
              factorial(j1 - m1 - t) * factorial(j2 + m2 - t) * \
              factorial(j3 - j2 + m1 + t) * factorial(j3 - j1 - m2 + t)
        sum_t += num / den
        
    triangle = factorial(j1 + j2 - j3) * factorial(j1 + j3 - j2) * factorial(j2 + j3 - j1) / factorial(j1 + j2 + j3 + 1)
    prefactor = math.sqrt(triangle) * math.sqrt(factorial(j1 + m1) * factorial(j1 - m1) * factorial(j2 + m2) * factorial(j2 - m2) * factorial(j3 + m3) * factorial(j3 - m3))
    
    return prefactor * sum_t * (-1)**(j1 - j2 - m3)

# Calculate the coefficient for exchange energy of a closed shell n'l' acting on orbital nl
# Formula for the exchange potential term acting on P_nl(r):
# V_X P_nl = - sum_{n'l' (closed)} (Occupation_{n'l'} / (4l+2) ?) ...
#
# Actually, let's use the standard result from Condon & Shortley or Cowan.
# For closed shells:
# E_exchange = - sum_{pairs} G^k(nl, n'l') * coeff
#
# But the user is asking about the Radial Fock Equation potential V_X(r).
# The term is: - sum_{n'l'} sum_k  C_{ll'k} Y^k(nl, n'l'; r) P_{n'l'}(r) / r
#
# Where C_{ll'k} is the coefficient we need.
# For closed shells, the standard result (e.g., from Froese Fischer's MCHF book, eq 1.15) is:
# C_{ll'k} = (2l'+1) * (3j(l, k, l', 0, 0, 0))^2
#
# Wait, let's verify if there is a 1/2 factor.
# The user's formula has 0.5 * n_j.
# For a closed shell, n_j = 2(2l'+1).
# So 0.5 * n_j = 2l'+1.
# If the user uses 0.5 * n_j as the prefactor, they effectively have (2l'+1).
#
# Let's compute the bare 3j^2 and the full (2l'+1)*3j^2 values again to be sure.

print("Checking coefficients...")
pairs = [("s-s", 0, 0), ("s-p", 0, 1), ("p-s", 1, 0), ("p-p", 1, 1), ("p-d", 1, 2), ("d-p", 2, 1), ("d-d", 2, 2)]

for label, l, lp in pairs:
    print(f"\nInteraction {label} (Target l={l}, Source l'={lp}):")
    # Possible k values: |l-l'| to l+l' in steps of 2 (parity conservation)
    ks = range(abs(l-lp), l+lp+1, 2)
    
    for k in ks:
        three_j = wigner_3j(l, k, lp, 0, 0, 0)
        bare = three_j**2
        
        # Standard coefficient for closed shell exchange in Fock equation
        # C_{ll'k} = (2l'+1) * 3j^2
        full = (2*lp + 1) * bare
        
        # User's potential formula: 0.5 * n_j * c_k_user
        # If user uses c_k_user = bare, and n_j = 2(2lp+1), then
        # Term = 0.5 * 2(2lp+1) * bare = (2lp+1) * bare.
        # This matches the standard result!
        
        print(f"  k={k}: Bare 3j^2 = {bare:.5f} ({float(bare).as_integer_ratio()})")
        print(f"       Full (2l'+1)*3j^2 = {full:.5f} ({float(full).as_integer_ratio()})")
