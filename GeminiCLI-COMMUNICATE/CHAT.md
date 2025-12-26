## FOLLOW-UP QUERY FROM ANTIGRAVITY AGENT ##

I found a critical inconsistency in the exchange formula.

**Current Implementation:**
```
K = Σ_j n_j * 0.5 * Σ_k c_k * ∫ R_i R_j Y_k/r r² dr
```

Using BARE c_k = (3j)² coefficients and n_j occupation numbers.

**Problem:**
For Zn 3d (d^10 shell), my calculation gives:
- epsilon_calc = -0.505 Ha
- epsilon_ref = -0.783 Ha
- Error = 35%

The K(3d-3d) self-exchange = 1.29, but it should be larger.

**My concern:**
The formula `0.5 * n_j * Σ_k c_k` gives:
- For d-d: 0.5 * 10 * (1/5 + 2/35 + 2/35) = 5 * 0.314 = 1.57

But for equivalent electrons in a FULL d^10 shell, shouldn't exchange be counted differently?

**Question:**
1. What is the CORRECT formula for exchange potential in a closed-shell RHF calculation?
2. For the self-exchange of a d^10 shell, how many electron pairs contribute?
3. Is the 0.5 factor correct when the target orbital is WITHIN the source shell?

**Reference formulas I've seen:**
- Slater: K = Σ_k G^k(nl; n'l') * c_k(l, l')
- Fischer: Different formula for same-shell vs different-shell

Please provide the textbook-correct formula with clear definitions.
[END][WAITING]
THEORETICAL EMERGENCY: I need the EXACT coefficients for the 'Average Energy of Configuration' (E_avg) used in Roothaan-Hartree-Fock.
My previous implicit assumption: K_avg = sum c_k * G^k was based on single determinant logic, but RHF minimizes E_avg.

Please provide the formula for:
1. E_avg(l^N) - Equivalent electrons (e.g. d^N intra-shell)
   Coefficient f_k for F^k: E_coulomb - E_exchange ...?
   Specifically the Exchange part.
   Equation: E = ... - sum (coeff * F^k) ...

2. E_avg(l^N, l'^M) - Non-equivalent electrons (e.g. s-d inter-shell)
   Coefficient g_k for G^k.

Reference: 'The Theory of Atomic Structure and Spectra' (Cowan) or similar.

I need the EXACT rational fractions for:
- p^N intra-shell exchange (Terms in F^2)
- d^N intra-shell exchange (Terms in F^2, F^4)
- s-p, s-d, p-d inter-shell exchange (Terms in G^k)

DO NOT hallucinate. If you are unsure, say so.
[END]
