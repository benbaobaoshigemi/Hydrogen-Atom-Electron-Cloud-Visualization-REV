
const fs = require('fs');
const path = require('path');

// Mock global environment
const globalScope = {};
global.window = globalScope;
global.self = globalScope;

// Load SlaterBasis
const slaterBasisPath = path.join(__dirname, '..', 'slater_basis.js');
const slaterBasisContent = fs.readFileSync(slaterBasisPath, 'utf8');
eval(slaterBasisContent);

const SlaterBasis = globalScope.SlaterBasis;

function factorial(n) {
    if (n <= 1) return 1;
    let res = 1;
    for (let i = 2; i <= n; i++) res *= i;
    return res;
}

function calculateOverlap(ti, tj) {
    const ni = ti.nStar;
    const zi = ti.zeta;
    const nj = tj.nStar;
    const zj = tj.zeta;

    // Normalization constant for STO: (2zeta)^(n+0.5) / sqrt( (2n)! )
    const Ni = Math.pow(2 * zi, ni + 0.5) / Math.sqrt(factorial(2 * ni));
    const Nj = Math.pow(2 * zj, nj + 0.5) / Math.sqrt(factorial(2 * nj));

    // Integral r^2 * R_i(r) * R_j(r) dr
    // R(r) ~ r^(n-1) * exp(-zeta*r)
    // Integrand ~ r^(ni-1 + nj-1 + 2) * exp(-(zi+zj)r) = r^(ni+nj) * exp(-(zi+zj)r)
    // Gamma integral: integral_0^inf x^k e^(-a x) dx = k! / a^(k+1)

    const k = ni + nj;
    const alpha = zi + zj;
    const integral = factorial(k) / Math.pow(alpha, k + 1);

    return Ni * Nj * integral;
}

function checkNormalization(elementSymbol, orbitalName) {
    const element = SlaterBasis[elementSymbol];
    if (!element) {
        console.error(`Element ${elementSymbol} not found.`);
        return false;
    }
    const orbitals = element.orbitals;
    if (!orbitals || !orbitals[orbitalName]) {
        console.error(`Orbital ${orbitalName} for ${elementSymbol} not found.`);
        return false;
    }

    const basis = orbitals[orbitalName];
    let norm = 0;

    for (let i = 0; i < basis.length; i++) {
        for (let j = 0; j < basis.length; j++) {
            const overlap = calculateOverlap(basis[i], basis[j]);
            norm += basis[i].coeff * basis[j].coeff * overlap;
        }
    }

    const diff = Math.abs(norm - 1.0);
    const passed = diff < 1e-4;
    console.log(`[${elementSymbol} ${orbitalName}] Norm: ${norm.toFixed(6)} | Diff: ${diff.toExponential(2)} | ${passed ? 'PASS' : 'FAIL'}`);
    return passed;
}

console.log("=== Verifying Koga Basis Normalization (Recent Additions) ===");

const elementsToCheck = ['V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'];
const orbitalTypes = ['1s', '2p', '3d', '4s', '4p'];

let allPassed = true;

for (const el of elementsToCheck) {
    if (!SlaterBasis[el]) {
        console.warn(`Skipping ${el} (Not found in basis file)`);
        continue;
    }
    console.log(`\n--- Element: ${el} ---`);
    for (const orb of Object.keys(SlaterBasis[el].orbitals)) {
        if (!checkNormalization(el, orb)) {
            allPassed = false;
        }
    }
}

if (allPassed) {
    console.log("\nAll checked orbitals are correctly normalized.");
} else {
    console.log("\nSome normalization checks failed.");
}
