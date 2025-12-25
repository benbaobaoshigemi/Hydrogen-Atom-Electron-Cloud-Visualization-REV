const PI = Math.PI;

function realYlm_value(l, m, type, x, y, z) {
    const norm_p = Math.sqrt(3 / (4 * PI));
    if (l === 0) return 1 / Math.sqrt(4 * PI);
    if (l === 1) {
        if (m === 0) return norm_p * z;
        if (type === 'c') return norm_p * x;
        if (type === 's') return norm_p * y;
    }
    return 0;
}

function jacobiSVD(A) {
    const m = A.length;
    const n = A[0].length;
    let U = A.map(row => [...row]);
    let V = [];
    for (let i = 0; i < n; i++) {
        const row = new Array(n).fill(0);
        row[i] = 1;
        V.push(row);
    }
    const EPSILON = 1e-15;
    for (let iter = 0; iter < 50; iter++) {
        for (let i = 0; i < n - 1; i++) {
            for (let j = i + 1; j < n; j++) {
                let alpha = 0, beta = 0, gamma = 0;
                for (let k = 0; k < m; k++) {
                    alpha += U[k][i] * U[k][i];
                    beta += U[k][j] * U[k][j];
                    gamma += U[k][i] * U[k][j];
                }
                if (Math.abs(gamma) < EPSILON) continue;

                let zeta = (beta - alpha) / (2 * gamma);
                let t;
                if (zeta === 0) {
                    t = 1;
                } else {
                    t = Math.sign(zeta) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
                }

                let c = 1 / Math.sqrt(1 + t * t);
                let s = c * t;
                for (let k = 0; k < m; k++) {
                    let t1 = U[k][i]; let t2 = U[k][j];
                    U[k][i] = c * t1 - s * t2; U[k][j] = s * t1 + c * t2;
                }
                for (let k = 0; k < n; k++) {
                    let t1 = V[k][i]; let t2 = V[k][j];
                    V[k][i] = c * t1 - s * t2; V[k][j] = s * t1 + c * t2;
                }
            }
        }
    }
    const S = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
        let sum = 0;
        for (let k = 0; k < m; k++) sum += U[k][i] * U[k][i];
        S[i] = Math.sqrt(sum);
        if (S[i] > EPSILON) {
            for (let k = 0; k < m; k++) U[k][i] /= S[i];
        }
    }
    return { U, S, V };
}

function matTranspose(A) {
    return A[0].map((_, colIndex) => A.map(row => row[colIndex]));
}

function matMul(A, B) {
    const m = A.length, n = A[0].length, p = B[0].length;
    const C = Array.from({ length: m }, () => new Array(p).fill(0));
    for (let i = 0; i < m; i++)
        for (let j = 0; j < p; j++)
            for (let k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

function buildDirectionMatrix(directions, orbitalParams) {
    return directions.map(d => {
        const x = d[0], y = d[1], z = d[2];
        return orbitalParams.map(p => realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, x, y, z));
    });
}

// Test case: 4px + 4py (pure p-orbital hybridization)
// With parity-aware direction generation, we should get orthogonal directions
const directions = [[1, 0, 0], [0, 1, 0]]; // Parity-aware output
const orbitalParams = [
    { angKey: { l: 1, m: 1, t: 'c' } }, // px
    { angKey: { l: 1, m: 1, t: 's' } }  // py
];

console.log("=== Testing Parity-Aware Direction Generation ===");
console.log("Directions:", directions);

const A = buildDirectionMatrix(directions, orbitalParams);
console.log("Matrix A:");
console.log(A);

const { U, S, V } = jacobiSVD(A);
console.log("Singular Values S:", S);

const VT = matTranspose(V);
const Q = matMul(U, VT);

console.log("\nCoefficients Q:");
console.log(Q);

// Verify orthogonality
const row1 = Q[0], row2 = Q[1];
const dot = row1[0] * row2[0] + row1[1] * row2[1];
const norm1 = Math.sqrt(row1[0] ** 2 + row1[1] ** 2);
const norm2 = Math.sqrt(row2[0] ** 2 + row2[1] ** 2);
console.log("\nRow 1 Norm (should be 1):", norm1);
console.log("Row 2 Norm (should be 1):", norm2);
console.log("Dot product (should be ~0):", dot);

// Check if the two hybrid orbitals are distinct
console.log("\n=== Verification Summary ===");
if (Math.abs(dot) < 0.01 && Math.abs(norm1 - 1) < 0.01 && Math.abs(norm2 - 1) < 0.01) {
    console.log("✓ SUCCESS: Hybrid orbitals are orthogonal and normalized!");
    console.log("  Hybrid 1 = " + Q[0][0].toFixed(3) + "*px + " + Q[0][1].toFixed(3) + "*py");
    console.log("  Hybrid 2 = " + Q[1][0].toFixed(3) + "*px + " + Q[1][1].toFixed(3) + "*py");
} else {
    console.log("✗ FAILURE: Hybrid orbitals are NOT properly orthogonal!");
}
