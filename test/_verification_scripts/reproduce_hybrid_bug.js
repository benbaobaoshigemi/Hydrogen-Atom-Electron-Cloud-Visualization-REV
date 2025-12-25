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
        let maxError = 0;
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
                let t = Math.sign(zeta) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
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

// Reproduction
const n = 2;
// optimizeThomson(2) with golden ratio gives [0, 1, 0] and [0, -1, 0]
const directions = [[0, 1, 0], [0, -1, 0]];
const orbitalParams = [
    { angKey: { l: 1, m: 1, t: 'c' } }, // px
    { angKey: { l: 1, m: 1, t: 's' } }  // py
];

const norm_p = Math.sqrt(3 / (4 * PI));
const A = directions.map(d => {
    const x = d[0], y = d[1], z = d[2];
    return [
        norm_p * x, // px
        norm_p * y  // py
    ];
});

console.log("Matrix A:");
console.log(A);

const { U, S, V } = jacobiSVD(A);
console.log("Singular Values S:", S);
const VT = matTranspose(V);
const Q = matMul(U, VT);

console.log("Coefficients Q:");
console.log(Q);
