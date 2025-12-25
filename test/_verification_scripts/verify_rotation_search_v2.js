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

                // Fixed Jacobi logic
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
    const norm_p = Math.sqrt(3 / (4 * PI));
    return directions.map(d => {
        const x = d[0], y = d[1], z = d[2];
        return orbitalParams.map(p => realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, x, y, z));
    });
}

const rotateVec = (v, axis, angle) => {
    const cosA = Math.cos(angle);
    const sinA = Math.sin(angle);
    const dot = v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2];
    const cross = [
        v[1] * axis[2] - v[2] * axis[1],
        v[2] * axis[0] - v[0] * axis[2],
        v[0] * axis[1] - v[1] * axis[0]
    ];
    return [
        v[0] * cosA + cross[0] * sinA + axis[0] * dot * (1 - cosA),
        v[1] * cosA + cross[1] * sinA + axis[1] * dot * (1 - cosA),
        v[2] * cosA + cross[2] * sinA + axis[2] * dot * (1 - cosA)
    ];
};

// Simulation
const baseDirections = [[0, 1, 0], [0, -1, 0]]; // Original bad points
const orbitalParams = [
    { angKey: { l: 1, m: 1, t: 'c' } }, // px
    { angKey: { l: 1, m: 1, t: 's' } }  // py
];

const searchAxes = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]];
const searchAngles = [0, PI / 4, PI / 3, PI / 6, 0.12345];

let bestDirections = baseDirections;
let maxMetric = -1;

for (const axis of searchAxes) {
    const norm = Math.sqrt(axis[0] ** 2 + axis[1] ** 2 + axis[2] ** 2);
    const nAxis = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
    for (const angle of searchAngles) {
        const rotatedDirs = baseDirections.map(d => rotateVec(d, nAxis, angle));
        const A = buildDirectionMatrix(rotatedDirs, orbitalParams);
        const { S } = jacobiSVD(A);
        let metric = 1.0;
        for (let s of S) metric *= s;
        if (metric > maxMetric) {
            maxMetric = metric;
            bestDirections = rotatedDirs;
        }
    }
}

console.log("Best Metric:", maxMetric);
const finalA = buildDirectionMatrix(bestDirections, orbitalParams);
const { U, S, V } = jacobiSVD(finalA);
const VT = matTranspose(V);
const Q = matMul(U, VT);

console.log("Coefficients Q:");
console.log(Q);

// Check orthogonality
const row1 = Q[0], row2 = Q[1];
const dot = row1[0] * row2[0] + row1[1] * row2[1];
const norm1 = Math.sqrt(row1[0] ** 2 + row1[1] ** 2);
const norm2 = Math.sqrt(row2[0] ** 2 + row2[1] ** 2);
console.log("Row 1 Norm:", norm1);
console.log("Row 2 Norm:", norm2);
console.log("Dot product (should be ~0):", dot);
