// 离屏测试杂化系数计算算法 - 与 physics.js 结构完全一致
// 运行: node test_hybrid_js.js

// ==================== 算法实现（从 physics.js 精确复制） ====================

const PI = Math.PI;
const THOMSON_CACHE = {};

function thomsonEnergy(points) {
    let energy = 0;
    const n = points.length;
    for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
            const dx = points[i][0] - points[j][0];
            const dy = points[i][1] - points[j][1];
            const dz = points[i][2] - points[j][2];
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (dist > 1e-10) {
                energy += 1 / dist;
            }
        }
    }
    return energy;
}

function optimizeThomson(n, maxIter = 200, lr = 0.1) {
    if (THOMSON_CACHE[n]) return THOMSON_CACHE[n];

    if (n === 1) {
        THOMSON_CACHE[n] = [[0, 0, 1]];
        return THOMSON_CACHE[n];
    }
    if (n === 2) {
        THOMSON_CACHE[n] = [[0, 0, 1], [0, 0, -1]];
        return THOMSON_CACHE[n];
    }

    let points = [];
    const goldenRatio = (1 + Math.sqrt(5)) / 2;
    for (let i = 0; i < n; i++) {
        const y = 1 - (i / (n - 1)) * 2;
        const radius = Math.sqrt(1 - y * y);
        const theta = 2 * Math.PI * i / goldenRatio;
        points.push([Math.cos(theta) * radius, y, Math.sin(theta) * radius]);
    }

    for (let iter = 0; iter < maxIter; iter++) {
        const gradients = points.map(() => [0, 0, 0]);

        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                const dx = points[i][0] - points[j][0];
                const dy = points[i][1] - points[j][1];
                const dz = points[i][2] - points[j][2];
                const dist2 = dx * dx + dy * dy + dz * dz;
                const dist = Math.sqrt(dist2);
                if (dist > 1e-10) {
                    const factor = 1 / (dist2 * dist);
                    gradients[i][0] -= factor * dx;
                    gradients[i][1] -= factor * dy;
                    gradients[i][2] -= factor * dz;
                    gradients[j][0] += factor * dx;
                    gradients[j][1] += factor * dy;
                    gradients[j][2] += factor * dz;
                }
            }
        }

        for (let i = 0; i < n; i++) {
            points[i][0] -= lr * gradients[i][0];
            points[i][1] -= lr * gradients[i][1];
            points[i][2] -= lr * gradients[i][2];

            const norm = Math.sqrt(points[i][0] ** 2 + points[i][1] ** 2 + points[i][2] ** 2);
            if (norm > 1e-10) {
                points[i][0] /= norm;
                points[i][1] /= norm;
                points[i][2] /= norm;
            }
        }

        if (iter > 50) lr *= 0.995;
    }

    THOMSON_CACHE[n] = points;
    return points;
}

// 【精确复制 physics.js 的 buildDirectionMatrix】
function buildDirectionMatrix(directions, orbitalParams) {
    const nDirs = directions.length;
    const nOrbitals = orbitalParams.length;
    const A = [];

    for (let i = 0; i < nDirs; i++) {
        const d = directions[i];
        const norm = Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        const x = d[0] / norm, y = d[1] / norm, z = d[2] / norm;

        const row = [];
        for (let j = 0; j < nOrbitals; j++) {
            // 【关键】使用 .angKey 结构，与 physics.js 一致
            const { l, m, t } = orbitalParams[j].angKey || { l: 0, m: 0, t: 'c' };

            let val = 0;
            if (l === 0) {
                val = 1 / Math.sqrt(4 * PI);
            } else if (l === 1) {
                const norm_p = Math.sqrt(3 / (4 * PI));
                if (m === 0) val = norm_p * z;            // pz
                else if (t === 'c') val = norm_p * x;      // px
                else val = norm_p * y;                     // py
            } else if (l === 2) {
                if (m === 0) {
                    val = Math.sqrt(5 / (16 * PI)) * (3 * z * z - 1);  // dz²
                } else if (m === 1) {
                    if (t === 'c') val = Math.sqrt(15 / (4 * PI)) * x * z;  // dxz
                    else val = Math.sqrt(15 / (4 * PI)) * y * z;            // dyz
                } else if (m === 2) {
                    if (t === 'c') val = Math.sqrt(15 / (16 * PI)) * (x * x - y * y);  // dx²-y²
                    else val = Math.sqrt(15 / (4 * PI)) * x * y;                       // dxy
                }
            }
            row.push(val);
        }
        A.push(row);
    }
    return A;
}

// Jacobi SVD
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
    const MAX_ITER = 50;

    for (let iter = 0; iter < MAX_ITER; iter++) {
        let maxError = 0;

        for (let i = 0; i < n - 1; i++) {
            for (let j = i + 1; j < n; j++) {
                let alpha = 0, beta = 0, gamma = 0;

                for (let k = 0; k < m; k++) {
                    alpha += U[k][i] * U[k][i];
                    beta += U[k][j] * U[k][j];
                    gamma += U[k][i] * U[k][j];
                }

                maxError = Math.max(maxError, Math.abs(gamma) / Math.sqrt(alpha * beta + EPSILON));

                if (Math.abs(gamma) < EPSILON) continue;

                let zeta = (beta - alpha) / (2 * gamma);
                let t = Math.sign(zeta) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
                let c = 1 / Math.sqrt(1 + t * t);
                let s = c * t;

                for (let k = 0; k < m; k++) {
                    let t1 = U[k][i];
                    let t2 = U[k][j];
                    U[k][i] = c * t1 - s * t2;
                    U[k][j] = s * t1 + c * t2;
                }

                for (let k = 0; k < n; k++) {
                    let t1 = V[k][i];
                    let t2 = V[k][j];
                    V[k][i] = c * t1 - s * t2;
                    V[k][j] = s * t1 + c * t2;
                }
            }
        }

        if (maxError < 1e-10) break;
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

function matMul(A, B) {
    const m = A.length;
    const n = A[0].length;
    const p = B[0].length;
    const C = [];
    for (let i = 0; i < m; i++) {
        const row = new Array(p).fill(0);
        for (let j = 0; j < p; j++) {
            for (let k = 0; k < n; k++) {
                row[j] += A[i][k] * B[k][j];
            }
        }
        C.push(row);
    }
    return C;
}

function matTranspose(A) {
    const m = A.length;
    const n = A[0].length;
    const AT = [];
    for (let i = 0; i < n; i++) {
        const row = new Array(m).fill(0);
        for (let j = 0; j < m; j++) {
            row[j] = A[j][i];
        }
        AT.push(row);
    }
    return AT;
}

function computeHybridCoefficients(directions, orbitalParams) {
    const A = buildDirectionMatrix(directions, orbitalParams);
    const { U, V } = jacobiSVD(A);
    const VT = matTranspose(V);
    const Q = matMul(U, VT);
    return Q;
}

function getHybridCoefficients(numOrbitals, orbitalParams) {
    if (!orbitalParams || orbitalParams.length !== numOrbitals) {
        throw new Error(`Missing orbitalParams`);
    }
    const directions = optimizeThomson(numOrbitals);
    return computeHybridCoefficients(directions, orbitalParams);
}

// 【精确复制 physics.js 的 sortOrbitalsForHybridization】
function sortOrbitalsForHybridization(paramsList) {
    if (!paramsList) return [];

    return [...paramsList].sort((a, b) => {
        // 1. 按 l 升序 (s < p < d < f)
        if (a.l !== b.l) return a.l - b.l;

        // 2. 按 n 升序
        if (a.n !== b.n) return a.n - b.n;

        // 3. 按 m 和 t 排序 (特定顺序)
        // 对于 p (l=1): px(m=1,c), py(m=1,s), pz(m=0)
        if (a.l === 1) {
            const score = (p) => {
                if (p.angKey.m === 1 && p.angKey.t === 'c') return 1; // px
                if (p.angKey.m === 1 && p.angKey.t === 's') return 2; // py
                if (p.angKey.m === 0) return 3; // pz
                return 4;
            };
            return score(a) - score(b);
        }

        // 对于 d (l=2): dx2-y2(m=2,c), dz2(m=0), ...
        if (a.l === 2) {
            const score = (p) => {
                if (p.angKey.m === 2 && p.angKey.t === 'c') return 1; // dx2-y2
                if (p.angKey.m === 0) return 2; // dz2
                return 3 + p.angKey.m;
            };
            return score(a) - score(b);
        }

        return 0;
    });
}

// ==================== 测试用例 ====================

// 【关键】轨道参数定义 - 使用 angKey 结构，与 physics.js 一致
function makeOrbital(n, l, m, t) {
    return {
        n: n,
        l: l,
        angKey: { l: l, m: m, t: t }
    };
}

const ORBITALS = {
    '2s': makeOrbital(2, 0, 0, 'c'),
    '2px': makeOrbital(2, 1, 1, 'c'),
    '2py': makeOrbital(2, 1, 1, 's'),
    '2pz': makeOrbital(2, 1, 0, 'c'),
    '3dz2': makeOrbital(3, 2, 0, 'c'),
    '3dx2y2': makeOrbital(3, 2, 2, 'c'),
    '3dxy': makeOrbital(3, 2, 2, 's'),
    '3dxz': makeOrbital(3, 2, 1, 'c'),
    '3dyz': makeOrbital(3, 2, 1, 's')
};

// 验证正交性
function checkOrthogonality(matrix) {
    const n = matrix.length;
    let maxError = 0;

    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            let dot = 0;
            for (let k = 0; k < n; k++) {
                dot += matrix[i][k] * matrix[j][k];
            }
            const expected = (i === j) ? 1 : 0;
            maxError = Math.max(maxError, Math.abs(dot - expected));
        }
    }
    return maxError;
}

// 计算 s 轨道贡献（第一列系数的平方和）
function computeSContributions(matrix) {
    const n = matrix.length;
    const contributions = [];
    for (let i = 0; i < n; i++) {
        contributions.push(matrix[i][0] * matrix[i][0]);
    }
    return contributions;
}

// 测试杂化类型
const TEST_CASES = [
    { name: 'sp', orbitals: ['2s', '2pz'] },
    { name: 'sp2', orbitals: ['2s', '2px', '2py'] },
    { name: 'sp3', orbitals: ['2s', '2px', '2py', '2pz'] },
    { name: 'sp3d', orbitals: ['2s', '2px', '2py', '2pz', '3dz2'] },
    { name: 'sp3d2', orbitals: ['2s', '2px', '2py', '2pz', '3dz2', '3dx2y2'] }
];

console.log('========================================');
console.log('杂化轨道系数算法测试 (与 physics.js 一致)');
console.log('========================================\n');

let allPassed = true;

for (const testCase of TEST_CASES) {
    console.log(`--- ${testCase.name} (${testCase.orbitals.join(', ')}) ---`);

    // 获取轨道参数（未排序）
    let orbitalParams = testCase.orbitals.map(name => ORBITALS[name]);

    // 【关键】调用排序函数，与 physics.js 一致
    const sortedParams = sortOrbitalsForHybridization(orbitalParams);
    console.log(`排序后顺序: ${sortedParams.map(p => {
        if (p.l === 0) return 's';
        if (p.l === 1) {
            if (p.angKey.m === 0) return 'pz';
            if (p.angKey.t === 'c') return 'px';
            return 'py';
        }
        if (p.l === 2) {
            if (p.angKey.m === 0) return 'dz2';
            if (p.angKey.m === 2 && p.angKey.t === 'c') return 'dx2-y2';
            if (p.angKey.m === 2 && p.angKey.t === 's') return 'dxy';
            if (p.angKey.m === 1 && p.angKey.t === 'c') return 'dxz';
            if (p.angKey.m === 1 && p.angKey.t === 's') return 'dyz';
        }
        return '?';
    }).join(', ')}`);

    const n = sortedParams.length;

    try {
        const coeffMatrix = getHybridCoefficients(n, sortedParams);

        // 检查正交性
        const orthError = checkOrthogonality(coeffMatrix);
        const orthOK = orthError < 1e-6;
        console.log(`正交性误差: ${orthError.toExponential(3)} ${orthOK ? '✓' : '✗'}`);

        // 检查 s 贡献
        const sContributions = computeSContributions(coeffMatrix);
        console.log(`s轨道贡献: ${sContributions.map(x => (x * 100).toFixed(1) + '%').join(', ')}`);

        // 检查 s 贡献总和是否为 1
        const sSum = sContributions.reduce((a, b) => a + b, 0);
        console.log(`s贡献总和: ${(sSum * 100).toFixed(1)}% ${Math.abs(sSum - 1) < 0.01 ? '✓' : '✗'}`);

        // 打印系数矩阵
        console.log('系数矩阵:');
        for (let i = 0; i < n; i++) {
            console.log(`  杂化轨道${i}: [${coeffMatrix[i].map(x => x.toFixed(4)).join(', ')}]`);
        }

        if (!orthOK) allPassed = false;
        if (Math.abs(sSum - 1) > 0.01) allPassed = false;

    } catch (e) {
        console.log(`错误: ${e.message}`);
        console.log(e.stack);
        allPassed = false;
    }

    console.log('');
}

// 检查 Thomson 方向
console.log('--- Thomson 方向验证 ---');
for (const testCase of TEST_CASES) {
    const n = testCase.orbitals.length;
    const dirs = optimizeThomson(n);
    const energy = thomsonEnergy(dirs);

    // 理论最优能量
    const optimalEnergies = { 2: 0.5, 3: 1.732, 4: 3.674, 5: 6.475, 6: 9.985 };
    const optimal = optimalEnergies[n] || 0;
    const diff = Math.abs(energy - optimal);

    console.log(`${testCase.name} (n=${n}): 能量=${energy.toFixed(4)}, 理论=${optimal.toFixed(3)}, 偏差=${diff.toFixed(4)} ${diff < 0.01 ? '✓' : '?'}`);
}

console.log('\n========================================');
console.log(allPassed ? '所有测试通过 ✓' : '存在测试失败 ✗');
console.log('========================================');
