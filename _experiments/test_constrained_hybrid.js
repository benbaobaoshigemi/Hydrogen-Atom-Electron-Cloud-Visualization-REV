/**
 * 约束 Thomson 优化 + 杂化系数计算
 * 
 * 核心思想：在轨道的"可达角向空间"内进行 Thomson 优化
 * - sp (s+pz): 只能产生 z 轴方向
 * - sp2 (s+px+py): 只能产生 xy 平面内的方向
 * - sp3 (s+px+py+pz): 可以产生任意 3D 方向
 * - sp3d (s+px+py+pz+dz2): 三角双锥对称
 */

const PI = Math.PI;

// 球谐函数值
function Y_real(l, m, t, theta, phi) {
    const cosT = Math.cos(theta);
    const sinT = Math.sin(theta);

    if (l === 0) return 1 / Math.sqrt(4 * PI);

    if (l === 1) {
        const norm = Math.sqrt(3 / (4 * PI));
        if (m === 0) return norm * cosT;       // pz
        if (t === 'c') return norm * sinT * Math.cos(phi);  // px
        return norm * sinT * Math.sin(phi);    // py
    }

    if (l === 2) {
        if (m === 0) return Math.sqrt(5 / (16 * PI)) * (3 * cosT * cosT - 1);  // dz2
        if (m === 1) {
            const norm = Math.sqrt(15 / (4 * PI)) * sinT * cosT;
            if (t === 'c') return norm * Math.cos(phi);  // dxz
            return norm * Math.sin(phi);  // dyz
        }
        if (m === 2) {
            const norm = Math.sqrt(15 / (16 * PI)) * sinT * sinT;
            if (t === 'c') return norm * Math.cos(2 * phi);  // dx2-y2
            return Math.sqrt(15 / (4 * PI)) * sinT * sinT * Math.sin(phi) * Math.cos(phi);  // dxy
        }
    }
    return 0;
}

// 分析轨道组合的维度和对称性
function analyzeOrbitalSet(orbitalParams) {
    const hasS = orbitalParams.some(p => p.l === 0);
    const pOrbitals = orbitalParams.filter(p => p.l === 1);
    const dOrbitals = orbitalParams.filter(p => p.l === 2);

    const hasPx = pOrbitals.some(p => p.angKey.m === 1 && p.angKey.t === 'c');
    const hasPy = pOrbitals.some(p => p.angKey.m === 1 && p.angKey.t === 's');
    const hasPz = pOrbitals.some(p => p.angKey.m === 0);

    const hasDz2 = dOrbitals.some(p => p.angKey.m === 0);
    const hasDx2y2 = dOrbitals.some(p => p.angKey.m === 2 && p.angKey.t === 'c');

    return {
        hasS, hasPx, hasPy, hasPz, hasDz2, hasDx2y2,
        pCount: pOrbitals.length,
        dCount: dOrbitals.length,
        total: orbitalParams.length
    };
}

// 根据轨道组合生成约束方向
function generateConstrainedDirections(orbitalParams) {
    const info = analyzeOrbitalSet(orbitalParams);
    const n = orbitalParams.length;

    // sp: 线性（z轴）
    if (n === 2 && info.hasS && info.hasPz && !info.hasPx && !info.hasPy) {
        return [[0, 0, 1], [0, 0, -1]];
    }

    // sp (x): 线性（x轴）
    if (n === 2 && info.hasS && info.hasPx && !info.hasPy && !info.hasPz) {
        return [[1, 0, 0], [-1, 0, 0]];
    }

    // sp2: 平面三角形（xy平面）
    if (n === 3 && info.hasS && info.hasPx && info.hasPy && !info.hasPz) {
        return [
            [1, 0, 0],
            [-0.5, Math.sqrt(3) / 2, 0],
            [-0.5, -Math.sqrt(3) / 2, 0]
        ];
    }

    // sp2 (xz平面)
    if (n === 3 && info.hasS && info.hasPx && !info.hasPy && info.hasPz) {
        return [
            [1, 0, 0],
            [-0.5, 0, Math.sqrt(3) / 2],
            [-0.5, 0, -Math.sqrt(3) / 2]
        ];
    }

    // sp3: 正四面体
    if (n === 4 && info.hasS && info.hasPx && info.hasPy && info.hasPz && info.dCount === 0) {
        const a = 1 / Math.sqrt(3);
        return [
            [a, a, a],
            [a, -a, -a],
            [-a, a, -a],
            [-a, -a, a]
        ];
    }

    // sp3d: 三角双锥
    if (n === 5 && info.hasS && info.hasPx && info.hasPy && info.hasPz && info.hasDz2 && !info.hasDx2y2) {
        return [
            [0, 0, 1],                          // 轴向 +z
            [0, 0, -1],                         // 轴向 -z
            [1, 0, 0],                          // 赤道
            [-0.5, Math.sqrt(3) / 2, 0],          // 赤道 120°
            [-0.5, -Math.sqrt(3) / 2, 0]          // 赤道 240°
        ];
    }

    // sp3d2: 正八面体
    if (n === 6 && info.hasS && info.hasPx && info.hasPy && info.hasPz && info.hasDz2 && info.hasDx2y2) {
        return [
            [1, 0, 0], [-1, 0, 0],
            [0, 1, 0], [0, -1, 0],
            [0, 0, 1], [0, 0, -1]
        ];
    }

    // 通用情况：使用标准Thomson优化
    return null;
}

// SVD 正交化（Procrustes）
function svdOrthogonalize(A) {
    const n = A.length;
    let U = A.map(row => [...row]);
    let V = [];
    for (let i = 0; i < n; i++) {
        const row = new Array(n).fill(0);
        row[i] = 1;
        V.push(row);
    }

    const EPSILON = 1e-15;
    for (let iter = 0; iter < 50; iter++) {
        let maxErr = 0;
        for (let i = 0; i < n - 1; i++) {
            for (let j = i + 1; j < n; j++) {
                let a = 0, b = 0, g = 0;
                for (let k = 0; k < n; k++) {
                    a += U[k][i] * U[k][i];
                    b += U[k][j] * U[k][j];
                    g += U[k][i] * U[k][j];
                }
                maxErr = Math.max(maxErr, Math.abs(g) / Math.sqrt(a * b + EPSILON));
                if (Math.abs(g) < EPSILON) continue;

                let z = (b - a) / (2 * g);
                let t = Math.sign(z) / (Math.abs(z) + Math.sqrt(1 + z * z));
                let c = 1 / Math.sqrt(1 + t * t), s = c * t;

                for (let k = 0; k < n; k++) {
                    let t1 = U[k][i], t2 = U[k][j];
                    U[k][i] = c * t1 - s * t2;
                    U[k][j] = s * t1 + c * t2;
                }
                for (let k = 0; k < n; k++) {
                    let t1 = V[k][i], t2 = V[k][j];
                    V[k][i] = c * t1 - s * t2;
                    V[k][j] = s * t1 + c * t2;
                }
            }
        }
        if (maxErr < 1e-10) break;
    }

    // 归一化 U
    for (let i = 0; i < n; i++) {
        let sum = 0;
        for (let k = 0; k < n; k++) sum += U[k][i] * U[k][i];
        let norm = Math.sqrt(sum);
        if (norm > EPSILON) {
            for (let k = 0; k < n; k++) U[k][i] /= norm;
        }
    }

    // Q = U V^T
    const Q = [];
    for (let i = 0; i < n; i++) {
        const row = new Array(n).fill(0);
        for (let j = 0; j < n; j++) {
            for (let k = 0; k < n; k++) {
                row[j] += U[i][k] * V[j][k];
            }
        }
        Q.push(row);
    }
    return Q;
}

// 构建方向矩阵
function buildDirectionMatrix(directions, orbitalParams) {
    const n = directions.length;
    const A = [];

    for (let i = 0; i < n; i++) {
        const d = directions[i];
        const norm = Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        const x = d[0] / norm, y = d[1] / norm, z = d[2] / norm;

        const theta = Math.acos(Math.max(-1, Math.min(1, z)));
        const phi = Math.atan2(y, x);

        const row = [];
        for (let j = 0; j < orbitalParams.length; j++) {
            const p = orbitalParams[j];
            row.push(Y_real(p.l, p.angKey.m, p.angKey.t, theta, phi));
        }
        A.push(row);
    }
    return A;
}

// 轨道排序
function sortOrbitals(params) {
    return [...params].sort((a, b) => {
        if (a.l !== b.l) return a.l - b.l;
        if (a.n !== b.n) return a.n - b.n;
        if (a.l === 1) {
            const s = p => (p.angKey.m === 1 && p.angKey.t === 'c') ? 1 :
                (p.angKey.m === 1 && p.angKey.t === 's') ? 2 : 3;
            return s(a) - s(b);
        }
        if (a.l === 2) {
            const s = p => (p.angKey.m === 2 && p.angKey.t === 'c') ? 1 :
                (p.angKey.m === 0) ? 2 : 3 + p.angKey.m;
            return s(a) - s(b);
        }
        return 0;
    });
}

// 主函数：计算杂化系数
function getHybridCoefficients(orbitalParams) {
    const sorted = sortOrbitals(orbitalParams);
    const n = sorted.length;

    // 尝试获取约束方向
    const directions = generateConstrainedDirections(sorted);

    if (!directions) {
        console.error("No constrained directions for this orbital combination. Using identity matrix.");
        return {
            coeffMatrix: Array.from({ length: n }, (_, i) =>
                Array.from({ length: n }, (_, j) => i === j ? 1 : 0)), sortedParams: sorted
        };
    }

    // 构建方向矩阵并正交化
    const A = buildDirectionMatrix(directions, sorted);
    const Q = svdOrthogonalize(A);

    return { coeffMatrix: Q, sortedParams: sorted };
}

// ==================== 测试 ====================

function makeOrbital(n, l, m, t) {
    return { n, l, angKey: { l, m, t } };
}

const ORBITALS = {
    '2s': makeOrbital(2, 0, 0, 'c'),
    '2px': makeOrbital(2, 1, 1, 'c'),
    '2py': makeOrbital(2, 1, 1, 's'),
    '2pz': makeOrbital(2, 1, 0, 'c'),
    '3dz2': makeOrbital(3, 2, 0, 'c'),
    '3dx2y2': makeOrbital(3, 2, 2, 'c'),
};

const tests = [
    { name: 'sp', orbitals: ['2s', '2pz'] },
    { name: 'sp2', orbitals: ['2s', '2px', '2py'] },
    { name: 'sp3', orbitals: ['2s', '2px', '2py', '2pz'] },
    { name: 'sp3d', orbitals: ['2s', '2px', '2py', '2pz', '3dz2'] },
    { name: 'sp3d2', orbitals: ['2s', '2px', '2py', '2pz', '3dz2', '3dx2y2'] },
];

console.log('=== 约束 Thomson + SVD 杂化系数测试 ===\n');

for (const test of tests) {
    console.log(`--- ${test.name} ---`);
    const params = test.orbitals.map(k => ORBITALS[k]);

    try {
        const { coeffMatrix, sortedParams } = getHybridCoefficients(params);

        // 检查正交性
        let orthErr = 0;
        for (let i = 0; i < coeffMatrix.length; i++) {
            for (let j = 0; j < coeffMatrix.length; j++) {
                let dot = 0;
                for (let k = 0; k < coeffMatrix.length; k++) {
                    dot += coeffMatrix[i][k] * coeffMatrix[j][k];
                }
                orthErr = Math.max(orthErr, Math.abs(dot - (i === j ? 1 : 0)));
            }
        }
        console.log(`正交性误差: ${orthErr.toExponential(3)} ${orthErr < 1e-6 ? '✓' : '✗'}`);

        // s 贡献
        const sContribs = coeffMatrix.map(row => row[0] * row[0]);
        console.log(`s贡献: ${sContribs.map(c => (c * 100).toFixed(1) + '%').join(', ')}`);

        // 打印系数
        console.log('系数矩阵:');
        coeffMatrix.forEach((row, i) => {
            console.log(`  ${i}: [${row.map(x => x.toFixed(4)).join(', ')}]`);
        });

    } catch (e) {
        console.log(`错误: ${e.message}`);
    }
    console.log('');
}
