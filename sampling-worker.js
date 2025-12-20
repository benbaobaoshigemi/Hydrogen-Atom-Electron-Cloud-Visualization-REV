/**
 * Web Worker for Monte Carlo Sampling
 * 将计算密集型的蒙特卡洛采样移到后台线程，避免阻塞主线程渲染
 */

// 导入 Slater 基组数据
importScripts('slater_basis.js');

// ==================== 物理计算函数（从 physics.js 复制，Worker 无法访问主线程代码）====================

const A0 = 1; // Bohr radius unit
const PI = Math.PI;

// ==================== Thomson + SVD 杂化系数算法 ====================

const THOMSON_CACHE_WORKER = {};

function optimizeThomsonWorker(n, maxIter = 200, lr = 0.1) {
    if (THOMSON_CACHE_WORKER[n]) return THOMSON_CACHE_WORKER[n];

    if (n === 1) {
        THOMSON_CACHE_WORKER[n] = [[0, 0, 1]];
        return THOMSON_CACHE_WORKER[n];
    }
    if (n === 2) {
        THOMSON_CACHE_WORKER[n] = [[0, 0, 1], [0, 0, -1]];
        return THOMSON_CACHE_WORKER[n];
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
    THOMSON_CACHE_WORKER[n] = points;
    return points;
}

function buildDirectionMatrixWorker(directions, orbitalParams) {
    const nDirs = directions.length;
    const nOrbitals = orbitalParams.length;
    const A = [];
    for (let i = 0; i < nDirs; i++) {
        const d = directions[i];
        const norm = Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        const x = d[0] / norm, y = d[1] / norm, z = d[2] / norm;
        const row = [];
        for (let j = 0; j < nOrbitals; j++) {
            const { l, m, t } = orbitalParams[j].angKey || { l: 0, m: 0, t: 'c' };
            let val = 0;
            if (l === 0) {
                val = 1 / Math.sqrt(4 * PI);
            } else if (l === 1) {
                const norm_p = Math.sqrt(3 / (4 * PI));
                if (m === 0) val = norm_p * z;
                else if (t === 'c') val = norm_p * x;
                else val = norm_p * y;
            } else if (l === 2) {
                if (m === 0) {
                    val = Math.sqrt(5 / (16 * PI)) * (3 * z * z - 1);
                } else if (m === 1) {
                    if (t === 'c') val = Math.sqrt(15 / (4 * PI)) * x * z;
                    else val = Math.sqrt(15 / (4 * PI)) * y * z;
                } else if (m === 2) {
                    if (t === 'c') val = Math.sqrt(15 / (16 * PI)) * (x * x - y * y);
                    else val = Math.sqrt(15 / (4 * PI)) * x * y;
                }
            }
            row.push(val);
        }
        A.push(row);
    }
    return A;
}

function jacobiSVDWorker(A) {
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

function matMulWorker(A, B) {
    const m = A.length, n = A[0].length, p = B[0].length;
    const C = [];
    for (let i = 0; i < m; i++) {
        const row = new Array(p).fill(0);
        for (let j = 0; j < p; j++) {
            for (let k = 0; k < n; k++) row[j] += A[i][k] * B[k][j];
        }
        C.push(row);
    }
    return C;
}

function matTransposeWorker(A) {
    const m = A.length, n = A[0].length;
    const AT = [];
    for (let i = 0; i < n; i++) {
        const row = new Array(m).fill(0);
        for (let j = 0; j < m; j++) row[j] = A[j][i];
        AT.push(row);
    }
    return AT;
}

function sortOrbitalsForHybridizationWorker(paramsList) {
    if (!paramsList) return [];
    return [...paramsList].sort((a, b) => {
        if (a.l !== b.l) return a.l - b.l;
        if (a.n !== b.n) return a.n - b.n;
        if (a.l === 1) {
            const score = (p) => {
                if (p.angKey.m === 1 && p.angKey.t === 'c') return 1;
                if (p.angKey.m === 1 && p.angKey.t === 's') return 2;
                if (p.angKey.m === 0) return 3;
                return 4;
            };
            return score(a) - score(b);
        }
        if (a.l === 2) {
            const score = (p) => {
                if (p.angKey.m === 2 && p.angKey.t === 'c') return 1;
                if (p.angKey.m === 0) return 2;
                return 3 + p.angKey.m;
            };
            return score(a) - score(b);
        }
        return 0;
    });
}

const HYBRID_COEFF_CACHE_WORKER = {};

// 根据轨道组合生成约束方向
function generateConstrainedDirectionsWorker(orbitalParams) {
    const n = orbitalParams.length;

    const hasS = orbitalParams.some(p => p.l === 0);
    const pOrbitals = orbitalParams.filter(p => p.l === 1);
    const dOrbitals = orbitalParams.filter(p => p.l === 2);

    const hasPx = pOrbitals.some(p => p.angKey.m === 1 && p.angKey.t === 'c');
    const hasPy = pOrbitals.some(p => p.angKey.m === 1 && p.angKey.t === 's');
    const hasPz = pOrbitals.some(p => p.angKey.m === 0);

    const hasDz2 = dOrbitals.some(p => p.angKey.m === 0);
    const hasDx2y2 = dOrbitals.some(p => p.angKey.m === 2 && p.angKey.t === 'c');

    // sp
    if (n === 2 && hasS && pOrbitals.length === 1) {
        if (hasPz) return [[0, 0, 1], [0, 0, -1]];
        if (hasPx) return [[1, 0, 0], [-1, 0, 0]];
        if (hasPy) return [[0, 1, 0], [0, -1, 0]];
    }

    // sp2
    if (n === 3 && hasS && pOrbitals.length === 2) {
        if (hasPx && hasPy) return [[1, 0, 0], [-0.5, Math.sqrt(3) / 2, 0], [-0.5, -Math.sqrt(3) / 2, 0]];
        if (hasPx && hasPz) return [[1, 0, 0], [-0.5, 0, Math.sqrt(3) / 2], [-0.5, 0, -Math.sqrt(3) / 2]];
        if (hasPy && hasPz) return [[0, 1, 0], [0, -0.5, Math.sqrt(3) / 2], [0, -0.5, -Math.sqrt(3) / 2]];
    }

    // sp3
    if (n === 4 && hasS && pOrbitals.length === 3 && hasPx && hasPy && hasPz && dOrbitals.length === 0) {
        const a = 1 / Math.sqrt(3);
        return [[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]];
    }

    // sp3d
    if (n === 5 && hasS && hasPx && hasPy && hasPz && hasDz2 && !hasDx2y2) {
        return [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-0.5, Math.sqrt(3) / 2, 0], [-0.5, -Math.sqrt(3) / 2, 0]];
    }

    // sp3d2
    if (n === 6 && hasS && hasPx && hasPy && hasPz && hasDz2 && hasDx2y2) {
        return [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
    }

    return optimizeThomsonWorker(n);
}

function getHybridCoefficientsWorker(orbitalParams) {
    if (!orbitalParams || orbitalParams.length === 0) return null;

    // 排序轨道
    const sortedParams = sortOrbitalsForHybridizationWorker(orbitalParams);
    const n = sortedParams.length;

    // 缓存键
    const cacheKey = sortedParams.map(p => `${p.n}_${p.l}_${p.angKey.m}_${p.angKey.t}`).join('|');
    if (HYBRID_COEFF_CACHE_WORKER[cacheKey]) {
        return { coeffMatrix: HYBRID_COEFF_CACHE_WORKER[cacheKey], sortedParams };
    }

    // 【核心修复】使用约束方向
    const directions = generateConstrainedDirectionsWorker(sortedParams);
    const A = buildDirectionMatrixWorker(directions, sortedParams);
    const { U, V } = jacobiSVDWorker(A);
    const VT = matTransposeWorker(V);
    const Q = matMulWorker(U, VT);

    HYBRID_COEFF_CACHE_WORKER[cacheKey] = Q;
    return { coeffMatrix: Q, sortedParams };
}

// 阶乘表（预计算）
const FACT = (() => {
    const f = [1];
    for (let i = 1; i <= 64; i++) f[i] = f[i - 1] * i;
    return f;
})();

function factorial(n) {
    return FACT[n] ?? Infinity;
}

function binomialInt(n, k) {
    if (k < 0 || k > n) return 0;
    return factorial(n) / (factorial(k) * factorial(n - k));
}

// 广义 Laguerre 多项式
function generalizedLaguerre(k, alpha, x) {
    let sum = 0;
    for (let i = 0; i <= k; i++) {
        const c = ((i % 2) ? -1 : 1) * binomialInt(k + alpha, k - i) / factorial(i);
        sum += c * Math.pow(x, i);
    }
    return sum;
}

// 连带 Legendre 多项式 (Condon-Shortley 相位)
function associatedLegendre(l, m, x) {
    const mm = Math.abs(m);
    if (l < mm) return 0;
    let pmm = 1.0;
    if (mm > 0) {
        const s = Math.sqrt(Math.max(0, 1 - x * x));
        let fact = 1.0;
        for (let i = 1; i <= mm; i++) {
            pmm *= -fact * s;
            fact += 2.0;
        }
    }
    if (l === mm) return pmm;
    let pmmp1 = x * (2 * mm + 1) * pmm;
    if (l === mm + 1) return pmmp1;
    let pll = 0;
    for (let ll = mm + 2; ll <= l; ll++) {
        pll = ((2 * ll - 1) * x * pmmp1 - (ll + mm - 1) * pmm) / (ll - mm);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

// 复球谐函数
function Ylm_complex(l, m, theta, phi) {
    const mm = Math.abs(m);
    const Plm = associatedLegendre(l, mm, Math.cos(theta));
    const N = Math.sqrt(((2 * l + 1) / (4 * PI)) * (factorial(l - mm) / factorial(l + mm)));
    const base = N * Plm;
    if (m === 0) {
        return { re: base, im: 0 };
    }
    const cos_mphi = Math.cos(mm * phi);
    const sin_mphi = Math.sin(mm * phi);
    if (m > 0) {
        return { re: base * cos_mphi, im: base * sin_mphi };
    } else {
        const sign = (mm % 2) ? -1 : 1;
        return { re: sign * base * cos_mphi, im: -sign * base * sin_mphi };
    }
}

function Ylm_abs2(l, m, theta, phi) {
    const y = Ylm_complex(l, m, theta, phi);
    return y.re * y.re + y.im * y.im;
}

// 实球谐函数（模方）
// 【关键修复】复用 realYlm_value 确保相位修正一致性
// 与 physics.js 同步，防止采样/可视化方向不一致
function realYlm_abs2(l, m, type, theta, phi) {
    const val = realYlm_value(l, m, type, theta, phi);
    return val * val;
}

// 实球谐函数（值）
// Chemistry Convention: Positive lobes align with positive axes
// 【关键修复】抵消 Condon-Shortley 相位，与 physics.js 保持同步
function realYlm_value(l, m, type, theta, phi) {
    if (m === 0) {
        const y = Ylm_complex(l, 0, theta, phi);
        return y.re;
    }
    const mm = Math.abs(m);
    const y = Ylm_complex(l, mm, theta, phi);

    // 【修复】抵消 CS 相位：对于奇数 m，翻转符号
    const csPhaseCorrection = (mm % 2 === 1) ? -1 : 1;

    if (type === 'c') {
        return csPhaseCorrection * Math.SQRT2 * y.re;
    } else {
        return csPhaseCorrection * Math.SQRT2 * y.im;
    }
}

// Helper: Convert n, l to orbital key (e.g., 2,1 -> '2p')
function getOrbitalKey(n, l) {
    const subshells = ['s', 'p', 'd', 'f', 'g'];
    const letter = subshells[l] || '';
    return n + letter;
}

// Slater Type Orbital Radial Function
function slaterRadialR(basis, r) {
    if (!basis) return 0;
    let sum = 0;
    for (let i = 0; i < basis.length; i++) {
        const { nStar, zeta, coeff } = basis[i];
        // Normalization constant N = (2zeta)^(n*+0.5) / sqrt((2n*)!)
        const nFact2 = factorial(2 * nStar);
        const norm = Math.pow(2 * zeta, nStar + 0.5) / Math.sqrt(nFact2);
        // Radial part: r^(n*-1) * e^(-zeta*r)
        const val = coeff * norm * Math.pow(r, nStar - 1) * Math.exp(-zeta * r);
        sum += val;
    }
    return sum;
}

// 归一化径向函数 R_nl(r)
function radialR(n, l, r, Z = 1, a0 = A0, atomType = 'H') {
    // Strategy A: Slater Type Orbitals (STO)
    // Strategy A: Slater Type Orbitals (STO)
    if (atomType && atomType !== 'H' && (self.SlaterBasis || (typeof window !== 'undefined' && window.SlaterBasis))) {
        const orbitalKey = getOrbitalKey(n, l);
        const slater = self.SlaterBasis || window.SlaterBasis;
        const atomData = slater[atomType];

        if (atomData && atomData.orbitals && atomData.orbitals[orbitalKey]) {
            const basis = atomData.orbitals[orbitalKey];
            const val = slaterRadialR(basis, r);

            // 【智能相位校正】
            // 目标：确保 STO 波函数在远核区域（大 r）的符号与氢原子解析解一致
            // 氢原子解析解 R_nl(r) 在大 r 处的符号由 Laguerre 多项式的最高次项决定 -> (-1)^(n-l-1)
            const analyticalSign = ((n - l - 1) % 2 === 0) ? 1 : -1;

            // 寻找 STO 在大 r 处的主导项（zeta 最小的项，即衰减最慢的项）
            let dominantTerm = basis[0];
            for (let i = 1; i < basis.length; i++) {
                const term = basis[i];
                // 优先找 zeta 更小的
                // 如果 zeta 相同（极少见），找 nStar 更大的（r 的幂次更高）
                if (term.zeta < dominantTerm.zeta ||
                    (Math.abs(term.zeta - dominantTerm.zeta) < 1e-9 && term.nStar > dominantTerm.nStar)) {
                    dominantTerm = term;
                }
            }

            const stoSign = Math.sign(dominantTerm.coeff);

            // 如果符号不一致，则翻转 STO 结果
            if (stoSign * analyticalSign < 0) {
                return -val;
            }
            return val;
        }
    }

    // Strategy B: Hydrogen analytical solution (Default)
    if (n <= 0 || l < 0 || l >= n) return 0;
    const rho = (2 * Z * r) / (n * a0);
    const k = n - l - 1;
    if (k < 0) return 0;
    const pref = Math.pow(2 * Z / (n * a0), 1.5) * Math.sqrt(factorial(n - l - 1) / (2 * n * factorial(n + l)));
    const poly = generalizedLaguerre(k, 2 * l + 1, rho);
    return pref * Math.exp(-rho / 2) * Math.pow(rho, l) * poly;
}

// 3D 概率密度
function density3D_real(angKey, n, l, r, theta, phi, Z = 1, a0 = A0, atomType = 'H') {
    const R = radialR(n, l, r, Z, a0, atomType);
    let Y2 = 1 / (4 * PI);
    if (angKey && typeof angKey === 'object') {
        Y2 = realYlm_abs2(angKey.l, angKey.m, angKey.t, theta, phi);
    }
    return (R * R) * Y2;
}

/**
 * 杂化轨道概率密度计算
 * 
 * 【物理原理】
 * 杂化轨道是多个原子轨道波函数的线性组合：Ψ_hybrid = Σ c_i × Ψ_i
 * 概率密度为 |Ψ_hybrid|² = |Σ c_i × R_i(r) × Y_i(θ,φ)|²
 * 
 * @param {Array} paramsList - 轨道参数列表
 * @param {number} r - 径向距离
 * @param {number} theta - 极角
 * @param {number} phi - 方位角
 * @returns {number} - 概率密度 |Ψ_hybrid|²
 */
// 杂化轨道概率密度计算
function hybridDensity3D(paramsList, r, theta, phi, Z = 1, a0 = A0, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return 0;

    let psiReal = 0;
    const defaultCoeff = 1.0 / Math.sqrt(paramsList.length);

    for (const p of paramsList) {
        const coeff = p.coefficient !== undefined ? p.coefficient : defaultCoeff;
        const R = radialR(p.n, p.l, r, Z, a0, atomType);
        const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
        psiReal += coeff * R * Y;
    }

    return psiReal * psiReal;
}

/**
 * 估算杂化轨道的最大概率密度（用于拒绝采样）
 */
function hybridEstimateMaxDensity(paramsList, rMax, numSamples = 1000, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return 1;

    let maxDensity = 0;

    for (let i = 0; i < numSamples; i++) {
        const r = Math.random() * rMax * Math.pow(Math.random(), 0.5);
        const cosTheta = 2 * Math.random() - 1;
        const theta = Math.acos(cosTheta);
        const phi = 2 * PI * Math.random();

        const density = hybridDensity3D(paramsList, r, theta, phi, 1, 1, atomType);
        if (density > maxDensity) {
            maxDensity = density;
        }
    }

    return maxDensity * 1.5;
}

// 轨道参数映射 - 包含 n=1 到 n=5 的所有轨道
function orbitalParamsFromKey(key) {
    const R = (n, l, m, t) => ({ n, l, angKey: { l, m, t } });
    switch (key) {
        // n=1
        case '1s': return R(1, 0, 0, 'c');
        // n=2
        case '2s': return R(2, 0, 0, 'c');
        case '2pz': return R(2, 1, 0, 'c');
        case '2px': return R(2, 1, 1, 'c');
        case '2py': return R(2, 1, 1, 's');
        // n=3
        case '3s': return R(3, 0, 0, 'c');
        case '3pz': return R(3, 1, 0, 'c');
        case '3px': return R(3, 1, 1, 'c');
        case '3py': return R(3, 1, 1, 's');
        case '3d_z2': return R(3, 2, 0, 'c');
        case '3d_xz': return R(3, 2, 1, 'c');
        case '3d_yz': return R(3, 2, 1, 's');
        case '3d_xy': return R(3, 2, 2, 's');
        case '3d_x2-y2': return R(3, 2, 2, 'c');
        // n=4
        case '4s': return R(4, 0, 0, 'c');
        case '4pz': return R(4, 1, 0, 'c');
        case '4px': return R(4, 1, 1, 'c');
        case '4py': return R(4, 1, 1, 's');
        case '4d_z2': return R(4, 2, 0, 'c');
        case '4d_xz': return R(4, 2, 1, 'c');
        case '4d_yz': return R(4, 2, 1, 's');
        case '4d_xy': return R(4, 2, 2, 's');
        case '4d_x2-y2': return R(4, 2, 2, 'c');
        case '4f_z3': return R(4, 3, 0, 'c');
        case '4f_xz2': return R(4, 3, 1, 'c');
        case '4f_yz2': return R(4, 3, 1, 's');
        case '4f_z(x2-y2)': return R(4, 3, 2, 'c');
        case '4f_xyz': return R(4, 3, 2, 's');
        case '4f_x(x2-3y2)': return R(4, 3, 3, 'c');
        case '4f_y(3x2-y2)': return R(4, 3, 3, 's');
        // n=5 轨道
        case '5s': return R(5, 0, 0, 'c');
        case '5pz': return R(5, 1, 0, 'c');
        case '5px': return R(5, 1, 1, 'c');
        case '5py': return R(5, 1, 1, 's');
        case '5d_z2': return R(5, 2, 0, 'c');
        case '5d_xz': return R(5, 2, 1, 'c');
        case '5d_yz': return R(5, 2, 1, 's');
        case '5d_xy': return R(5, 2, 2, 's');
        case '5d_x2-y2': return R(5, 2, 2, 'c');
        case '5f_z3': return R(5, 3, 0, 'c');
        case '5f_xz2': return R(5, 3, 1, 'c');
        case '5f_yz2': return R(5, 3, 1, 's');
        case '5f_z(x2-y2)': return R(5, 3, 2, 'c');
        case '5f_xyz': return R(5, 3, 2, 's');
        case '5f_x(x2-3y2)': return R(5, 3, 3, 'c');
        case '5f_y(3x2-y2)': return R(5, 3, 3, 's');
        // 5g 轨道
        case '5g_z4': return R(5, 4, 0, 'c');
        case '5g_z3x': return R(5, 4, 1, 'c');
        case '5g_z3y': return R(5, 4, 1, 's');
        case '5g_z2xy': return R(5, 4, 2, 's');
        case '5g_z2(x2-y2)': return R(5, 4, 2, 'c');
        case '5g_zx(x2-3y2)': return R(5, 4, 3, 'c');
        case '5g_zy(3x2-y2)': return R(5, 4, 3, 's');
        case '5g_xy(x2-y2)': return R(5, 4, 4, 's');
        case '5g_x4-6x2y2+y4': return R(5, 4, 4, 'c');
        default: return R(1, 0, 0, 'c');
    }
}

// ==================== 重要性采样函数 ====================

// ==================== 精确逆CDF采样 ====================

// 径向概率密度函数 P(r) = r² |R_nl(r)|²
function radialPDF(n, l, r, Z = 1, a0 = A0, atomType = 'H') {
    const R = radialR(n, l, r, Z, a0, atomType);
    return r * r * R * R;
}

// CDF 表缓存
const _cdfCache = {};

/**
 * 构建径向分布的累积分布函数表
 */
function buildRadialCDF(n, l, numPoints = 2000, atomType = 'H') {
    const key = `${n}_${l}_${atomType}`;
    if (_cdfCache[key]) {
        return _cdfCache[key];
    }

    const rMax = Math.max(4 * n * n * A0, 50 * A0);
    const dr = rMax / numPoints;

    const r = new Float64Array(numPoints + 1);
    const cdf = new Float64Array(numPoints + 1);

    r[0] = 0;
    cdf[0] = 0;
    let cumulative = 0;

    for (let i = 1; i <= numPoints; i++) {
        r[i] = i * dr;
        const P_prev = radialPDF(n, l, r[i - 1], 1, A0, atomType);
        const P_curr = radialPDF(n, l, r[i], 1, A0, atomType);
        cumulative += (P_prev + P_curr) * dr / 2;
        cdf[i] = cumulative;
    }

    const totalProb = cdf[numPoints];
    if (totalProb > 0) {
        for (let i = 0; i <= numPoints; i++) {
            cdf[i] /= totalProb;
        }
    }

    const result = { r, cdf, rMax, dr, numPoints, totalProb };
    _cdfCache[key] = result;
    return result;
}

/**
 * 精确逆CDF采样：从径向分布 P(r) 精确采样
 */
function sampleRadialExact(n, l, atomType = 'H') {
    const { r, cdf, numPoints } = buildRadialCDF(n, l, 2000, atomType);

    const u = Math.random();

    // 二分查找
    let lo = 0, hi = numPoints;
    while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (cdf[mid] < u) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    // 线性插值
    const i = Math.max(0, lo - 1);
    const j = Math.min(numPoints, lo);

    if (i === j || cdf[j] === cdf[i]) {
        return r[i];
    }

    const t = (u - cdf[i]) / (cdf[j] - cdf[i]);
    return r[i] + t * (r[j] - r[i]);
}

// ==================== 保留旧的重要性采样代码（兼容性）====================

/**
 * 计算轨道的所有径向峰位置（数值方法）
 */
function findRadialPeaks(n, l, a0 = A0) {
    const peaks = [];
    const rmax = n * n * 3 * a0;
    const dr = 0.1 * a0;

    let prevVal = 0;
    let prevPrevVal = 0;

    for (let r = dr; r < rmax; r += dr) {
        const R = radialR(n, l, r);
        const val = r * r * R * R;

        if (r > 2 * dr && prevVal > prevPrevVal && prevVal > val) {
            peaks.push({ r: r - dr, pdf: prevVal });
        }

        prevPrevVal = prevVal;
        prevVal = val;
    }

    return peaks;
}

// 缓存峰位置
const _peakCache = {};
function getCachedPeaks(n, l) {
    const key = `${n}_${l}`;
    if (!_peakCache[key]) {
        _peakCache[key] = findRadialPeaks(n, l);
    }
    return _peakCache[key];
}

/**
 * 多峰混合提议分布参数
 * 
 * 【物理优化】使用积分贡献估计权重，而非峰高度
 * 这样可以确保每个峰按其实际概率贡献获得采样机会
 */
function getMultiPeakMixtureParams(n, l, a0 = A0) {
    const peaks = getCachedPeaks(n, l);

    if (peaks.length === 0) {
        const r_peak = n * n * a0;
        return { components: [{ alpha: 2.0 / r_peak, weight: 1.0 }] };
    }

    // 【改进】计算每个峰的权重（基于估计的积分贡献）
    // 积分贡献 ≈ 峰高 × 峰宽，峰宽与 r_peak 成正比
    let totalWeight = 0;
    const components = [];

    for (const peak of peaks) {
        // 使用 pdf × r_peak 作为积分贡献的估计（更接近真实积分）
        const estimatedIntegral = peak.pdf * peak.r;
        const weight = estimatedIntegral;
        totalWeight += weight;

        const alpha = 2.0 / peak.r;
        components.push({ alpha, weight, r_peak: peak.r, pdf: peak.pdf });
    }

    for (const comp of components) {
        comp.weight /= totalWeight;
    }

    return { components };
}

/**
 * 计算多峰混合提议分布的 PDF
 */
function radialProposalPDF(r, n, l) {
    if (r <= 0) return 0;

    const params = getMultiPeakMixtureParams(n, l);
    let pdf = 0;

    for (const comp of params.components) {
        const alpha = comp.alpha;
        const norm = alpha * alpha * alpha / 2.0;
        pdf += comp.weight * norm * r * r * Math.exp(-alpha * r);
    }

    return pdf;
}

/**
 * 从 Gamma(3) 分布采样
 */
function sampleGamma3(alpha) {
    let sum = 0;
    for (let i = 0; i < 3; i++) {
        sum -= Math.log(Math.random()) / alpha;
    }
    return sum;
}

/**
 * 从多峰混合提议分布采样
 */
function sampleRadialProposal(n, l) {
    const params = getMultiPeakMixtureParams(n, l);

    const u = Math.random();
    let cumWeight = 0;

    for (const comp of params.components) {
        cumWeight += comp.weight;
        if (u < cumWeight) {
            return sampleGamma3(comp.alpha);
        }
    }

    return sampleGamma3(params.components[params.components.length - 1].alpha);
}

/**
 * 从均匀球面分布采样角度
 */
function sampleUniformSphere() {
    const phi = 2 * PI * Math.random();
    const cosTheta = 2 * Math.random() - 1;
    const theta = Math.acos(cosTheta);
    return { theta, phi };
}

// 缓存权重上界计算结果
const _maxWeightCache = {};

/**
 * 计算径向权重的理论上界
 * 【物理优化】使用数值搜索找到精确的最大权重
 */
function getMaxRadialWeight(n, l) {
    const key = `${n}_${l}`;
    if (_maxWeightCache[key]) {
        return _maxWeightCache[key];
    }

    const peaks = getCachedPeaks(n, l);
    let maxWeight = 1.0;

    // 在每个峰附近搜索最大权重
    for (const peak of peaks) {
        const rPeak = peak.r;
        const rMin = Math.max(0.1, rPeak * 0.3);
        const rMax = rPeak * 2.5;
        const numSamples = 200;
        const dr = (rMax - rMin) / numSamples;

        for (let i = 0; i <= numSamples; i++) {
            const r = rMin + i * dr;
            const R = radialR(n, l, r);
            const P_radial = r * r * R * R;
            const q_r = radialProposalPDF(r, n, l);

            if (q_r > 1e-300) {
                const w = P_radial / q_r;
                if (w > maxWeight) {
                    maxWeight = w;
                }
            }
        }
    }

    // 添加 20% 安全边际
    const result = maxWeight * 1.2;
    _maxWeightCache[key] = result;

    return result;
}

/**
 * 重要性采样：生成一个采样点
 * 
 * 【物理准确性保证】
 * 采用"分离变量"策略：
 * - 径向：使用精确逆CDF采样（100%接受率，无偏差）
 * - 角向：均匀球面采样 + |Y|² 权重接受-拒绝
 */
function importanceSample(n, l, angKey, samplingBoundary, atomType = 'H') {
    // ==================== 第一步：径向采样（精确逆CDF） ====================
    let r = sampleRadialExact(n, l, atomType);

    if (r > samplingBoundary * 2) {
        return null;
    }

    // ==================== 第二步：角向采样 ====================
    const { theta, phi } = sampleUniformSphere();

    const Y2 = realYlm_abs2(angKey.l, angKey.m, angKey.t, theta, phi);
    const w_angular = 4 * PI * Y2;
    // 【物理修正】角向权重上界的精确计算
    const maxAngularWeight = (angKey.m === 0) ? (2 * angKey.l + 1.2) : (2 * (2 * angKey.l + 1) + 0.5);

    if (Math.random() * maxAngularWeight > w_angular) {
        return { x: 0, y: 0, z: 0, r, theta, phi, weight: w_angular, accepted: false };
    }

    const sinTheta = Math.sin(theta);
    const x = r * sinTheta * Math.cos(phi);
    const y = r * sinTheta * Math.sin(phi);
    const z = r * Math.cos(theta);

    // 摩尔纹抑制
    const dither = 0.01;
    const dx = (Math.random() - 0.5) * dither;
    const dy = (Math.random() - 0.5) * dither;
    const dz = (Math.random() - 0.5) * dither;

    return { x: x + dx, y: y + dy, z: z + dz, r, theta, phi, weight: 1, accepted: true };
}

// ==================== 自适应网格采样模块 ====================
// 该模块提供通用的自适应采样能力，可处理任意形状的概率分布
// 包括标准原子轨道、杂化轨道、以及任意线性组合

const ADAPTIVE_GRID_SIZE = 24;  // 24³ = 13824 个网格单元，平衡精度与内存
const WARMUP_SAMPLES = 800;     // 预热阶段采样点数
const EXPLORATION_RATIO = 0.15; // 15% 用于探索新区域

// 自适应采样状态（每次新采样任务时重置）
let adaptiveGrid = {
    cells: null,           // Float32Array，存储每个网格单元的累积概率
    cellCounts: null,      // Uint32Array，每个单元被采中的次数
    totalHits: 0,          // 总命中次数
    boundary: 0,           // 当前网格覆盖的空间边界
    isWarmedUp: false,     // 是否已完成预热
    densityFunc: null      // 当前使用的密度函数
};

/**
 * 初始化自适应网格
 * @param {number} boundary - 采样空间边界（半径）
 */
function initAdaptiveGrid(boundary) {
    const totalCells = ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE;
    adaptiveGrid.cells = new Float32Array(totalCells);
    adaptiveGrid.cellCounts = new Uint32Array(totalCells);
    adaptiveGrid.totalHits = 0;
    adaptiveGrid.boundary = boundary;
    adaptiveGrid.isWarmedUp = false;

    // 初始化为均匀分布（每个格子概率相等）
    const uniformProb = 1.0 / totalCells;
    for (let i = 0; i < totalCells; i++) {
        adaptiveGrid.cells[i] = uniformProb;
    }
}

/**
 * 将空间坐标转换为网格索引
 */
function coordToGridIndex(x, y, z, boundary) {
    const halfSize = boundary;
    const cellSize = (2 * boundary) / ADAPTIVE_GRID_SIZE;

    let ix = Math.floor((x + halfSize) / cellSize);
    let iy = Math.floor((y + halfSize) / cellSize);
    let iz = Math.floor((z + halfSize) / cellSize);

    // 边界裁剪
    ix = Math.max(0, Math.min(ADAPTIVE_GRID_SIZE - 1, ix));
    iy = Math.max(0, Math.min(ADAPTIVE_GRID_SIZE - 1, iy));
    iz = Math.max(0, Math.min(ADAPTIVE_GRID_SIZE - 1, iz));

    return ix + iy * ADAPTIVE_GRID_SIZE + iz * ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE;
}

/**
 * 将网格索引转换为单元格中心坐标
 */
function gridIndexToCoord(index, boundary) {
    const cellSize = (2 * boundary) / ADAPTIVE_GRID_SIZE;
    const halfSize = boundary;

    const iz = Math.floor(index / (ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE));
    const iy = Math.floor((index % (ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE)) / ADAPTIVE_GRID_SIZE);
    const ix = index % ADAPTIVE_GRID_SIZE;

    return {
        x: (ix + 0.5) * cellSize - halfSize,
        y: (iy + 0.5) * cellSize - halfSize,
        z: (iz + 0.5) * cellSize - halfSize,
        cellSize: cellSize
    };
}

/**
 * 更新网格概率分布（基于新采样到的点）
 */
function updateAdaptiveGrid(x, y, z, density) {
    const index = coordToGridIndex(x, y, z, adaptiveGrid.boundary);

    // 累加密度值
    adaptiveGrid.cells[index] += density;
    adaptiveGrid.cellCounts[index]++;
    adaptiveGrid.totalHits++;
}

/**
 * 归一化网格概率（预热结束后调用）
 */
function normalizeAdaptiveGrid() {
    const totalCells = ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE;
    let totalProb = 0;

    for (let i = 0; i < totalCells; i++) {
        totalProb += adaptiveGrid.cells[i];
    }

    if (totalProb > 0) {
        for (let i = 0; i < totalCells; i++) {
            adaptiveGrid.cells[i] /= totalProb;
        }
    }

    adaptiveGrid.isWarmedUp = true;
}

/**
 * 从自适应网格中采样一个单元格
 * 使用逆变换采样（Inverse Transform Sampling）
 */
function sampleFromAdaptiveGrid() {
    const totalCells = ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE * ADAPTIVE_GRID_SIZE;
    const u = Math.random();

    let cumulative = 0;
    for (let i = 0; i < totalCells; i++) {
        cumulative += adaptiveGrid.cells[i];
        if (u <= cumulative) {
            return i;
        }
    }

    return totalCells - 1; // 边界情况
}

/**
 * 通用密度函数接口 - 计算任意轨道组合在某点的概率密度
 * @param {Array} orbitalConfigs - 轨道配置数组 [{params, coefficient}, ...]
 * @param {number} r - 径向距离
 * @param {number} theta - 极角
 * @param {number} phi - 方位角
 * @returns {number} - 概率密度 |Ψ|²
 */
function computeGeneralDensity(orbitalConfigs, r, theta, phi) {
    if (!orbitalConfigs || orbitalConfigs.length === 0) return 0;

    // 如果只有一个轨道且系数为1，直接使用原有公式（最优路径）
    if (orbitalConfigs.length === 1 && orbitalConfigs[0].coefficient === 1) {
        const p = orbitalConfigs[0].params;
        return density3D_real(p.angKey, p.n, p.l, r, theta, phi);
    }

    // 计算波函数的线性叠加：Ψ = Σ c_i * Ψ_i
    let psiReal = 0;
    let psiImag = 0; // 预留，目前只用实球谐函数

    for (const config of orbitalConfigs) {
        const { params, coefficient } = config;
        const R = radialR(params.n, params.l, r);
        const Y = realYlm_value(params.angKey.l, params.angKey.m, params.angKey.t, theta, phi);
        psiReal += coefficient * R * Y;
    }

    // 返回概率密度 |Ψ|²
    return psiReal * psiReal + psiImag * psiImag;
}

/**
 * 自适应采样 - 主入口函数
 * 适用于任意形状的概率分布，包括杂化轨道
 * @param {Array} orbitalConfigs - 轨道配置
 * @param {number} samplingBoundary - 采样边界
 * @returns {Object|null} - 采样结果
 */
function adaptiveSample(orbitalConfigs, samplingBoundary) {
    const boundary = samplingBoundary;
    const cellSize = (2 * boundary) / ADAPTIVE_GRID_SIZE;

    let x, y, z;

    // 决定是探索还是利用
    if (!adaptiveGrid.isWarmedUp || Math.random() < EXPLORATION_RATIO) {
        // 探索模式：均匀随机采样
        x = (Math.random() - 0.5) * 2 * boundary;
        y = (Math.random() - 0.5) * 2 * boundary;
        z = (Math.random() - 0.5) * 2 * boundary;
    } else {
        // 利用模式：根据网格概率采样
        const cellIndex = sampleFromAdaptiveGrid();
        const center = gridIndexToCoord(cellIndex, boundary);

        // 在选中的单元格内均匀采样
        x = center.x + (Math.random() - 0.5) * cellSize;
        y = center.y + (Math.random() - 0.5) * cellSize;
        z = center.z + (Math.random() - 0.5) * cellSize;
    }

    // 计算球坐标
    const r = Math.sqrt(x * x + y * y + z * z);
    if (r < 1e-10) return null;

    const theta = Math.acos(z / r);
    const phi = Math.atan2(y, x);

    // 计算该点的概率密度
    const density = computeGeneralDensity(orbitalConfigs, r, theta, phi);

    // 更新网格（无论是否接受，都要记录密度信息）
    updateAdaptiveGrid(x, y, z, density);

    // 接受-拒绝判定
    // 使用动态缩放因子，基于当前已知的最大密度
    const maxDensity = Math.max(1e-10, density * 2); // 保守估计
    const acceptProb = density / maxDensity;

    if (Math.random() < acceptProb) {
        // 摩尔纹抑制
        const dither = 0.01;
        return {
            x: x + (Math.random() - 0.5) * dither,
            y: y + (Math.random() - 0.5) * dither,
            z: z + (Math.random() - 0.5) * dither,
            r, theta, phi,
            weight: 1,
            accepted: true
        };
    }

    return { x: 0, y: 0, z: 0, r, theta, phi, weight: 0, accepted: false };
}

/**
 * 判断是否为标准氢原子轨道（可使用高效重要性采样）
 * @param {Object} config - 采样配置
 * @returns {boolean}
 */
function isStandardHydrogenMode(config) {
    // 条件1：必须有轨道配置
    if (!config.orbitals || config.orbitals.length === 0) return false;

    // 条件2：不是杂化模式（杂化模式使用拒绝采样）
    if (config.isHybridMode) return false;

    // 条件3：不是线性组合模式
    if (config.isLinearCombination) return false;

    // 条件4：所有轨道都是已知的标准氢原子轨道
    for (const key of config.orbitals) {
        const params = orbitalParamsFromKey(key);
        if (!params) return false;
    }

    return true;
}

// ==================== 采样逻辑 ====================

// 是否启用重要性采样
let useImportanceSampling = true;

// 缓存杂化轨道的最大密度估计值
let hybridMaxDensityCache = null;

/**
 * 批量采样函数 - 在 Worker 中执行
 * 
 * 【重要修复】多选/比照模式使用轮流采样（Round-Robin）策略：
 * - 为每个轨道平均分配采样名额
 * - 每个轨道使用自己的重要性采样提议分布
 * - 避免了随机选择轨道导致的密度偏差
 * 
 * @param {Object} config - 采样配置
 * @returns {Object} - 采样结果
 */
function performSampling(config) {
    const {
        orbitals,           // 轨道键数组
        samplingBoundary,   // 采样边界
        maxAttempts,        // 最大尝试次数
        targetPoints,       // 目标点数
        isIndependentMode,  // 是否为独立模式（多选/比照）
        isMultiselectMode,  // 是否为多选模式
        isCompareMode,      // 【新增】是否为比照模式
        isHybridMode,       // 是否为杂化模式
        phaseOn,            // 是否显示相位
        compareColors,      // 比照模式颜色
        atomType,           // 原子类型
        slotConfigs         // 【重构】比照模式slot配置数组 [{orbital, atom, slotIndex}]
    } = config;

    // 【重构】比照模式使用slotConfigs构建paramsList
    let paramsList;
    let slotAtomTypes = [];  // 每个slot对应的原子类型
    let effectiveOrbitals; // 【关键修复】与paramsList一一对应的轨道键数组

    if (isCompareMode && slotConfigs && slotConfigs.length > 0) {
        // 比照模式：每个slot独立
        paramsList = slotConfigs.map(s => orbitalParamsFromKey(s.orbital)).filter(Boolean);
        slotAtomTypes = slotConfigs.map(s => s.atom || 'H');
        // 【关键修复】使用slotConfigs中的orbital作为轨道键，确保长度与paramsList一致
        effectiveOrbitals = slotConfigs.map(s => s.orbital);
    } else {
        // 其他模式：使用orbitals数组
        paramsList = orbitals.map(k => orbitalParamsFromKey(k)).filter(Boolean);
        slotAtomTypes = paramsList.map(() => atomType || 'H');
        effectiveOrbitals = orbitals;
    }

    if (!paramsList.length) {
        return { points: [], samples: [], farthestDistance: 0 };
    }

    // 【杂化模式】使用专门的杂化采样逻辑
    if (isHybridMode) {
        return performHybridSampling(paramsList, samplingBoundary, maxAttempts, targetPoints, phaseOn, atomType);
    }

    const samplingVolumeSize = samplingBoundary * 2;
    const points = [];       // 采样到的点 [{x,y,z,r,g,b,orbitalIndex}]
    const samples = [];      // 用于图表的样本数据 [{r,theta,orbitalKey}]
    let farthestDistance = 0;
    let attempts = 0;

    // 【关键修复】多轨道模式下，确保每个轨道获得相同数量的最终采样点
    const numOrbitals = paramsList.length;
    const orbitalPointCounts = new Array(numOrbitals).fill(0); // 每个轨道已采样的点数
    const targetPointsPerOrbital = Math.ceil(targetPoints / numOrbitals); // 每个轨道的目标点数

    // 找到当前点数最少的轨道进行采样
    function getOrbitalWithFewestPoints() {
        let minIdx = 0;
        let minCount = orbitalPointCounts[0];
        for (let i = 1; i < numOrbitals; i++) {
            if (orbitalPointCounts[i] < minCount) {
                minCount = orbitalPointCounts[i];
                minIdx = i;
            }
        }
        return minIdx;
    }

    while (attempts < maxAttempts && points.length < targetPoints) {
        attempts++;

        let accepted = false;
        let point = null;
        let x, y, z, r, theta, phi;
        let orbitalIndex = 0;

        // 【关键修复】选择当前点数最少的轨道进行采样，确保各轨道点数均衡
        if (isIndependentMode && numOrbitals > 1) {
            orbitalIndex = getOrbitalWithFewestPoints();
        }
        const p = paramsList[orbitalIndex];

        // 【重构】使用slotAtomTypes数组获取该slot对应的原子类型
        const currentAtomType = slotAtomTypes[orbitalIndex] || atomType || 'H';

        // 使用重要性采样（针对当前轨道优化）
        if (useImportanceSampling) {
            const result = importanceSample(p.n, p.l, p.angKey, samplingBoundary, currentAtomType);

            if (result && result.accepted) {
                x = result.x;
                y = result.y;
                z = result.z;
                r = result.r;
                theta = result.theta;
                phi = result.phi;
                accepted = true;
            }
        } else {
            // 传统拒绝采样（降级方案）
            x = (Math.random() - 0.5) * samplingVolumeSize;
            y = (Math.random() - 0.5) * samplingVolumeSize;
            z = (Math.random() - 0.5) * samplingVolumeSize;
            r = Math.sqrt(x * x + y * y + z * z);

            if (r === 0) continue;

            theta = Math.acos(z / r);
            phi = Math.atan2(y, x);

            if (isIndependentMode) {
                const density = density3D_real(p.angKey, p.n, p.l, r, theta, phi, 1, 1, currentAtomType);
                const scaleFactor = Math.min(200, 3.0 * Math.pow(p.n, 3));
                accepted = Math.random() < density * scaleFactor;
            } else {
                let probability = 0;
                let minN = paramsList[0].n;
                for (let i = 0; i < paramsList.length; i++) {
                    const params = paramsList[i];
                    probability += density3D_real(params.angKey, params.n, params.l, r, theta, phi, 1, 1, atomType);
                    if (params.n < minN) minN = params.n;
                }
                const scaleFactor = Math.min(200, 3.0 * Math.pow(minN, 3) / Math.sqrt(paramsList.length));
                accepted = Math.random() < probability * scaleFactor;
            }
        }

        if (accepted) {
            let r_color, g_color, b_color;

            if (isIndependentMode && !isMultiselectMode) {
                // 比照模式：固定颜色
                if (compareColors && orbitalIndex < compareColors.length) {
                    const color = compareColors[orbitalIndex];
                    r_color = color[0];
                    g_color = color[1];
                    b_color = color[2];
                } else {
                    r_color = 1; g_color = 1; b_color = 1;
                }
            } else {
                // 单选模式或多选模式：支持相位
                if (phaseOn) {
                    let psi = 0;
                    if (isMultiselectMode) {
                        // 多选模式：只计算当前轨道的相位
                        const R = radialR(p.n, p.l, r, 1, 1, atomType);
                        const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                        psi = R * Y;
                    } else {
                        // 单选模式：计算所有轨道的叠加相位
                        for (const params of paramsList) {
                            const R = radialR(params.n, params.l, r, 1, 1, atomType);
                            const Y = realYlm_value(params.angKey.l, params.angKey.m, params.angKey.t, theta, phi);
                            psi += R * Y;
                        }
                    }
                    // 相位颜色：与主线程保持一致（正=蓝、负=红）
                    if (psi > 0) {
                        r_color = 0; g_color = 0; b_color = 1;
                    } else if (psi < 0) {
                        r_color = 1; g_color = 0; b_color = 0;
                    } else {
                        r_color = 1; g_color = 1; b_color = 1;
                    }
                } else {
                    r_color = 1; g_color = 1; b_color = 1;
                }
            }

            point = { x, y, z, r: r_color, g: g_color, b: b_color, orbitalIndex: isIndependentMode ? orbitalIndex : -1 };
            // 【关键修复】添加orbitalIndex到samples，用于比照模式构建唯一键
            // 【关键修复】使用effectiveOrbitals确保比照模式下每个slot都有正确的轨道键
            // 【新增】添加phi用于φ角向分布图表
            samples.push({ r, theta, phi, orbitalKey: isIndependentMode ? effectiveOrbitals[orbitalIndex] : null, orbitalIndex: isIndependentMode ? orbitalIndex : -1 });
            points.push(point);

            // 【关键】更新该轨道的已采样点数计数
            if (isIndependentMode && numOrbitals > 1) {
                orbitalPointCounts[orbitalIndex]++;
            }

            if (r > farthestDistance) {
                farthestDistance = r;
            }
        }
    }

    return {
        points,
        samples,
        farthestDistance,
        attempts
    };
}

// ==================== 杂化模式采样 ====================

// CDF 缓存（用于杂化轨道高效采样）
// CDF 缓存（用于杂化轨道高效采样）
const _hybridCdfCache = {};

/**
 * 构建杂化轨道的径向CDF
 */
function buildHybridRadialCDF(paramsList, numPoints = 2000, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return null;

    const key = paramsList.map(p => `${p.n}_${p.l}`).join('|') + `_${atomType}`;
    if (_hybridCdfCache[key]) {
        return _hybridCdfCache[key];
    }

    const numOrbitals = paramsList.length;
    const defaultCoeff = 1.0 / Math.sqrt(numOrbitals);

    // 确定最大积分范围
    let rMax = 0;
    for (const p of paramsList) {
        rMax = Math.max(rMax, 4 * p.n * p.n * A0);
    }
    rMax = Math.max(rMax, 50 * A0);

    const dr = rMax / numPoints;
    const r = new Float64Array(numPoints + 1);
    const cdf = new Float64Array(numPoints + 1);

    // 数值积分
    r[0] = 0;
    cdf[0] = 0;
    let cumulative = 0;

    for (let i = 1; i <= numPoints; i++) {
        r[i] = i * dr;

        let P_prev = 0, P_curr = 0;
        const r_prev = r[i - 1];
        const r_curr = r[i];

        // 【数学严谨性修复】为支持非正交轨道组合（如 1s+2s），引入双重循环计算交叉项
        // P(r) = r² * Σ_i Σ_j c_i c_j R_i R_j * ⟨Y_i|Y_j⟩

        let densityPrev = 0;
        let densityCurr = 0;

        for (let j = 0; j < numOrbitals; j++) {
            const p1 = paramsList[j];
            const c1 = p1.coefficient !== undefined ? p1.coefficient : defaultCoeff;
            const R1_prev = radialR(p1.n, p1.l, r_prev, 1, A0, atomType);
            const R1_curr = radialR(p1.n, p1.l, r_curr, 1, A0, atomType);

            for (let k = 0; k < numOrbitals; k++) {
                const p2 = paramsList[k];
                const c2 = p2.coefficient !== undefined ? p2.coefficient : defaultCoeff;

                // 检查角向正交性
                const theta1 = p1.angKey ? p1.angKey.t : '';
                const theta2 = p2.angKey ? p2.angKey.t : '';
                const isNonOrthogonal = (p1.angKey.l === p2.angKey.l) &&
                    (p1.angKey.m === p2.angKey.m) &&
                    (theta1 === theta2);

                if (isNonOrthogonal) {
                    const R2_prev = radialR(p2.n, p2.l, r_prev, 1, A0, atomType);
                    const R2_curr = radialR(p2.n, p2.l, r_curr, 1, A0, atomType);

                    densityPrev += c1 * c2 * R1_prev * R2_prev;
                    densityCurr += c1 * c2 * R1_curr * R2_curr;
                }
            }
        }

        P_prev += r_prev * r_prev * densityPrev;
        P_curr += r_curr * r_curr * densityCurr;

        cumulative += (P_prev + P_curr) * dr / 2;
        cdf[i] = cumulative;
    }

    // 归一化
    const totalProb = cdf[numPoints];
    if (totalProb > 0) {
        for (let i = 0; i <= numPoints; i++) {
            cdf[i] /= totalProb;
        }
    }

    const result = { r, cdf, rMax, dr, numPoints, totalProb };
    _hybridCdfCache[key] = result;
    return result;
}

/**
 * 从均匀球面分布采样角度
 */
function sampleUniformSphere() {
    const phi = 2 * PI * Math.random();
    const cosTheta = 2 * Math.random() - 1;
    const theta = Math.acos(cosTheta);
    return { theta, phi };
}

/**
 * 高效杂化轨道采样（基于精确径向CDF + 角向接受-拒绝）
 */
function hybridPreciseSample(paramsList, samplingBoundary, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return null;

    const numOrbitals = paramsList.length;
    const defaultCoeff = 1.0 / Math.sqrt(numOrbitals);

    // 构建/获取杂化CDF
    const cdfData = buildHybridRadialCDF(paramsList, 2000, atomType);
    if (!cdfData) return null;

    // 第一步：从杂化径向CDF精确采样
    const { r: rTable, cdf, numPoints } = cdfData;
    const u = Math.random();

    // 二分查找
    let lo = 0, hi = numPoints;
    while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (cdf[mid] < u) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    const i = Math.max(0, lo - 1);
    const j = Math.min(numPoints, lo);
    let r;
    if (i === j || cdf[j] === cdf[i]) {
        r = rTable[i];
    } else {
        const t = (u - cdf[i]) / (cdf[j] - cdf[i]);
        r = rTable[i] + t * (rTable[j] - rTable[i]);
    }

    // 边界检查
    if (r > samplingBoundary * 2 || r < 1e-10) {
        return { x: 0, y: 0, z: 0, r, theta: 0, phi: 0, psi: 0, accepted: false };
    }

    // 第二步：均匀球面采样角度
    const { theta, phi } = sampleUniformSphere();

    // 第三步：计算角向权重
    let psi = 0;
    let radialSum2 = 0;

    for (const p of paramsList) {
        const coeff = p.coefficient !== undefined ? p.coefficient : defaultCoeff;
        const R = radialR(p.n, p.l, r, 1, A0, atomType);
        const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
        psi += coeff * R * Y;
        radialSum2 += coeff * coeff * R * R;
    }

    const angularFactor = psi * psi;
    const expectedFactor = radialSum2 / (4 * PI);

    if (expectedFactor < 1e-300) {
        return { x: 0, y: 0, z: 0, r, theta, phi, psi, accepted: false };
    }

    const angularWeight = angularFactor / expectedFactor;
    // 【关键修复】角向权重上界应与轨道数量成正比
    // 当多个轨道的球谐函数相干叠加时，|Σ c_i Y_i|² 可能远大于单轨道情况
    // 旧值 8π 导致高密度区域样本被错误拒绝
    const maxAngularWeight = 4 * PI * (numOrbitals + 2);

    // 角向接受-拒绝
    if (Math.random() * maxAngularWeight > angularWeight) {
        return { x: 0, y: 0, z: 0, r, theta, phi, psi, accepted: false };
    }

    // 采样成功
    const sinTheta = Math.sin(theta);
    const x = r * sinTheta * Math.cos(phi);
    const y = r * sinTheta * Math.sin(phi);
    const z = r * Math.cos(theta);

    const dither = 0.01;
    return {
        x: x + (Math.random() - 0.5) * dither,
        y: y + (Math.random() - 0.5) * dither,
        z: z + (Math.random() - 0.5) * dither,
        r, theta, phi,
        psi,
        accepted: true
    };
}

/**
 * 杂化模式采样函数
 * 
 * 【物理原理】
 * 杂化轨道是多个原子轨道波函数的线性组合：Ψ_hybrid = Σ c_i × Ψ_i
 * 概率密度为 |Ψ_hybrid|² = |Σ c_i × R_i(r) × Y_i(θ,φ)|²
 * 
 * 【采样方法】
 * 优先使用高效的精确CDF方法，回退到基础拒绝采样
 * 
 * @param {Array} paramsList - 轨道参数列表
 * @param {number} samplingBoundary - 采样边界
 * @param {number} maxAttempts - 最大尝试次数
 * @param {number} targetPoints - 目标点数
 * @param {boolean} phaseOn - 是否显示相位
 * @returns {Object} - 采样结果
 */
function performHybridSampling(paramsList, samplingBoundary, maxAttempts, targetPoints, phaseOn, atomType = 'H') {
    const points = [];
    const samples = [];
    let farthestDistance = 0;
    let attempts = 0;

    // 【关键修复】使用 Thomson + SVD 算法计算杂化系数
    const hybridResult = getHybridCoefficientsWorker(paramsList);
    if (!hybridResult) {
        console.error("Worker: Failed to compute hybrid coefficients");
        return { points, samples, farthestDistance, attempts };
    }

    const { coeffMatrix, sortedParams } = hybridResult;
    const numOrbitals = sortedParams.length;

    // 轮流采样每条杂化轨道（Round-Robin）
    let hybridIndex = 0;

    // 尝试使用高效采样方法
    while (attempts < maxAttempts && points.length < targetPoints) {
        attempts++;

        // 选择当前杂化轨道的系数
        const k = hybridIndex % numOrbitals;
        hybridIndex++;
        const coeffs = coeffMatrix[k];

        // 将系数注入到排序后的参数列表
        const paramsWithCoeffs = sortedParams.map((p, i) => ({
            ...p,
            coefficient: coeffs[i]
        }));

        const result = hybridPreciseSample(paramsWithCoeffs, samplingBoundary, atomType);

        if (!result || !result.accepted) continue;

        // 计算颜色
        let r_color, g_color, b_color;
        if (phaseOn) {
            // 相位颜色：与主线程保持一致（正=蓝、负=红）
            if (result.psi > 0) {
                r_color = 0; g_color = 0; b_color = 1;
            } else if (result.psi < 0) {
                r_color = 1; g_color = 0; b_color = 0;
            } else {
                r_color = 1; g_color = 1; b_color = 1;
            }
        } else {
            r_color = 1; g_color = 1; b_color = 1;
        }

        points.push({
            x: result.x,
            y: result.y,
            z: result.z,
            r: r_color,
            g: g_color,
            b: b_color,
            orbitalIndex: -1
        });

        samples.push({ r: result.r, theta: result.theta, orbitalKey: null });

        if (result.r > farthestDistance) {
            farthestDistance = result.r;
        }
    }

    return {
        points,
        samples,
        farthestDistance,
        attempts
    };
}

// ==================== Worker 消息处理 ====================

self.onmessage = function (e) {
    const { type, data, taskId, sessionId } = e.data;

    switch (type) {
        case 'sample':
            // 执行采样任务
            const result = performSampling(data);
            self.postMessage({
                type: 'sample-result',
                taskId: taskId,
                sessionId: sessionId, // 返回会话 ID 以便主线程验证
                result: result
            });
            break;

        case 'reset-hybrid-cache':
            // 重置杂化缓存
            hybridMaxDensityCache = null;
            self.postMessage({ type: 'hybrid-cache-reset', taskId: taskId });
            break;

        case 'ping':
            // 健康检查
            self.postMessage({ type: 'pong', taskId: taskId });
            break;

        default:
            console.warn('Worker: Unknown message type', type);
    }
};

// 通知主线程 Worker 已就绪
self.postMessage({ type: 'ready' });
