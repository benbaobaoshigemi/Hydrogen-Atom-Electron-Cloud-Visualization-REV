/**
 * 共享物理核心模块 (v2.0)
 * 
 * 此文件包含 physics.js 和 sampling-worker.js 共享的核心物理计算函数。
 * 解决了“三重冗余”问题，成为项目唯一的物理逻辑来源。
 * 
 * 【重要】修改此文件将同时影响主线程和 Worker，请务必进行全量回归测试。
 */

(function (root) {
    'use strict';

    const A0 = 1; // Bohr radius unit
    const PI = Math.PI;
    const TWO_PI = 2 * Math.PI;

    // 缓存容器
    const _thomsonCache = {};
    const _cdfCache = {};
    const _peakCache = {};

    // 阶乘表（预计算）
    const FACT = (() => {
        const f = [1];
        for (let i = 1; i <= 64; i++) f[i] = f[i - 1] * i;
        return f;
    })();

    // ==================== 1. 基础数学与特殊函数 ====================

    function factorial(n) {
        if (n < 0) return 0;
        if (n >= FACT.length) return Infinity;
        return FACT[n];
    }

    function binomialInt(n, k) {
        if (k < 0 || k > n) return 0;
        return factorial(n) / (factorial(k) * factorial(n - k));
    }

    function generalizedLaguerre(k, alpha, x) {
        let sum = 0;
        for (let i = 0; i <= k; i++) {
            const c = ((i % 2) ? -1 : 1) * binomialInt(k + alpha, k - i) / factorial(i);
            sum += c * Math.pow(x, i);
        }
        return sum;
    }

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

    // ==================== 2. 线性代数工具 ====================

    function matMul(A, B) {
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

    function matTranspose(A) {
        const m = A.length, n = A[0].length;
        const AT = [];
        for (let i = 0; i < n; i++) {
            const row = new Array(m).fill(0);
            for (let j = 0; j < m; j++) row[j] = A[j][i];
            AT.push(row);
        }
        return AT;
    }

    function jacobiSVD(A) {
        const m = A.length, n = A[0].length;
        let U = A.map(row => [...row]);
        let V = [];
        for (let i = 0; i < n; i++) {
            const row = new Array(n).fill(0);
            row[i] = 1;
            V.push(row);
        }
        const EPSILON = 1e-15, MAX_ITER = 50;
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
                    let c = 1 / Math.sqrt(1 + t * t), s = c * t;
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

    // ==================== 3. 核心物理计算 (Ylm, Radial, STO) ====================

    function Ylm_complex(l, m, theta, phi) {
        const mm = Math.abs(m);
        const Plm = associatedLegendre(l, mm, Math.cos(theta));
        const N = Math.sqrt(((2 * l + 1) / (4 * PI)) * (factorial(l - mm) / factorial(l + mm)));
        const base = N * Plm;
        if (m === 0) return { re: base, im: 0 };
        const cos_mphi = Math.cos(mm * phi), sin_mphi = Math.sin(mm * phi);
        if (m > 0) return { re: base * cos_mphi, im: base * sin_mphi };
        const sign = (mm % 2) ? -1 : 1;
        return { re: sign * base * cos_mphi, im: -sign * base * sin_mphi };
    }

    function realYlm_value(l, m, type, theta, phi) {
        if (m === 0) return Ylm_complex(l, 0, theta, phi).re;
        const mm = Math.abs(m);
        const y = Ylm_complex(l, mm, theta, phi);
        const csPhaseCorrection = (mm % 2 === 1) ? -1 : 1;
        if (type === 'c' || type === 'cos') return csPhaseCorrection * Math.SQRT2 * y.re;
        return csPhaseCorrection * Math.SQRT2 * y.im;
    }

    function realYlm_abs2(l, m, type, theta, phi) {
        const val = realYlm_value(l, m, type, theta, phi);
        return val * val;
    }

    function getOrbitalKey(n, l) {
        const subshells = ['s', 'p', 'd', 'f', 'g'];
        return n + (subshells[l] || '');
    }

    function slaterRadialR(basis, r) {
        if (!basis) return 0;
        let sum = 0;
        for (let i = 0; i < basis.length; i++) {
            const { nStar, zeta, coeff } = basis[i];
            const norm = Math.pow(2 * zeta, nStar + 0.5) / Math.sqrt(factorial(2 * nStar));
            sum += coeff * norm * Math.pow(r, nStar - 1) * Math.exp(-zeta * r);
        }
        return sum;
    }

    function getSlaterBasis() {
        if (typeof self !== 'undefined' && self.SlaterBasis) return self.SlaterBasis;
        if (typeof window !== 'undefined' && window.SlaterBasis) return window.SlaterBasis;
        return null;
    }

    function radialR(n, l, r, Z = 1, a0 = A0, atomType = 'H') {
        const SlaterBasis = getSlaterBasis();
        if (atomType && atomType !== 'H' && SlaterBasis) {
            const atomData = SlaterBasis[atomType];
            const basis = atomData && atomData.orbitals && atomData.orbitals[getOrbitalKey(n, l)];
            if (basis) {
                const val = slaterRadialR(basis, r);
                const analyticalSign = ((n - l - 1) % 2 === 0) ? 1 : -1;
                let dom = basis[0];
                for (let i = 1; i < basis.length; i++) {
                    if (basis[i].zeta < dom.zeta || (Math.abs(basis[i].zeta - dom.zeta) < 1e-9 && basis[i].nStar > dom.nStar)) dom = basis[i];
                }
                if (Math.sign(dom.coeff) * analyticalSign < 0) return -val;
                return val;
            }
            // 【物理严谨性强制】对于非氢原子，如果没有 STO 基组数据，严禁退化到氢原子公式。
            return 0;
        }
        if (n <= 0 || l < 0 || l >= n) return 0;
        const rho = (2 * Z * r) / (n * a0);
        const k = n - l - 1;
        const pref = Math.pow(2 * Z / (n * a0), 1.5) * Math.sqrt(factorial(k) / (2 * n * factorial(n + l)));
        return pref * Math.exp(-rho / 2) * Math.pow(rho, l) * generalizedLaguerre(k, 2 * l + 1, rho);
    }

    function radialPDF(n, l, r, Z = 1, a0 = A0, atomType = 'H') {
        const R = radialR(n, l, r, Z, a0, atomType);
        return r * r * R * R;
    }

    // ==================== 4. 杂化轨道逻辑 (Thomson, SVD) ====================

    function optimizeThomson(n, maxIter = 200, lr = 0.1) {
        if (n <= 1) return [[0, 0, 1]];
        if (_thomsonCache[n]) return JSON.parse(JSON.stringify(_thomsonCache[n]));
        let points = [];
        const gr = (1 + Math.sqrt(5)) / 2;
        for (let i = 0; i < n; i++) {
            const y = 1 - (i / (n - 1)) * 2, rad = Math.sqrt(1 - y * y), t = TWO_PI * i / gr;
            points.push([Math.cos(t) * rad, y, Math.sin(t) * rad]);
        }
        for (let iter = 0; iter < maxIter; iter++) {
            const grad = points.map(() => [0, 0, 0]);
            for (let i = 0; i < n; i++) {
                for (let j = i + 1; j < n; j++) {
                    const dx = points[i][0] - points[j][0], dy = points[i][1] - points[j][1], dz = points[i][2] - points[j][2];
                    const d2 = dx * dx + dy * dy + dz * dz, d3 = d2 * Math.sqrt(d2) + 1e-10;
                    grad[i][0] -= dx / d3; grad[i][1] -= dy / d3; grad[i][2] -= dz / d3;
                    grad[j][0] += dx / d3; grad[j][1] += dy / d3; grad[j][2] += dz / d3;
                }
            }
            for (let i = 0; i < n; i++) {
                points[i][0] -= lr * grad[i][0]; points[i][1] -= lr * grad[i][1]; points[i][2] -= lr * grad[i][2];
                const norm = Math.sqrt(points[i][0] ** 2 + points[i][1] ** 2 + points[i][2] ** 2);
                points[i] = [points[i][0] / norm, points[i][1] / norm, points[i][2] / norm];
            }
            if (iter % 50 === 0) lr *= 0.8;
        }
        _thomsonCache[n] = JSON.parse(JSON.stringify(points));
        return points;
    }

    function generateConstrainedDirections(orbitalParams) {
        // 彻底移除所有硬编码的几何对齐逻辑 (如 sp3 -> 四面体)
        // 回退到纯粹的数值推演：解 Thomson 问题 (极外层电子斥能最小化)
        return optimizeThomson(orbitalParams.length);
    }

    function sortOrbitalsForHybridization(paramsList) {
        return [...paramsList].sort((a, b) => {
            if (a.l !== b.l) return a.l - b.l;
            if (a.n !== b.n) return a.n - b.n;
            const score = p => {
                if (p.l === 1) return (p.angKey.m === 1 && p.angKey.t === 'c') ? 1 : (p.angKey.m === 1 && p.angKey.t === 's' ? 2 : 3);
                if (p.l === 2) return (p.angKey.m === 2 && p.angKey.t === 'c') ? 1 : (p.angKey.m === 0 ? 2 : 3 + p.angKey.m);
                return 0;
            };
            return score(a) - score(b);
        });
    }

    function buildDirectionMatrix(directions, orbitalParams) {
        return directions.map(d => {
            const norm = Math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2);
            const x = d[0] / norm, y = d[1] / norm, z = d[2] / norm;

            // 转换为球面坐标
            const theta = Math.acos(Math.max(-1, Math.min(1, z)));
            const phi = Math.atan2(y, x);

            return orbitalParams.map(p => {
                const { l, m, t } = p.angKey;
                // 使用通用的实球谐函数，针对方向 (theta, phi) 计算波函数振幅
                return realYlm_value(l, m, t, theta, phi);
            });
        });
    }

    function getHybridCoefficients(orbitalParams) {
        const sorted = sortOrbitalsForHybridization(orbitalParams);
        const directions = generateConstrainedDirections(sorted);
        const A = buildDirectionMatrix(directions, sorted);
        const { U, V } = jacobiSVD(A);
        return matMul(U, matTranspose(V));
    }

    // ==================== 5. 高级分布计算 ====================
    // 计算所有杂化轨道叠加的总密度（对于 sp³ 为 4 条轨道叠加，结果应趋于各向同性）
    function allHybridOrbitalsDensity3D(paramsList, r, theta, phi, Z = 1, a0 = A0, atomType = 'H') {
        if (!paramsList || paramsList.length === 0) return 0;
        const numOrbitals = paramsList.length;
        const coeffMatrix = getHybridCoefficients(paramsList);
        let totalDensity = 0;
        for (let i = 0; i < numOrbitals; i++) {
            totalDensity += singleHybridDensity3D(paramsList, i, r, theta, phi, Z, a0, atomType);
        }
        return totalDensity;
    }

    // 计算特定一条杂化轨道的概率密度
    function singleHybridDensity3D(paramsList, hybridIndex, r, theta, phi, Z = 1, a0 = A0, atomType = 'H') {
        if (!paramsList || paramsList.length === 0) return 0;
        const coeffMatrix = getHybridCoefficients(paramsList);
        const idx = hybridIndex % paramsList.length;
        const coeffs = coeffMatrix[idx];
        let psi = 0;
        for (let i = 0; i < paramsList.length; i++) {
            const p = paramsList[i];
            psi += coeffs[i] * radialR(p.n, p.l, r, Z, a0, atomType) * realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
        }
        return psi * psi;
    }

    function singleHybridWavefunction(paramsList, hybridIndex, r, theta, phi, Z = 1, a0 = A0, atomType = 'H') {
        if (!paramsList || paramsList.length === 0) return 0;
        const coeffMatrix = getHybridCoefficients(paramsList);
        const idx = hybridIndex % paramsList.length;
        const coeffs = coeffMatrix[idx];
        let psi = 0;
        for (let i = 0; i < paramsList.length; i++) {
            const p = paramsList[i];
            psi += coeffs[i] * radialR(p.n, p.l, r, Z, a0, atomType) * realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
        }
        return psi;
    }

    function hybridEstimateMaxDensity(paramsList, rMax, numSamples = 1000) {
        if (!paramsList || paramsList.length === 0) return 0;
        let maxDensity = 0;
        for (let i = 0; i < numSamples; i++) {
            const r = Math.random() * rMax, ct = 2 * Math.random() - 1, t = Math.acos(ct), p = TWO_PI * Math.random();
            const d = singleHybridDensity3D(paramsList, 0, r, t, p); // 假设各瓣等效
            if (d > maxDensity) maxDensity = d;
        }
        return maxDensity * 1.5;
    }

    function orbitalParamsFromKey(key) {
        if (!key) return null;
        const R = (n, l, m, t) => ({ n, l, angKey: { l, m, t } });

        // Regex 匹配模式: (主量子数n)(角量子数l符号)(可选下划线)(标签)
        // 例如: '1s', '2pz', '3d_xz', '4f_xyz'
        const match = key.match(/^(\d+)([spdfghi])_?(.+)?$/);
        if (!match) {
            // 兼容无标签格式，如 '1s', '2p'
            const simpleMatch = key.match(/^(\d+)([spdfghi])$/);
            if (!simpleMatch) return null;
            const n = parseInt(simpleMatch[1]), l = ['s', 'p', 'd', 'f', 'g', 'h', 'i'].indexOf(simpleMatch[2]);
            return R(n, l, 0, 'c');
        }

        const n = parseInt(match[1]);
        const lChar = match[2];
        const label = match[3] || '';
        const l = ['s', 'p', 'd', 'f', 'g', 'h', 'i'].indexOf(lChar);

        // 标签到 (m, t) 的通用映射表
        const labelMap = {
            's': { m: 0, t: 'c' },
            'z': { m: 0, t: 'c' }, 'x': { m: 1, t: 'c' }, 'y': { m: 1, t: 's' }, // p 轨道
            'z2': { m: 0, t: 'c' }, 'xz': { m: 1, t: 'c' }, 'yz': { m: 1, t: 's' }, 'xy': { m: 2, t: 's' }, 'x2-y2': { m: 2, t: 'c' }, // d 轨道
            'z3': { m: 0, t: 'c' }, 'xz2': { m: 1, t: 'c' }, 'yz2': { m: 1, t: 's' }, // f 轨道
            'z(x2-y2)': { m: 2, t: 'c' }, 'xyz': { m: 2, t: 's' }, 'x(x2-3y2)': { m: 3, t: 'c' }, 'y(3x2-y2)': { m: 3, t: 's' },
            'z4': { m: 0, t: 'c' } // g 轨道初步支持
        };

        // 处理 pz, px, py 这种带前缀的情况
        let cleanLabel = label;
        if (l === 1 && label.length > 1 && label.startsWith('p')) cleanLabel = label.substring(1);

        const config = labelMap[cleanLabel] || labelMap[label];
        if (!config) return R(n, l, 0, 'c');
        return R(n, l, config.m, config.t);
    }

    const PhysicsCore = {
        factorial, binomialInt, generalizedLaguerre, associatedLegendre,
        matMul, matTranspose, jacobiSVD,
        Ylm_complex, realYlm_value, realYlm_abs2, getOrbitalKey, slaterRadialR, radialR, radialPDF,
        optimizeThomson, generateConstrainedDirections, sortOrbitalsForHybridization, getHybridCoefficients,
        allHybridOrbitalsDensity3D, singleHybridDensity3D, singleHybridWavefunction, hybridEstimateMaxDensity, orbitalParamsFromKey
    };

    if (typeof module !== 'undefined' && module.exports) module.exports = PhysicsCore;
    else if (typeof define === 'function' && define.amd) define(() => PhysicsCore);
    else root.PhysicsCore = PhysicsCore;

})(typeof self !== 'undefined' ? self : this);
