/**
 * 共享物理核心模块 (v2.0)
 *
 * 本文件包含 physics.js 与 sampling-worker.js 共用的核心物理计算函数，
 * 作为项目物理逻辑的唯一来源。
 *
 * 重要：修改本文件会同时影响主线程与 Worker，请务必进行全量回归测试。
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

    // ==================== 2. 线性代数工具 & 四元数 ====================

    // 四元数工具
    function quatMultiply(p, q) {
        return [
            p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3], // w
            p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2], // x
            p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1], // y
            p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0]  // z
        ];
    }

    function quatRotateVector(q, v) {
        // v' = q * v * q_inv
        // vector v as quaternion [0, v]
        const vq = [0, v[0], v[1], v[2]];
        const qInv = [q[0], -q[1], -q[2], -q[3]]; // q must be normalized
        const temp = quatMultiply(q, vq);
        const res = quatMultiply(temp, qInv);
        return [res[1], res[2], res[3]];
    }

    function quatNormalize(q) {
        const norm = Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
        if (norm < 1e-9) return [1, 0, 0, 0];
        return [q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm];
    }

    function quatRandom() {
        // Uniform random quaternion
        const u1 = Math.random(), u2 = Math.random(), u3 = Math.random();
        const sq10 = Math.sqrt(1 - u1), squ1 = Math.sqrt(u1);
        return [
            sq10 * Math.sin(TWO_PI * u2),
            sq10 * Math.cos(TWO_PI * u2),
            squ1 * Math.sin(TWO_PI * u3),
            squ1 * Math.cos(TWO_PI * u3)
        ];
    }

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

    // 杂化方向缓存：基于轨道参数的完整 key
    const _hybridDirCache = {};

    /**
     * 从轨道角向量子数 (l, m, t) 数值计算极值方向
     * 通过在单位球面上搜索 realYlm_value 的最大值
     */
    function computeOrbitalPeakDirection(l, m, t) {
        if (l === 0) return null; // s 轨道各向同性，无主轴

        let maxVal = -Infinity, bestT = 0, bestP = 0;
        // 粗搜索
        for (let theta = 0.05; theta < Math.PI; theta += 0.1) {
            for (let phi = 0; phi < TWO_PI; phi += 0.1) {
                const val = realYlm_value(l, m, t, theta, phi);
                if (val > maxVal) { maxVal = val; bestT = theta; bestP = phi; }
            }
        }
        // 细搜索
        for (let theta = Math.max(0.01, bestT - 0.15); theta < Math.min(Math.PI - 0.01, bestT + 0.15); theta += 0.02) {
            for (let phi = bestP - 0.15; phi < bestP + 0.15; phi += 0.02) {
                const val = realYlm_value(l, m, t, theta, phi);
                if (val > maxVal) { maxVal = val; bestT = theta; bestP = phi; }
            }
        }
        return [
            Math.sin(bestT) * Math.cos(bestP),
            Math.sin(bestT) * Math.sin(bestP),
            Math.cos(bestT)
        ];
    }

    /**
     * 计算旋转四元数，将向量 from 旋转到向量 to
     */
    function quatFromVectors(from, to) {
        const dot = from[0] * to[0] + from[1] * to[1] + from[2] * to[2];
        if (dot > 0.9999) return [1, 0, 0, 0]; // 无需旋转
        if (dot < -0.9999) {
            // 180度旋转，找一个垂直轴
            let axis = [0, 1, 0];
            if (Math.abs(from[1]) > 0.9) axis = [1, 0, 0];
            return quatNormalize([0, axis[0], axis[1], axis[2]]);
        }
        const cross = [
            from[1] * to[2] - from[2] * to[1],
            from[2] * to[0] - from[0] * to[2],
            from[0] * to[1] - from[1] * to[0]
        ];
        return quatNormalize([1 + dot, cross[0], cross[1], cross[2]]);
    }

    /**
     * 检验杂化轨道组合是否为支持的标准构型
     * 只允许: sp, sp2, sp3, sp3d, sp3d2
     * @param {Array} orbitalParams - 轨道参数数组
     * @returns {boolean} - 是否为有效的杂化组合
     */
    function isValidHybridization(orbitalParams) {
        const n = orbitalParams.length;
        const counts = { s: 0, p: 0, d: 0, f: 0 };
        orbitalParams.forEach(p => {
            if (p.angKey.l === 0) counts.s++;
            else if (p.angKey.l === 1) counts.p++;
            else if (p.angKey.l === 2) counts.d++;
            else counts.f++;
        });

        // Allowed: sp, sp2, sp3, sp3d, sp3d2
        if (n === 2 && counts.s === 1 && counts.p === 1) return true; // sp
        if (n === 3 && counts.s === 1 && counts.p === 2) return true; // sp2
        if (n === 4 && counts.s === 1 && counts.p === 3) return true; // sp3
        if (n === 5 && counts.s === 1 && counts.p === 3 && counts.d === 1) return true; // sp3d
        if (n === 6 && counts.s === 1 && counts.p === 3 && counts.d === 2) return true; // sp3d2

        return false; // 其他组合不支持
    }

    function generateConstrainedDirections(orbitalParams) {
        // 【v11.0 终极方案】硬编码优先 + 通用兜底
        // 用户指令：允许硬编码以保证教科书式的完美构型。

        const n = orbitalParams.length;

        // 1. 尝试匹配标准杂化 (Hardcoded Lookup)
        const counts = { s: 0, p: 0, d: 0, f: 0 };
        const p_axes = []; // To track which p orbitals are used (for sp/sp2 alignment)

        orbitalParams.forEach(p => {
            if (p.angKey.l === 0) counts.s++;
            else if (p.angKey.l === 1) {
                counts.p++;
                // Identify axis: pz(m=0), px(m=1,c), py(m=1,s)
                // Simplified check:
                if (p.angKey.m === 0) p_axes.push([0, 0, 1]);
                else if (p.angKey.t === 'c') p_axes.push([1, 0, 0]);
                else p_axes.push([0, 1, 0]);
            }
            else if (p.angKey.l === 2) counts.d++;
            else counts.f++;
        });

        // sp3 (tetrahedral) -> [1,1,1] family
        if (n === 4 && counts.s === 1 && counts.p === 3) {
            const sqrt3 = 1 / Math.sqrt(3);
            return [
                [sqrt3, sqrt3, sqrt3],
                [sqrt3, -sqrt3, -sqrt3],
                [-sqrt3, sqrt3, -sqrt3],
                [-sqrt3, -sqrt3, sqrt3]
            ];
        }

        // sp2 (trigonal planar) -> Aligned to the plane of the two p-orbitals
        if (n === 3 && counts.s === 1 && counts.p === 2) {
            // Determine plane normal (cross product of the two p axes)
            // But usually it's xy plane (px, py).
            // Let's perform a simple robust generation:
            // If we have px and py, generate in XY plane.
            // If px and pz, generate in XZ plane.
            // Standard Reference: x-axis aligned.

            // Default to XY if ambig, or use p-axes to define plane.
            // But for simplicity, let's just return standard XY plane if typical case.
            // Setup standard trigonal: 0, 120, 240.
            const base = [[1, 0, 0], [-0.5, 0.866, 0], [-0.5, -0.866, 0]];

            // If p-orbitals are NOT px,py, we might need rotation. 
            // But let's trust the general algorithm for weird sp2 (like s+py+pz).
            // Only hardcode s+px+py.
            const hasPx = p_axes.some(v => Math.abs(v[0]) > 0.9);
            const hasPy = p_axes.some(v => Math.abs(v[1]) > 0.9);
            if (hasPx && hasPy) return base;
        }

        // sp (linear) -> Aligned to the p-orbital axis
        if (n === 2 && counts.s === 1 && counts.p === 1) {
            const axis = p_axes[0]; // The p-orbital direction
            return [axis, [-axis[0], -axis[1], -axis[2]]];
        }

        // sp3d (trigonal bipyramidal)
        if (n === 5 && counts.s === 1 && counts.p === 3 && counts.d === 1) {
            // Axial (z), Equatorial (xy)
            return [
                [0, 0, 1], [0, 0, -1],
                [1, 0, 0], [-0.5, 0.866, 0], [-0.5, -0.866, 0]
            ];
        }

        // sp3d2 (octahedral)
        if (n === 6 && counts.s === 1 && counts.p === 3 && counts.d === 2) {
            return [
                [0, 0, 1], [0, 0, -1],
                [1, 0, 0], [-1, 0, 0],
                [0, 1, 0], [0, -1, 0]
            ];
        }

        // px + py (Orthogonal)
        if (n === 2 && counts.s === 0 && counts.p === 2) {
            // Just return the p-axes themselves
            return [p_axes[0], p_axes[1]];
        }


        // 生成缓存 key (通用算法兜底)
        const cacheKey = orbitalParams.map(p =>
            `${p.angKey.l}_${p.angKey.m}_${p.angKey.t}`
        ).join('|');

        if (_hybridDirCache[cacheKey]) {
            return JSON.parse(JSON.stringify(_hybridDirCache[cacheKey]));
        }

        // 1. 计算锚点（物理极值方向）
        const anchors = [];
        for (const p of orbitalParams) {
            const dir = computeOrbitalPeakDirection(p.angKey.l, p.angKey.m, p.angKey.t);
            if (dir) anchors.push(dir);
        }

        // 2. 尝试生成 Thomson 几何（理想斥力模型）
        const thomsonPoints = optimizeThomson(n);

        // 3. 【核心判据】物理可行性检验 (Rank Check)
        // 检查选定的轨道基组是否在数学上能够支撑起 Thomson 几何
        const A = [];
        for (let i = 0; i < n; i++) {
            const p = thomsonPoints[i];
            const theta = Math.acos(p[2]);
            const phi = Math.atan2(p[1], p[0]);
            const row = orbitalParams.map(param =>
                realYlm_value(param.angKey.l, param.angKey.m, param.angKey.t, theta, phi)
            );
            A.push(row);
        }

        const { S } = jacobiSVD(A);
        const minSingularValue = Math.min(...S); // 最小奇异值
        const isGeometricallyCompatible = minSingularValue > 1e-4;

        if (!isGeometricallyCompatible) {
            // [CASE A] 几何不兼容 -> 回退到锚点（自然正交基）
            if (anchors.length !== n) {
                _hybridDirCache[cacheKey] = thomsonPoints; // 兜底
                return thomsonPoints;
            }
            _hybridDirCache[cacheKey] = JSON.parse(JSON.stringify(anchors));
            return anchors;
        }

        // [CASE B] 几何兼容 -> 使用 Thomson + 对齐优化
        let bestPoints = thomsonPoints;

        if (anchors.length > 0) {
            let bestScore = -Infinity;
            let bestR = [1, 0, 0, 0];

            // 扩展锚点
            const fullAnchors = [];
            anchors.forEach(a => {
                fullAnchors.push(a);
                fullAnchors.push([-a[0], -a[1], -a[2]]);
            });

            // 全局搜索对齐
            for (let i = 0; i < n; i++) {
                const P = thomsonPoints[i];
                for (const A of fullAnchors) {
                    const qBase = quatFromVectors(P, A);

                    for (let angle = 0; angle < TWO_PI; angle += 0.174) { // 10 degree step
                        const half = angle / 2;
                        const s = Math.sin(half);
                        const qRot = [Math.cos(half), A[0] * s, A[1] * s, A[2] * s];
                        const q = quatMultiply(qRot, qBase);

                        const rotated = thomsonPoints.map(pt => quatRotateVector(q, pt));

                        // v10.0: 对称性优先评分
                        const dots = [];
                        let sumMaxDot = 0;
                        for (const rp of rotated) {
                            let maxDot = 0;
                            for (const fa of fullAnchors) {
                                const d = rp[0] * fa[0] + rp[1] * fa[1] + rp[2] * fa[2];
                                if (d > maxDot) maxDot = d;
                            }
                            dots.push(maxDot);
                            sumMaxDot += maxDot;
                        }

                        // 计算标准差
                        const mean = sumMaxDot / n;
                        let sqDiffSum = 0;
                        for (let k = 0; k < n; k++) sqDiffSum += (dots[k] - mean) * (dots[k] - mean);
                        const stdev = Math.sqrt(sqDiffSum / n);

                        // Score = Sum - 100 * Stdev
                        let score = sumMaxDot - 100 * stdev;

                        if (score > bestScore) {
                            bestScore = score;
                            bestR = q;
                        }
                    }
                }
            }
            bestPoints = thomsonPoints.map(p => quatRotateVector(bestR, p));
        }

        _hybridDirCache[cacheKey] = JSON.parse(JSON.stringify(bestPoints));
        return bestPoints;
    }
    // 【核心优化器】基于爬山法的四元数旋转优化
    function optimizeHybridAlignment(initialPoints, orbitalParams) {
        const N = initialPoints.length;
        const M = orbitalParams.length; // usually N=M

        // 目标函数：计算当前旋转下的 Vol Index
        function calculateScore(q) {
            // Apply rotation
            const rotated = initialPoints.map(p => quatRotateVector(q, p));

            // Build Matrix A (N x M)
            const A = [];
            for (let i = 0; i < N; i++) {
                const p = rotated[i];
                const r = Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                // Robust acos
                const theta = Math.acos(Math.max(-1, Math.min(1, p[2] / r)));
                const phi = Math.atan2(p[1], p[0]);

                const row = orbitalParams.map(param => {
                    const { l, m, t } = param.angKey;
                    return realYlm_value(l, m, t, theta, phi);
                });
                A.push(row);
            }

            // SVD to get Singular Values
            const { S } = jacobiSVD(A);

            // Score = Sum(log(sigma + eps))
            // 避免奇异值接近0导致的数值不稳定
            let score = 0;
            for (let k = 0; k < S.length; k++) {
                score += Math.log(S[k] + 1e-12);
            }
            return score;
        }

        let bestScore = -Infinity;
        let bestPoints = initialPoints; // Fallback

        // Random Restart Hill Climbing
        // N=20 restarts seems sufficient for low dim
        const NUM_RESTARTS = 20;
        const STEPS_PER_RUN = 50;

        for (let run = 0; run < NUM_RESTARTS; run++) {
            // Random start
            let currentQ = quatRandom();
            let currentScore = calculateScore(currentQ);

            // Simple Hill Climbing
            let stepSize = 0.2;
            for (let step = 0; step < STEPS_PER_RUN; step++) {
                // Try perturbing quaternion
                const perturb = [
                    (Math.random() - 0.5) * stepSize,
                    (Math.random() - 0.5) * stepSize,
                    (Math.random() - 0.5) * stepSize,
                    (Math.random() - 0.5) * stepSize
                ];
                const trialQ = quatNormalize([
                    currentQ[0] + perturb[0], currentQ[1] + perturb[1],
                    currentQ[2] + perturb[2], currentQ[3] + perturb[3]
                ]);

                const trialScore = calculateScore(trialQ);

                if (trialScore > currentScore) {
                    currentScore = trialScore;
                    currentQ = trialQ;
                } else {
                    // Decay step size if no improvement
                    stepSize *= 0.95;
                }
            }

            if (currentScore > bestScore) {
                bestScore = currentScore;
                bestPoints = initialPoints.map(p => quatRotateVector(currentQ, p));
            }
        }

        return bestPoints;
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

    // 杂化系数矩阵缓存
    const _hybridCoeffCache = {};

    function getHybridCoefficients(orbitalParams) {
        // 生成缓存 key
        const cacheKey = orbitalParams.map(p =>
            `${p.n}_${p.l}_${p.angKey.l}_${p.angKey.m}_${p.angKey.t}`
        ).join('|');

        if (_hybridCoeffCache[cacheKey]) {
            return _hybridCoeffCache[cacheKey];
        }

        const sorted = sortOrbitalsForHybridization(orbitalParams);
        const directions = generateConstrainedDirections(sorted);
        const A = buildDirectionMatrix(directions, sorted);
        const { U, V } = jacobiSVD(A);
        const result = matMul(U, matTranspose(V));

        _hybridCoeffCache[cacheKey] = result;
        return result;
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
        allHybridOrbitalsDensity3D, singleHybridDensity3D, singleHybridWavefunction, hybridEstimateMaxDensity, orbitalParamsFromKey,
        computeOrbitalPeakDirection,
        quatFromVectors, quatRotateVector, quatNormalize, quatMultiply,
        isValidHybridization
    };

    if (typeof module !== 'undefined' && module.exports) module.exports = PhysicsCore;
    else if (typeof define === 'function' && define.amd) define(() => PhysicsCore);
    else root.PhysicsCore = PhysicsCore;

})(typeof self !== 'undefined' ? self : this);
