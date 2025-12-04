/**
 * Web Worker for Monte Carlo Sampling
 * 将计算密集型的蒙特卡洛采样移到后台线程，避免阻塞主线程渲染
 */

// ==================== 物理计算函数（从 physics.js 复制，Worker 无法访问主线程代码）====================

const A0 = 1; // Bohr radius unit
const PI = Math.PI;

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
function realYlm_abs2(l, m, type, theta, phi) {
    if (m === 0) {
        return Ylm_abs2(l, 0, theta, phi);
    }
    const mm = Math.abs(m);
    const y = Ylm_complex(l, mm, theta, phi);
    if (type === 'c') {
        const v = Math.SQRT2 * y.re;
        return v * v;
    } else {
        const v = Math.SQRT2 * y.im;
        return v * v;
    }
}

// 实球谐函数（值）
function realYlm_value(l, m, type, theta, phi) {
    if (m === 0) {
        const y = Ylm_complex(l, 0, theta, phi);
        return y.re;
    }
    const mm = Math.abs(m);
    const y = Ylm_complex(l, mm, theta, phi);
    if (type === 'c') {
        return Math.SQRT2 * y.re;
    } else {
        return Math.SQRT2 * y.im;
    }
}

// 归一化径向函数 R_nl(r)
function radialR(n, l, r, Z = 1, a0 = A0) {
    if (n <= 0 || l < 0 || l >= n) return 0;
    const rho = (2 * Z * r) / (n * a0);
    const k = n - l - 1;
    if (k < 0) return 0;
    const pref = Math.pow(2 * Z / (n * a0), 1.5) * Math.sqrt(factorial(n - l - 1) / (2 * n * factorial(n + l)));
    const poly = generalizedLaguerre(k, 2 * l + 1, rho);
    return pref * Math.exp(-rho / 2) * Math.pow(rho, l) * poly;
}

// 3D 概率密度
function density3D_real(angKey, n, l, r, theta, phi, Z = 1, a0 = A0) {
    const R = radialR(n, l, r, Z, a0);
    let Y2 = 1 / (4 * PI);
    if (angKey && typeof angKey === 'object') {
        Y2 = realYlm_abs2(angKey.l, angKey.m, angKey.t, theta, phi);
    }
    return (R * R) * Y2;
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
function radialPDF(n, l, r, Z = 1, a0 = A0) {
    const R = radialR(n, l, r, Z, a0);
    return r * r * R * R;
}

// CDF 表缓存
const _cdfCache = {};

/**
 * 构建径向分布的累积分布函数表
 */
function buildRadialCDF(n, l, numPoints = 2000) {
    const key = `${n}_${l}`;
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
        const P_prev = radialPDF(n, l, r[i - 1]);
        const P_curr = radialPDF(n, l, r[i]);
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
function sampleRadialExact(n, l) {
    const { r, cdf, numPoints } = buildRadialCDF(n, l);

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
function importanceSample(n, l, angKey, samplingBoundary) {
    // ==================== 第一步：径向采样（精确逆CDF） ====================
    let r = sampleRadialExact(n, l);

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

// ==================== 采样逻辑 ====================

// 是否启用重要性采样
let useImportanceSampling = true;

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
        phaseOn,            // 是否显示相位
        compareColors       // 比照模式颜色
    } = config;

    const paramsList = orbitals.map(k => orbitalParamsFromKey(k)).filter(Boolean);
    if (!paramsList.length) {
        return { points: [], samples: [], farthestDistance: 0 };
    }

    const samplingVolumeSize = samplingBoundary * 2;
    const points = [];       // 采样到的点 [{x,y,z,r,g,b,orbitalIndex}]
    const samples = [];      // 用于图表的样本数据 [{r,theta,orbitalKey}]
    let farthestDistance = 0;
    let attempts = 0;

    // 【关键修复】多轨道模式下，确保每个轨道获得相同数量的最终采样点
    // 而不是相同数量的尝试次数（因为不同轨道的接受率不同）
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

        // 使用重要性采样（针对当前轨道优化）
        if (useImportanceSampling) {
            const result = importanceSample(p.n, p.l, p.angKey, samplingBoundary);

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
                const density = density3D_real(p.angKey, p.n, p.l, r, theta, phi);
                const scaleFactor = Math.min(200, 3.0 * Math.pow(p.n, 3));
                accepted = Math.random() < density * scaleFactor;
            } else {
                let probability = 0;
                let minN = paramsList[0].n;
                for (let i = 0; i < paramsList.length; i++) {
                    const params = paramsList[i];
                    probability += density3D_real(params.angKey, params.n, params.l, r, theta, phi);
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
                        const R = radialR(p.n, p.l, r);
                        const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                        psi = R * Y;
                    } else {
                        // 单选模式：计算所有轨道的叠加相位
                        for (const params of paramsList) {
                            const R = radialR(params.n, params.l, r);
                            const Y = realYlm_value(params.angKey.l, params.angKey.m, params.angKey.t, theta, phi);
                            psi += R * Y;
                        }
                    }
                    if (psi > 0) {
                        r_color = 1; g_color = 0.2; b_color = 0.2;
                    } else if (psi < 0) {
                        r_color = 0.2; g_color = 0.2; b_color = 1;
                    } else {
                        r_color = 1; g_color = 1; b_color = 1;
                    }
                } else {
                    r_color = 1; g_color = 1; b_color = 1;
                }
            }

            point = { x, y, z, r: r_color, g: g_color, b: b_color, orbitalIndex: isIndependentMode ? orbitalIndex : -1 };
            samples.push({ r, theta, orbitalKey: isIndependentMode ? orbitals[orbitalIndex] : null });
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
