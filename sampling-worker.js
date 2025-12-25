/**
 * Web Worker for Monte Carlo Sampling
 * 将计算密集型的蒙特卡洛采样移到后台线程，避免阻塞主线程渲染
 */

// 核心物理数据和逻辑
importScripts('slater_basis.js');
importScripts('physics-core.js');
const core = PhysicsCore;

const A0 = 1; // Bohr radius unit
const PI = Math.PI;
const TWO_PI = 2 * PI;

// 别名：保持与 Worker 内部旧代码兼容
const optimizeThomsonWorker = core.optimizeThomson;
const buildDirectionMatrixWorker = core.buildDirectionMatrix;
const jacobiSVDWorker = core.jacobiSVD;
const matMulWorker = core.matMul;
const matTransposeWorker = core.matTranspose;
const sortOrbitalsForHybridizationWorker = core.sortOrbitalsForHybridization;
const getHybridCoefficientsWorker = (paramsList) => core.getHybridCoefficients(paramsList);

const factorial = core.factorial;
const binomialInt = core.binomialInt;
const generalizedLaguerre = core.generalizedLaguerre;
const associatedLegendre = core.associatedLegendre;
const Ylm_complex = core.Ylm_complex;
const Ylm_abs2 = core.Ylm_abs2;
const realYlm_abs2 = core.realYlm_abs2;
const realYlm_value = core.realYlm_value;
const getOrbitalKey = core.getOrbitalKey;
const slaterRadialR = core.slaterRadialR;
const radialR = core.radialR;
const radialPDF = core.radialPDF;
const density3D_real = core.density3D_real;
const hybridDensity3D = core.allHybridOrbitalsDensity3D;
const hybridEstimateMaxDensity = core.hybridEstimateMaxDensity;
const orbitalParamsFromKey = core.orbitalParamsFromKey;

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
