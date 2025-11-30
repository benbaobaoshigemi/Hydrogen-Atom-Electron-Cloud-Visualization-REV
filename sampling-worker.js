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

// 轨道参数映射
function orbitalParamsFromKey(key) {
    const R = (n, l, m, t) => ({ n, l, angKey: { l, m, t } });
    switch (key) {
        case '1s': return R(1, 0, 0, 'c');
        case '2s': return R(2, 0, 0, 'c');
        case '2pz': return R(2, 1, 0, 'c');
        case '2px': return R(2, 1, 1, 'c');
        case '2py': return R(2, 1, 1, 's');
        case '3s': return R(3, 0, 0, 'c');
        case '3pz': return R(3, 1, 0, 'c');
        case '3px': return R(3, 1, 1, 'c');
        case '3py': return R(3, 1, 1, 's');
        case '3d_z2': return R(3, 2, 0, 'c');
        case '3d_xz': return R(3, 2, 1, 'c');
        case '3d_yz': return R(3, 2, 1, 's');
        case '3d_xy': return R(3, 2, 2, 's');
        case '3d_x2-y2': return R(3, 2, 2, 'c');
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
        default: return R(1, 0, 0, 'c');
    }
}

// ==================== 采样逻辑 ====================

/**
 * 批量采样函数 - 在 Worker 中执行
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

    while (attempts < maxAttempts && points.length < targetPoints) {
        attempts++;

        const x = (Math.random() - 0.5) * samplingVolumeSize;
        const y = (Math.random() - 0.5) * samplingVolumeSize;
        const z = (Math.random() - 0.5) * samplingVolumeSize;
        const r = Math.sqrt(x * x + y * y + z * z);

        if (r === 0) continue;

        const theta = Math.acos(z / r);
        const phi = Math.atan2(y, x);

        let accepted = false;
        let point = null;

        if (isIndependentMode) {
            // 独立模式：随机选择一个轨道
            const randomOrbitalIndex = Math.floor(Math.random() * paramsList.length);
            const p = paramsList[randomOrbitalIndex];
            const density = density3D_real(p.angKey, p.n, p.l, r, theta, phi);
            const scaleFactor = Math.min(200, 3.0 * Math.pow(p.n, 3));

            if (Math.random() < density * scaleFactor) {
                accepted = true;
                let r_color, g_color, b_color;

                if (isMultiselectMode) {
                    // 多选模式：支持相位
                    if (phaseOn) {
                        const R = radialR(p.n, p.l, r);
                        const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                        const psi = R * Y;
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
                } else {
                    // 比照模式：固定颜色
                    if (compareColors && randomOrbitalIndex < compareColors.length) {
                        const color = compareColors[randomOrbitalIndex];
                        r_color = color[0];
                        g_color = color[1];
                        b_color = color[2];
                    } else {
                        r_color = 1; g_color = 1; b_color = 1;
                    }
                }

                point = { x, y, z, r: r_color, g: g_color, b: b_color, orbitalIndex: randomOrbitalIndex };
                samples.push({ r, theta, orbitalKey: orbitals[randomOrbitalIndex] });
            }
        } else {
            // 正常模式：概率叠加
            let probability = 0;
            let minN = paramsList[0].n;

            for (let i = 0; i < paramsList.length; i++) {
                const p = paramsList[i];
                const density = density3D_real(p.angKey, p.n, p.l, r, theta, phi);
                probability += density;
                if (p.n < minN) minN = p.n;
            }

            const scaleFactor = Math.min(200, 3.0 * Math.pow(minN, 3) / Math.sqrt(paramsList.length));

            if (Math.random() < probability * scaleFactor) {
                accepted = true;
                let r_color, g_color, b_color;

                if (phaseOn) {
                    let psi = 0;
                    for (const p of paramsList) {
                        const R = radialR(p.n, p.l, r);
                        const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                        psi += R * Y;
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

                point = { x, y, z, r: r_color, g: g_color, b: b_color, orbitalIndex: -1 };
                samples.push({ r, theta, orbitalKey: null });
            }
        }

        if (accepted && point) {
            points.push(point);
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

self.onmessage = function(e) {
    const { type, data, taskId } = e.data;

    switch (type) {
        case 'sample':
            // 执行采样任务
            const result = performSampling(data);
            self.postMessage({
                type: 'sample-result',
                taskId: taskId,
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
