/**
 * 等值面计算 Web Worker
 * 异步执行 Marching Cubes 算法，避免阻塞主线程
 */

'use strict';

// =============== 物理计算函数（从physics.js复制） ===============

const PI = Math.PI;
const A0 = 1.0; // 玻尔半径（原子单位）

// Slater 规则参数（从physics.js复制）
const SlaterConfig = {
    H: { Z: 1, Z_eff: { '1s': 1.0 } },
    He: { Z: 2, Z_eff: { '1s': 1.6875 } },
    C: { Z: 6, Z_eff: { '1s': 5.6727, '2s': 3.2166, '2p': 3.1358 } },
    N: { Z: 7, Z_eff: { '1s': 6.6651, '2s': 3.8474, '2p': 3.8340 } },
    O: { Z: 8, Z_eff: { '1s': 7.6579, '2s': 4.4916, '2p': 4.4532 } },
    // 更多元素...
};

// 获取有效核电荷
// 【修复】H 原子的所有轨道 Z_eff = Z = 1（氢原子解析解）
function getZeff(atomType, n, l) {
    // H 原子是精确解析解，Z_eff = Z = 1 对所有轨道
    if (atomType === 'H') return 1.0;

    const config = SlaterConfig[atomType];
    if (!config) return 1.0;

    const shellKey = l === 0 ? `${n}s` : l === 1 ? `${n}p` : l === 2 ? `${n}d` : `${n}f`;
    return config.Z_eff[shellKey] || 1.0; // 默认返回 1.0 而非 Z/n
}

// 阶乘
const FACT = (() => {
    const f = [1];
    for (let i = 1; i <= 64; i++) f[i] = f[i - 1] * i;
    return f;
})();

function factorial(n) {
    return n < 0 ? 1 : (FACT[n] || 1);
}

// Laguerre 多项式（递推）
function laguerreL(n, alpha, x) {
    if (n === 0) return 1;
    if (n === 1) return 1 + alpha - x;
    let L0 = 1, L1 = 1 + alpha - x;
    for (let k = 2; k <= n; k++) {
        const L2 = ((2 * k - 1 + alpha - x) * L1 - (k - 1 + alpha) * L0) / k;
        L0 = L1; L1 = L2;
    }
    return L1;
}

// Legendre 多项式
function legendreP(l, m, x) {
    const absM = Math.abs(m);
    let pmm = 1.0;
    if (absM > 0) {
        const somx2 = Math.sqrt((1 - x) * (1 + x));
        let fact = 1.0;
        for (let i = 1; i <= absM; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l === absM) return pmm;
    let pmmp1 = x * (2 * absM + 1) * pmm;
    if (l === absM + 1) return pmmp1;
    let pll = 0;
    for (let ll = absM + 2; ll <= l; ll++) {
        pll = ((2 * ll - 1) * x * pmmp1 - (ll + absM - 1) * pmm) / (ll - absM);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

// 归一化的实球谐函数
function realYlm_value(l, m, type, theta, phi) {
    const cosT = Math.cos(theta);
    const Plm = legendreP(l, Math.abs(m), cosT);

    // 归一化系数
    let norm;
    if (m === 0) {
        norm = Math.sqrt((2 * l + 1) / (4 * PI));
    } else {
        norm = Math.sqrt((2 * l + 1) / (2 * PI) * factorial(l - Math.abs(m)) / factorial(l + Math.abs(m)));
    }

    let Y = norm * Plm;
    if (m !== 0) {
        if (type === 'c') {
            Y *= Math.cos(Math.abs(m) * phi);
        } else {
            Y *= Math.sin(Math.abs(m) * phi);
        }
    }

    // Condon-Shortley 相位
    if (m > 0 && m % 2 !== 0) Y = -Y;

    return Y;
}

// 获取轨道键名 (如 "2s", "3p")
function getOrbitalKey(n, l) {
    const lNames = ['s', 'p', 'd', 'f', 'g'];
    return `${n}${lNames[l] || 's'}`;
}

// Slater Type Orbital 径向函数
function slaterRadialR(basis, r) {
    if (!basis) return 0;
    let sum = 0;
    for (let i = 0; i < basis.length; i++) {
        const { nStar, zeta, coeff } = basis[i];
        const nFact2 = factorial(2 * nStar);
        const norm = Math.pow(2 * zeta, nStar + 0.5) / Math.sqrt(nFact2);
        const val = coeff * norm * Math.pow(r, nStar - 1) * Math.exp(-zeta * r);
        sum += val;
    }
    return sum;
}

// 全局变量：存储 SlaterBasis 数据（由主线程传入）
let workerSlaterBasis = null;

// 径向波函数 R(r) - 支持 STO 策略
function radialR(n, l, r, Z, a0, atomType, slaterBasis) {
    // 策略 A: Slater Type Orbitals (STO) - 用于非氢原子
    if (atomType && atomType !== 'H' && slaterBasis) {
        const orbitalKey = getOrbitalKey(n, l);
        const orbitals = slaterBasis.orbitals;
        if (orbitals && orbitals[orbitalKey]) {
            const basis = orbitals[orbitalKey];
            const val = slaterRadialR(basis, r);

            // 智能相位校正
            const analyticalSign = ((n - l - 1) % 2 === 0) ? 1 : -1;
            let dominantTerm = basis[0];
            for (let i = 1; i < basis.length; i++) {
                if (basis[i].zeta < dominantTerm.zeta) dominantTerm = basis[i];
            }
            const stoSign = Math.sign(dominantTerm.coeff);
            if (stoSign * analyticalSign < 0) return -val;
            return val;
        }
    }

    // 策略 B: 氢原子解析解（默认）
    const Zeff = getZeff(atomType, n, l);
    const rho = 2 * Zeff * r / (n * a0);
    const normFactor = Math.sqrt(
        Math.pow(2 * Zeff / (n * a0), 3) * factorial(n - l - 1) / (2 * n * factorial(n + l))
    );
    const lag = laguerreL(n - l - 1, 2 * l + 1, rho);
    return normFactor * Math.exp(-rho / 2) * Math.pow(rho, l) * lag;
}

// =============== Marching Cubes 算法（从marching_cubes.js复制） ===============

const TRI_TABLE = [
    [],
    [0, 8, 3],
    [0, 1, 9],
    [1, 8, 3, 9, 8, 1],
    [1, 2, 10],
    [0, 8, 3, 1, 2, 10],
    [9, 2, 10, 0, 2, 9],
    [2, 8, 3, 2, 10, 8, 10, 9, 8],
    [3, 11, 2],
    [0, 11, 2, 8, 11, 0],
    [1, 9, 0, 2, 3, 11],
    [1, 11, 2, 1, 9, 11, 9, 8, 11],
    [3, 10, 1, 11, 10, 3],
    [0, 10, 1, 0, 8, 10, 8, 11, 10],
    [3, 9, 0, 3, 11, 9, 11, 10, 9],
    [9, 8, 10, 10, 8, 11],
    [4, 7, 8],
    [4, 3, 0, 7, 3, 4],
    [0, 1, 9, 8, 4, 7],
    [4, 1, 9, 4, 7, 1, 7, 3, 1],
    [1, 2, 10, 8, 4, 7],
    [3, 4, 7, 3, 0, 4, 1, 2, 10],
    [9, 2, 10, 9, 0, 2, 8, 4, 7],
    [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4],
    [8, 4, 7, 3, 11, 2],
    [11, 4, 7, 11, 2, 4, 2, 0, 4],
    [9, 0, 1, 8, 4, 7, 2, 3, 11],
    [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1],
    [3, 10, 1, 3, 11, 10, 7, 8, 4],
    [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4],
    [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3],
    [4, 7, 11, 4, 11, 9, 9, 11, 10],
    [9, 5, 4],
    [9, 5, 4, 0, 8, 3],
    [0, 5, 4, 1, 5, 0],
    [8, 5, 4, 8, 3, 5, 3, 1, 5],
    [1, 2, 10, 9, 5, 4],
    [3, 0, 8, 1, 2, 10, 4, 9, 5],
    [5, 2, 10, 5, 4, 2, 4, 0, 2],
    [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8],
    [9, 5, 4, 2, 3, 11],
    [0, 11, 2, 0, 8, 11, 4, 9, 5],
    [0, 5, 4, 0, 1, 5, 2, 3, 11],
    [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5],
    [10, 3, 11, 10, 1, 3, 9, 5, 4],
    [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10],
    [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3],
    [5, 4, 8, 5, 8, 10, 10, 8, 11],
    [9, 7, 8, 5, 7, 9],
    [9, 3, 0, 9, 5, 3, 5, 7, 3],
    [0, 7, 8, 0, 1, 7, 1, 5, 7],
    [1, 5, 3, 3, 5, 7],
    [9, 7, 8, 9, 5, 7, 10, 1, 2],
    [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3],
    [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2],
    [2, 10, 5, 2, 5, 3, 3, 5, 7],
    [7, 9, 5, 7, 8, 9, 3, 11, 2],
    [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11],
    [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7],
    [11, 2, 1, 11, 1, 7, 7, 1, 5],
    [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11],
    [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0],
    [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0],
    [11, 10, 5, 7, 11, 5],
    [10, 6, 5],
    [0, 8, 3, 5, 10, 6],
    [9, 0, 1, 5, 10, 6],
    [1, 8, 3, 1, 9, 8, 5, 10, 6],
    [1, 6, 5, 2, 6, 1],
    [1, 6, 5, 1, 2, 6, 3, 0, 8],
    [9, 6, 5, 9, 0, 6, 0, 2, 6],
    [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8],
    [2, 3, 11, 10, 6, 5],
    [11, 0, 8, 11, 2, 0, 10, 6, 5],
    [0, 1, 9, 2, 3, 11, 5, 10, 6],
    [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11],
    [6, 3, 11, 6, 5, 3, 5, 1, 3],
    [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6],
    [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9],
    [6, 5, 9, 6, 9, 11, 11, 9, 8],
    [5, 10, 6, 4, 7, 8],
    [4, 3, 0, 4, 7, 3, 6, 5, 10],
    [1, 9, 0, 5, 10, 6, 8, 4, 7],
    [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4],
    [6, 1, 2, 6, 5, 1, 4, 7, 8],
    [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7],
    [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6],
    [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9],
    [3, 11, 2, 7, 8, 4, 10, 6, 5],
    [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11],
    [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6],
    [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6],
    [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6],
    [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11],
    [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7],
    [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9],
    [10, 4, 9, 6, 4, 10],
    [4, 10, 6, 4, 9, 10, 0, 8, 3],
    [10, 0, 1, 10, 6, 0, 6, 4, 0],
    [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10],
    [1, 4, 9, 1, 2, 4, 2, 6, 4],
    [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4],
    [0, 2, 4, 4, 2, 6],
    [8, 3, 2, 8, 2, 4, 4, 2, 6],
    [10, 4, 9, 10, 6, 4, 11, 2, 3],
    [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6],
    [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10],
    [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1],
    [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3],
    [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1],
    [3, 11, 6, 3, 6, 0, 0, 6, 4],
    [6, 4, 8, 11, 6, 8],
    [7, 10, 6, 7, 8, 10, 8, 9, 10],
    [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10],
    [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0],
    [10, 6, 7, 10, 7, 1, 1, 7, 3],
    [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7],
    [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9],
    [7, 8, 0, 7, 0, 6, 6, 0, 2],
    [7, 3, 2, 6, 7, 2],
    [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7],
    [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7],
    [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11],
    [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1],
    [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6],
    [0, 9, 1, 11, 6, 7],
    [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0],
    [7, 11, 6],
    [7, 6, 11],
    [3, 0, 8, 11, 7, 6],
    [0, 1, 9, 11, 7, 6],
    [8, 1, 9, 8, 3, 1, 11, 7, 6],
    [10, 1, 2, 6, 11, 7],
    [1, 2, 10, 3, 0, 8, 6, 11, 7],
    [2, 9, 0, 2, 10, 9, 6, 11, 7],
    [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8],
    [7, 2, 3, 6, 2, 7],
    [7, 0, 8, 7, 6, 0, 6, 2, 0],
    [2, 7, 6, 2, 3, 7, 0, 1, 9],
    [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6],
    [10, 7, 6, 10, 1, 7, 1, 3, 7],
    [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8],
    [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7],
    [7, 6, 10, 7, 10, 8, 8, 10, 9],
    [6, 8, 4, 11, 8, 6],
    [3, 6, 11, 3, 0, 6, 0, 4, 6],
    [8, 6, 11, 8, 4, 6, 9, 0, 1],
    [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6],
    [6, 8, 4, 6, 11, 8, 2, 10, 1],
    [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6],
    [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9],
    [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3],
    [8, 2, 3, 8, 4, 2, 4, 6, 2],
    [0, 4, 2, 4, 6, 2],
    [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8],
    [1, 9, 4, 1, 4, 2, 2, 4, 6],
    [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1],
    [10, 1, 0, 10, 0, 6, 6, 0, 4],
    [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3],
    [10, 9, 4, 6, 10, 4],
    [4, 9, 5, 7, 6, 11],
    [0, 8, 3, 4, 9, 5, 11, 7, 6],
    [5, 0, 1, 5, 4, 0, 7, 6, 11],
    [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5],
    [9, 5, 4, 10, 1, 2, 7, 6, 11],
    [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5],
    [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2],
    [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6],
    [7, 2, 3, 7, 6, 2, 5, 4, 9],
    [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7],
    [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0],
    [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8],
    [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7],
    [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4],
    [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10],
    [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10],
    [6, 9, 5, 6, 11, 9, 11, 8, 9],
    [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5],
    [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11],
    [6, 11, 3, 6, 3, 5, 5, 3, 1],
    [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6],
    [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10],
    [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5],
    [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3],
    [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2],
    [9, 5, 6, 9, 6, 0, 0, 6, 2],
    [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8],
    [1, 5, 6, 2, 1, 6],
    [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6],
    [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0],
    [0, 3, 8, 5, 6, 10],
    [10, 5, 6],
    [11, 5, 10, 7, 5, 11],
    [11, 5, 10, 11, 7, 5, 8, 3, 0],
    [5, 11, 7, 5, 10, 11, 1, 9, 0],
    [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1],
    [11, 1, 2, 11, 7, 1, 7, 5, 1],
    [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11],
    [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7],
    [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2],
    [2, 5, 10, 2, 3, 5, 3, 7, 5],
    [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5],
    [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2],
    [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2],
    [1, 3, 5, 3, 7, 5],
    [0, 8, 7, 0, 7, 1, 1, 7, 5],
    [9, 0, 3, 9, 3, 5, 5, 3, 7],
    [9, 8, 7, 5, 9, 7],
    [5, 8, 4, 5, 10, 8, 10, 11, 8],
    [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0],
    [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5],
    [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4],
    [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8],
    [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11],
    [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5],
    [9, 4, 5, 2, 11, 3],
    [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4],
    [5, 10, 2, 5, 2, 4, 4, 2, 0],
    [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9],
    [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2],
    [8, 4, 5, 8, 5, 3, 3, 5, 1],
    [0, 4, 5, 1, 0, 5],
    [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5],
    [9, 4, 5],
    [4, 11, 7, 4, 9, 11, 9, 10, 11],
    [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11],
    [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11],
    [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4],
    [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2],
    [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3],
    [11, 7, 4, 11, 4, 2, 2, 4, 0],
    [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4],
    [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9],
    [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7],
    [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10],
    [1, 10, 2, 8, 7, 4],
    [4, 9, 1, 4, 1, 7, 7, 1, 3],
    [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1],
    [4, 0, 3, 7, 4, 3],
    [4, 8, 7],
    [9, 10, 8, 10, 11, 8],
    [3, 0, 9, 3, 9, 11, 11, 9, 10],
    [0, 1, 10, 0, 10, 8, 8, 10, 11],
    [3, 1, 10, 11, 3, 10],
    [1, 2, 11, 1, 11, 9, 9, 11, 8],
    [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9],
    [0, 2, 11, 8, 0, 11],
    [3, 2, 11],
    [2, 3, 8, 2, 8, 10, 10, 8, 9],
    [9, 10, 2, 0, 9, 2],
    [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8],
    [1, 10, 2],
    [1, 3, 8, 9, 1, 8],
    [0, 9, 1],
    [0, 3, 8],
    []
];

const EDGE_TABLE = [
    [0, 1], [1, 2], [2, 3], [3, 0],
    [4, 5], [5, 6], [6, 7], [7, 4],
    [0, 4], [1, 5], [2, 6], [3, 7]
];

const VERTEX_OFFSETS = [
    [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]
];

function marchingCubes(calcPsi, bounds, resolution, isovalue, onProgress) {
    const [minX, minY, minZ] = bounds.min;
    const [maxX, maxY, maxZ] = bounds.max;
    const dx = (maxX - minX) / resolution;
    const dy = (maxY - minY) / resolution;
    const dz = (maxZ - minZ) / resolution;

    // 计算网格值
    const grid = new Float32Array((resolution + 1) ** 3);
    let idx = 0;
    const totalPoints = (resolution + 1) ** 3;

    for (let iz = 0; iz <= resolution; iz++) {
        for (let iy = 0; iy <= resolution; iy++) {
            for (let ix = 0; ix <= resolution; ix++) {
                grid[idx++] = calcPsi(minX + ix * dx, minY + iy * dy, minZ + iz * dz);
            }
        }
        // 每完成一层 z，报告进度
        if (onProgress && iz % 10 === 0) {
            onProgress(idx / totalPoints * 0.8); // 网格计算占 80%
        }
    }

    const positiveTriangles = [];
    const negativeTriangles = [];

    for (let iz = 0; iz < resolution; iz++) {
        for (let iy = 0; iy < resolution; iy++) {
            for (let ix = 0; ix < resolution; ix++) {
                const baseX = minX + ix * dx;
                const baseY = minY + iy * dy;
                const baseZ = minZ + iz * dz;

                const values = [];
                for (let v = 0; v < 8; v++) {
                    const [vx, vy, vz] = VERTEX_OFFSETS[v];
                    const gi = (iz + vz) * (resolution + 1) * (resolution + 1) +
                        (iy + vy) * (resolution + 1) + (ix + vx);
                    values[v] = grid[gi];
                }

                processCube(values, isovalue, baseX, baseY, baseZ, dx, dy, dz, positiveTriangles);
                processCube(values, -isovalue, baseX, baseY, baseZ, dx, dy, dz, negativeTriangles);
            }
        }
    }

    if (onProgress) onProgress(1.0);

    return { positive: positiveTriangles, negative: negativeTriangles };
}

function processCube(values, threshold, baseX, baseY, baseZ, dx, dy, dz, output) {
    let cubeIndex = 0;
    for (let i = 0; i < 8; i++) {
        if (values[i] > threshold) cubeIndex |= (1 << i);
    }
    if (cubeIndex === 0 || cubeIndex === 255) return;

    const edgeVertices = [];
    for (let e = 0; e < 12; e++) {
        const [v0, v1] = EDGE_TABLE[e];
        const val0 = values[v0], val1 = values[v1];
        if ((val0 > threshold) !== (val1 > threshold)) {
            const t = (threshold - val0) / (val1 - val0);
            const [ox0, oy0, oz0] = VERTEX_OFFSETS[v0];
            const [ox1, oy1, oz1] = VERTEX_OFFSETS[v1];
            edgeVertices[e] = [
                baseX + (ox0 + t * (ox1 - ox0)) * dx,
                baseY + (oy0 + t * (oy1 - oy0)) * dy,
                baseZ + (oz0 + t * (oz1 - oz0)) * dz
            ];
        }
    }

    const triList = TRI_TABLE[cubeIndex];
    for (let i = 0; i < triList.length; i += 3) {
        const v0 = edgeVertices[triList[i]];
        const v1 = edgeVertices[triList[i + 1]];
        const v2 = edgeVertices[triList[i + 2]];
        if (v0 && v1 && v2) output.push(...v0, ...v1, ...v2);
    }
}

function separateComponents(triangles) {
    if (triangles.length < 9) return [triangles];
    const n = Math.floor(triangles.length / 9);
    const parent = new Int32Array(n);
    for (let i = 0; i < n; i++) parent[i] = i;

    function find(x) { if (parent[x] !== x) parent[x] = find(parent[x]); return parent[x]; }
    function union(a, b) { const ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; }

    const vertexToTris = new Map();
    const eps = 0.001;
    function key(i) {
        return `${Math.round(triangles[i] / eps)},${Math.round(triangles[i + 1] / eps)},${Math.round(triangles[i + 2] / eps)}`;
    }

    for (let t = 0; t < n; t++) {
        for (let v = 0; v < 3; v++) {
            const k = key(t * 9 + v * 3);
            if (!vertexToTris.has(k)) vertexToTris.set(k, []);
            vertexToTris.get(k).push(t);
        }
    }

    for (const tris of vertexToTris.values()) {
        for (let i = 1; i < tris.length; i++) union(tris[0], tris[i]);
    }

    const groups = new Map();
    for (let t = 0; t < n; t++) {
        const root = find(t);
        if (!groups.has(root)) groups.set(root, []);
        const arr = groups.get(root);
        for (let i = 0; i < 9; i++) arr.push(triangles[t * 9 + i]);
    }

    return Array.from(groups.values());
}

// =============== Worker 消息处理 ===============

self.onmessage = function (e) {
    const { type, taskId, data } = e.data;

    if (type === 'compute-isosurface') {
        const { orbitalParams, coeffs, bound, resolution, isovalue, atomType, slaterBasis, color } = data;

        // 波函数计算器（使用传入的 slaterBasis 数据）
        function calcPsi(x, y, z) {
            const r = Math.sqrt(x * x + y * y + z * z);
            if (r < 1e-10) return 0;
            const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
            const phi = Math.atan2(y, x);

            let psi = 0;
            for (let j = 0; j < orbitalParams.length; j++) {
                const op = orbitalParams[j];
                const R = radialR(op.n, op.l, r, 1, A0, atomType, slaterBasis);
                const Y = realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                psi += coeffs[j] * R * Y;
            }
            return psi;
        }

        // 进度回调
        function onProgress(progress) {
            self.postMessage({ type: 'progress', taskId, progress });
        }

        try {
            const bounds = { min: [-bound, -bound, -bound], max: [bound, bound, bound] };
            const result = marchingCubes(calcPsi, bounds, resolution, isovalue, onProgress);

            // 分离连通分量
            const positiveComponents = separateComponents(result.positive);
            const negativeComponents = separateComponents(result.negative);

            self.postMessage({
                type: 'isosurface-result',
                taskId,
                result: {
                    positive: positiveComponents,
                    negative: negativeComponents,
                    color
                }
            });
        } catch (error) {
            self.postMessage({
                type: 'error',
                taskId,
                error: error.message
            });
        }
    }
};

// Worker 就绪
self.postMessage({ type: 'ready' });
