/**
 * 共享物理核心模块
 * 
 * 此文件包含 physics.js 和 sampling-worker.js 共享的核心物理计算函数
 * 通过 importScripts 在 Worker 中加载，或通过 script 标签在主线程加载
 * 
 * 【重要】修改此文件时，两边的行为会同时改变，确保测试充分
 */

(function (root) {
    'use strict';

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
        const val = realYlm_value(l, m, type, theta, phi);
        return val * val;
    }

    // 实球谐函数（值）
    // Chemistry Convention: Positive lobes align with positive axes
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
            const nFact2 = factorial(2 * nStar);
            const norm = Math.pow(2 * zeta, nStar + 0.5) / Math.sqrt(nFact2);
            const val = coeff * norm * Math.pow(r, nStar - 1) * Math.exp(-zeta * r);
            sum += val;
        }
        return sum;
    }

    // 获取 SlaterBasis（兼容主线程和 Worker）
    function getSlaterBasis() {
        if (typeof window !== 'undefined' && window.SlaterBasis) {
            return window.SlaterBasis;
        }
        if (typeof self !== 'undefined' && self.SlaterBasis) {
            return self.SlaterBasis;
        }
        return null;
    }

    // 归一化径向函数 R_nl(r)
    function radialR(n, l, r, Z, a0, atomType) {
        Z = Z !== undefined ? Z : 1;
        a0 = a0 !== undefined ? a0 : A0;
        atomType = atomType || 'H';

        // Strategy A: Slater Type Orbitals (STO)
        const SlaterBasis = getSlaterBasis();
        if (atomType && atomType !== 'H' && SlaterBasis) {
            const orbitalKey = getOrbitalKey(n, l);
            const atomData = SlaterBasis[atomType];
            if (atomData && atomData.orbitals && atomData.orbitals[orbitalKey]) {
                // 【关键修复】Clementi-Roetti STO 基组的相位约定与氢原子解析解相反
                return -slaterRadialR(atomData.orbitals[orbitalKey], r);
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

    // 径向概率密度函数
    function radialPDF(n, l, r, Z, a0, atomType) {
        Z = Z !== undefined ? Z : 1;
        a0 = a0 !== undefined ? a0 : A0;
        atomType = atomType || 'H';
        const R = radialR(n, l, r, Z, a0, atomType);
        return r * r * R * R;
    }

    // 导出到全局
    const PhysicsCore = {
        A0,
        PI,
        factorial,
        binomialInt,
        generalizedLaguerre,
        associatedLegendre,
        Ylm_complex,
        Ylm_abs2,
        realYlm_abs2,
        realYlm_value,
        getOrbitalKey,
        slaterRadialR,
        radialR,
        radialPDF
    };

    // 兼容主线程和 Worker
    if (typeof module !== 'undefined' && module.exports) {
        module.exports = PhysicsCore;
    } else if (typeof define === 'function' && define.amd) {
        define(function () { return PhysicsCore; });
    } else {
        root.PhysicsCore = PhysicsCore;
    }

})(typeof self !== 'undefined' ? self : this);
