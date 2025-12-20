// Hydrogen atom physics model (real orbitals). Attaches to window.Hydrogen
(function () {
  const A0 = 1; // Bohr radius unit scale for visualization
  const PI = Math.PI;
  const TWO_PI = 2 * Math.PI;

  const FACT = (() => { const f = [1]; for (let i = 1; i <= 64; i++) f[i] = f[i - 1] * i; return f; })();
  function factorial(n) { return FACT[n] ?? Infinity; }
  function binomialInt(n, k) { if (k < 0 || k > n) return 0; return factorial(n) / (factorial(k) * factorial(n - k)); }

  // Generalized Laguerre L_k^{alpha}(x), alpha integer >=0, small k
  function generalizedLaguerre(k, alpha, x) {
    let sum = 0;
    for (let i = 0; i <= k; i++) {
      const c = ((i % 2) ? -1 : 1) * binomialInt(k + alpha, k - i) / factorial(i);
      sum += c * Math.pow(x, i);
    }
    return sum;
  }

  // Associated Legendre P_l^m(x), Condon–Shortley
  function associatedLegendre(l, m, x) {
    const mm = Math.abs(m);
    if (l < mm) return 0;
    let pmm = 1.0;
    if (mm > 0) {
      const s = Math.sqrt(Math.max(0, 1 - x * x));
      let fact = 1.0;
      for (let i = 1; i <= mm; i++) { pmm *= -fact * s; fact += 2.0; }
    }
    if (l === mm) return pmm;
    let pmmp1 = x * (2 * mm + 1) * pmm;
    if (l === mm + 1) return pmmp1;
    let pll = 0;
    for (let ll = mm + 2; ll <= l; ll++) {
      pll = ((2 * ll - 1) * x * pmmp1 - (ll + mm - 1) * pmm) / (ll - mm);
      pmm = pmmp1; pmmp1 = pll;
    }
    return pll;
  }

  // Complex spherical harmonic Y_l^m(θ,φ)
  // MOVED PHASE CORRECTION TO realYlm_value TO MATCH CHEMISTRY CONVENTION
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
    // Standard definition without CS phase for internal calculation, handled in real combo
    // We keep basic math correct here
    if (m > 0) {
      return { re: base * cos_mphi, im: base * sin_mphi };
    } else {
      // m < 0
      const sign = (mm % 2) ? -1 : 1;
      // Note: We will handle the chemistry sign convention in realYlm_value
      return { re: sign * base * cos_mphi, im: -sign * base * sin_mphi };
    }
  }

  function Ylm_abs2(l, m, theta, phi) {
    const y = Ylm_complex(l, m, theta, phi);
    return y.re * y.re + y.im * y.im;
  }

  // Real spherical harmonics (linear combos of complex):
  // Chemistry convention: Positive lobes align with positive axes
  // m=0: Y_l0 real
  // m>0: cos-type: √2 Re Y_l^m
  // m>0: sin-type: √2 Im Y_l^m
  // Removed Condon-Shortley phase factor (-1)^m to match Chemistry convention
  function realYlm_abs2(l, m, type, theta, phi) {
    const val = realYlm_value(l, m, type, theta, phi); // Re-use value function to ensure consistency
    return val * val;
  }

  // Real spherical harmonics value (not squared), normalized so ∫|Y|^2 dΩ = 1
  // Real spherical harmonics value (not squared), normalized
  // Chemistry Convention Enforced here:
  // p_x (l=1, m=1, c) ~ x
  // p_y (l=1, m=1, s) ~ y
  // 
  // 【关键修复】抵消 Condon-Shortley 相位
  // associatedLegendre 包含 (-1)^m 因子，对于奇数 m（如 p 轨道 m=1）会产生符号翻转
  // 我们需要抵消这个相位，确保 p_x 在 +x 轴方向为正值
  function realYlm_value(l, m, type, theta, phi) {
    if (m === 0) {
      const y = Ylm_complex(l, 0, theta, phi);
      return y.re;
    }
    const mm = Math.abs(m);
    const y = Ylm_complex(l, mm, theta, phi);

    // 【修复】抵消 CS 相位：对于奇数 m，翻转符号
    // 这确保了化学约定：正瓣指向正轴
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

  /**
   * 基于 Slater 规则计算有效核电荷 Z_eff
   * 
   * Slater 屏蔽规则（简化版）：
   * - 同组 (n, l) 电子：屏蔽常数 0.35（1s 为 0.30）
   * - (n-1) 壳层：屏蔽常数 0.85（s/p→d）或 1.00（d/f）
   * - 更内层壳层：屏蔽常数 1.00
   * 
   * @param {number} Z - 核电荷数
   * @param {number} n - 主量子数
   * @param {number} l - 角量子数
   * @returns {number} 有效核电荷 Z_eff
   */
  function getEffectiveZ(Z, n, l) {
    // 简化的电子构型屏蔽估计
    // 基于 Slater 规则的近似实现

    // 轨道分组：(1s), (2s,2p), (3s,3p), (3d), (4s,4p), (4d), (4f), ...
    // 屏蔽常数：
    // - 同组电子：0.35（1s 为 0.30）
    // - 相邻 (n-1) 的 s/p：0.85
    // - 更内层或 d/f：1.00

    let shielding = 0;

    // 1s 轨道
    if (n === 1) {
      // 另一个 1s 电子的屏蔽
      shielding = (Z > 1) ? 0.30 : 0;
    }
    // 2s, 2p 轨道
    else if (n === 2 && l <= 1) {
      // 1s 完全屏蔽
      shielding = 2 * 0.85;
      // 同层电子（2s, 2p 共 8 个，减去自己）
      const sameShell = Math.min(Z - 1, 8) - 2; // 减去 1s 的 2 个
      if (sameShell > 0) shielding += sameShell * 0.35;
    }
    // 3s, 3p 轨道
    else if (n === 3 && l <= 1) {
      // 1s: 1.00, 2s2p: 0.85
      shielding = 2 * 1.00 + 8 * 0.85;
      const sameShell = Math.max(0, Math.min(Z - 10, 8) - 1);
      shielding += sameShell * 0.35;
    }
    // 3d 轨道
    else if (n === 3 && l === 2) {
      // 1s, 2s2p: 1.00, 3s3p: 1.00（对 d 轨道屏蔽更强）
      shielding = 2 * 1.00 + 8 * 1.00 + 8 * 1.00;
      const sameShell = Math.max(0, Math.min(Z - 18, 10) - 1);
      shielding += sameShell * 0.35;
    }
    // 4s, 4p 轨道
    else if (n === 4 && l <= 1) {
      // 内层完全屏蔽
      shielding = 2 * 1.00 + 8 * 1.00 + 8 * 0.85 + 10 * 1.00;
      const sameShell = Math.max(0, Math.min(Z - 28, 8) - 1);
      shielding += sameShell * 0.35;
    }
    // 更高轨道的简化处理
    else {
      // 简单估计：内层电子数减去穿透效应补偿
      const innerElectrons = 2 * n * (n - 1); // 粗略估计
      shielding = Math.min(innerElectrons, Z - 1) * 0.85;
    }

    const Zeff = Math.max(1, Z - shielding);
    return Zeff;
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

  // Normalized radial R_nl(r) with ∫|R|^2 r^2 dr = 1
  // Modified to support STO strategy for non-Hydrogen atoms
  function radialR(n, l, r, Z = 1, a0 = A0, atomType = 'H') {
    // 策略 A: Slater Type Orbitals (STO)
    if (atomType && atomType !== 'H' && window.SlaterBasis) {
      const orbitalKey = getOrbitalKey(n, l);
      const atomData = window.SlaterBasis[atomType];
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
        // 注意：我们忽略极小数值的符号震荡，只看主导符号
        if (stoSign * analyticalSign < 0) {
          return -val;
        }
        return val;
      }
      // If orbital not defined for this atom (e.g. C-3s), fall back to Hydrogen-like or 0?
      // Fallback to Hydrogen-like with Z_eff might be better, or just pure Hydrogen Z=1 for now.
      // Let's fall back to Hydrogen Z=1 to avoid breaking.
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

  function radialPDF(n, l, r, Z = 1, a0 = A0, atomType = 'H') {
    const R = radialR(n, l, r, Z, a0, atomType);
    return r * r * (R * R);
  }

  /**
   * 计算杂化轨道的径向概率密度函数
   * 
   * 【物理原理】
   * 杂化轨道波函数：Ψ_hybrid = Σ cᵢ Rᵢ(r) Yᵢ(θ,φ)
   * 由于球谐函数的正交性，对角向积分后交叉项消失
   * 
   * P_hybrid(r) = r² × Σ |cᵢ|² × |Rᵢ(r)|²
   * 
   * 注意：这与简单地平均各轨道的 radialPDF 不同！
   * 
   * @param {Array} paramsList - 轨道参数列表
   * @param {number} r - 径向距离
   * @param {number} Z - 核电荷数（默认1）
   * @param {number} a0 - 玻尔半径（默认1）
   * @param {string} atomType - 原子类型
   * @returns {number} - 杂化轨道的径向概率密度 P(r)
   */
  function hybridRadialPDF(paramsList, r, Z = 1, a0 = A0, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return 0;

    // 【数学严谨性修复】
    // 原公式 P = r² * Σ |c_i|² * |R_i|² 仅在所有轨道角向正交时成立
    // 对于角向相同轨道（如 1s + 2s 或 2pz + 3pz）的混合，必须包含交叉项（干涉项）
    // 通用公式：P(r) = r² * ∫ |Σ c_i R_i Y_i|² dΩ
    //                = r² * Σ_i Σ_j c_i c_j R_i R_j * ⟨Y_i|Y_j⟩
    // 由于实球谐函数正交性，⟨Y_i|Y_j⟩ = δ_{ang_i, ang_j}

    const numOrbitals = paramsList.length;
    const defaultCoeff = 1.0 / Math.sqrt(numOrbitals);

    // 预计算所有 R 值，避免重复计算
    const R_values = new Float32Array(numOrbitals);
    for (let i = 0; i < numOrbitals; i++) {
      const p = paramsList[i];
      R_values[i] = radialR(p.n, p.l, r, Z, a0, atomType);
    }

    let densitySum = 0;

    for (let i = 0; i < numOrbitals; i++) {
      const p1 = paramsList[i];
      const c1 = p1.coefficient !== undefined ? p1.coefficient : defaultCoeff;

      for (let j = 0; j < numOrbitals; j++) {
        const p2 = paramsList[j];
        const c2 = p2.coefficient !== undefined ? p2.coefficient : defaultCoeff;

        // 检查角向部分是否正交
        // 实球谐函数正交条件：l, m, type 均相同
        // 注意：type (t) 可能是 undefined (m=0时), 需要处理
        const theta1 = p1.angKey ? p1.angKey.t : '';
        const theta2 = p2.angKey ? p2.angKey.t : '';

        const isOrthogonal = (p1.angKey.l !== p2.angKey.l) ||
          (p1.angKey.m !== p2.angKey.m) ||
          (theta1 !== theta2);

        if (!isOrthogonal) {
          // 角向相同，积分结果为 1，贡献交叉项
          densitySum += c1 * c2 * R_values[i] * R_values[j];
        }
        // 如果角向不同，积分结果为 0，交叉项消失
      }
    }

    return r * r * densitySum;
  }

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
   * 【重要说明】
   * 此函数计算的是"一条"杂化轨道的密度（使用默认等权系数）
   * 对于 sp³ 杂化，共有 4 条杂化轨道，此函数只返回其中一条的密度
   * 
   * 如需计算所有杂化轨道的总密度，请使用 allHybridOrbitalsDensity3D
   * 如需计算第 i 条杂化轨道的密度，请使用 singleHybridDensity3D
   * 
   * @param {Array} paramsList - 轨道参数列表
   * @param {number} r - 径向距离
   * @param {number} theta - 极角
   * @param {number} phi - 方位角
   * @param {string} atomType - 原子类型
   * @returns {number} - 概率密度 |Ψ_hybrid|²
   */
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
   * 计算所有杂化轨道的总概率密度
   * 
   * 【物理原理】
   * N 个原子轨道产生 N 条杂化轨道
   * "全部"模式下显示所有 N 条杂化轨道的叠加：Σ|ψ_hybrid_i|²
   * 
   * 例如 sp³ 杂化：4 个原子轨道 → 4 条杂化轨道，指向四面体的四个顶点
   * 
   * @param {Array} paramsList - 原子轨道参数列表
   * @param {number} r - 径向距离
   * @param {number} theta - 极角
   * @param {number} phi - 方位角
   * @param {string} atomType - 原子类型
   * @returns {number} - 所有杂化轨道的总概率密度
   */
  function allHybridOrbitalsDensity3D(paramsList, r, theta, phi, Z = 1, a0 = A0, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return 0;

    const numOrbitals = paramsList.length;
    const coeffMatrix = getHybridCoefficients(numOrbitals, paramsList);

    let totalDensity = 0;

    // 遍历每条杂化轨道
    for (let hybridIdx = 0; hybridIdx < numOrbitals; hybridIdx++) {
      const coeffs = coeffMatrix[hybridIdx];
      let psiReal = 0;

      // 计算该杂化轨道的波函数
      for (let i = 0; i < numOrbitals; i++) {
        const p = paramsList[i];
        const coeff = coeffs[i];
        const R = radialR(p.n, p.l, r, Z, a0, atomType);
        const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
        psiReal += coeff * R * Y;
      }

      // 累加该杂化轨道的概率密度
      totalDensity += psiReal * psiReal;
    }

    return totalDensity;
  }

  /**
   * 对轨道列表进行排序，以符合杂化系数矩阵的预期顺序
   * 顺序：s, p (x, y, z), d (z2, x2-y2, ...)
   */
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

      // 对于 d (l=2): dz2(m=0), dx2-y2(m=2,c), ...
      // sp3d2 需要: dx2-y2, dz2
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

  /**
   * 获取杂化轨道的系数矩阵
   * 
   * 【核心算法：Thomson 方向优化】
   * 使用 Thomson 问题（球面电荷排斥）找到 N 个最优分布的方向，
   * 然后用几何优化计算杂化系数，使杂化轨道指向这些方向。
   * 
   * @param {number} numOrbitals - 参与杂化的轨道数量
   * @param {Array} orbitalParams - 轨道参数列表（用于确定角量子数）
   * @returns {Array} - 系数矩阵，每行是一个杂化轨道的系数
   */

  // Thomson 方向缓存
  const _thomsonCache = {};

  /**
   * 计算 N 个球面点的静电势能
   */
  function thomsonEnergy(points) {
    const n = points.length;
    let energy = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const dx = points[i][0] - points[j][0];
        const dy = points[i][1] - points[j][1];
        const dz = points[i][2] - points[j][2];
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        energy += 1 / (dist + 1e-10);
      }
    }
    return energy;
  }

  /**
   * 归一化点到球面
   */
  function normalizeToSphere(points) {
    return points.map(p => {
      const norm = Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
      if (norm < 1e-10) return [0, 0, 1];
      return [p[0] / norm, p[1] / norm, p[2] / norm];
    });
  }

  /**
   * 用梯度下降优化 Thomson 问题
   */
  function optimizeThomson(n, maxIter = 200, lr = 0.1) {
    if (n <= 1) return [[0, 0, 1]];

    // 缓存检查
    if (_thomsonCache[n]) {
      return JSON.parse(JSON.stringify(_thomsonCache[n]));
    }

    // 【关键修复】使用确定性初始化（黄金角度分布）
    // 确保每次结果相同，且与 Worker 中的算法一致
    let points = [];
    const goldenRatio = (1 + Math.sqrt(5)) / 2;
    for (let i = 0; i < n; i++) {
      const y = 1 - (i / (n - 1)) * 2;
      const radius = Math.sqrt(1 - y * y);
      const t = 2 * Math.PI * i / goldenRatio;
      points.push([Math.cos(t) * radius, y, Math.sin(t) * radius]);
    }

    // 梯度下降
    for (let iter = 0; iter < maxIter; iter++) {
      // 计算梯度
      const grad = points.map(() => [0, 0, 0]);

      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          const dx = points[i][0] - points[j][0];
          const dy = points[i][1] - points[j][1];
          const dz = points[i][2] - points[j][2];
          const dist2 = dx * dx + dy * dy + dz * dz;
          const dist = Math.sqrt(dist2);
          const dist3 = dist2 * dist + 1e-10;

          // 梯度: d(1/r)/dr = -1/r² * (r_i - r_j)/|r_i - r_j|
          const factor = 1 / dist3;
          grad[i][0] -= factor * dx;
          grad[i][1] -= factor * dy;
          grad[i][2] -= factor * dz;
          grad[j][0] += factor * dx;
          grad[j][1] += factor * dy;
          grad[j][2] += factor * dz;
        }
      }

      // 更新位置
      for (let i = 0; i < n; i++) {
        points[i][0] -= lr * grad[i][0];
        points[i][1] -= lr * grad[i][1];
        points[i][2] -= lr * grad[i][2];
      }

      // 投影回球面
      points = normalizeToSphere(points);

      // 自适应学习率
      if (iter % 50 === 0) lr *= 0.8;
    }

    // 缓存结果
    _thomsonCache[n] = JSON.parse(JSON.stringify(points));
    return points;
  }

  /**
   * 根据轨道组合生成约束方向
   * 这些方向在轨道的"可达角向空间"内，确保物理正确性
   */
  function generateConstrainedDirections(orbitalParams) {
    const n = orbitalParams.length;

    // 分析轨道组合
    const hasS = orbitalParams.some(p => p.l === 0);
    const pOrbitals = orbitalParams.filter(p => p.l === 1);
    const dOrbitals = orbitalParams.filter(p => p.l === 2);

    const hasPx = pOrbitals.some(p => p.angKey.m === 1 && p.angKey.t === 'c');
    const hasPy = pOrbitals.some(p => p.angKey.m === 1 && p.angKey.t === 's');
    const hasPz = pOrbitals.some(p => p.angKey.m === 0);

    const hasDz2 = dOrbitals.some(p => p.angKey.m === 0);
    const hasDx2y2 = dOrbitals.some(p => p.angKey.m === 2 && p.angKey.t === 'c');

    // sp: 线性（z轴或x轴）
    if (n === 2 && hasS && pOrbitals.length === 1) {
      if (hasPz) return [[0, 0, 1], [0, 0, -1]];
      if (hasPx) return [[1, 0, 0], [-1, 0, 0]];
      if (hasPy) return [[0, 1, 0], [0, -1, 0]];
    }

    // sp2: 平面三角形
    if (n === 3 && hasS && pOrbitals.length === 2) {
      if (hasPx && hasPy) {
        // xy平面
        return [[1, 0, 0], [-0.5, Math.sqrt(3) / 2, 0], [-0.5, -Math.sqrt(3) / 2, 0]];
      }
      if (hasPx && hasPz) {
        // xz平面
        return [[1, 0, 0], [-0.5, 0, Math.sqrt(3) / 2], [-0.5, 0, -Math.sqrt(3) / 2]];
      }
      if (hasPy && hasPz) {
        // yz平面
        return [[0, 1, 0], [0, -0.5, Math.sqrt(3) / 2], [0, -0.5, -Math.sqrt(3) / 2]];
      }
    }

    // sp3: 正四面体
    if (n === 4 && hasS && pOrbitals.length === 3 && hasPx && hasPy && hasPz && dOrbitals.length === 0) {
      const a = 1 / Math.sqrt(3);
      return [[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]];
    }

    // sp3d: 三角双锥
    if (n === 5 && hasS && hasPx && hasPy && hasPz && hasDz2 && !hasDx2y2) {
      return [
        [0, 0, 1],                          // 轴向 +z
        [0, 0, -1],                         // 轴向 -z
        [1, 0, 0],                          // 赤道
        [-0.5, Math.sqrt(3) / 2, 0],          // 赤道 120°
        [-0.5, -Math.sqrt(3) / 2, 0]          // 赤道 240°
      ];
    }

    // sp3d2: 正八面体
    if (n === 6 && hasS && hasPx && hasPy && hasPz && hasDz2 && hasDx2y2) {
      return [
        [1, 0, 0], [-1, 0, 0],
        [0, 1, 0], [0, -1, 0],
        [0, 0, 1], [0, 0, -1]
      ];
    }

    // 其他情况：使用通用Thomson优化（适用于全3D空间可达的组合）
    return optimizeThomson(n);
  }

  /**
   * 计算方向矩阵：每行是一个方向，每列是一个轨道的角向值
   */
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
        const { l, m, t } = orbitalParams[j].angKey || { l: 0, m: 0, t: 'c' };

        // 计算实球谐函数在该方向的值
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

  /**
   * 简单的 Jacobi SVD 算法 (A = U S V^T)
   * 用于小矩阵 (N <= 10) 的奇异值分解
   * 返回 { U, S, V }
   */
  function jacobiSVD(A) {
    const m = A.length;
    const n = A[0].length;
    const minDim = Math.min(m, n);

    // 初始化
    let U = A.map(row => [...row]); // 复制 A
    let V = [];
    for (let i = 0; i < n; i++) {
      const row = new Array(n).fill(0);
      row[i] = 1;
      V.push(row);
    }

    const EPSILON = 1e-15;
    const MAX_ITER = 50;

    // Jacobi 迭代
    for (let iter = 0; iter < MAX_ITER; iter++) {
      let maxError = 0;

      for (let i = 0; i < n - 1; i++) {
        for (let j = i + 1; j < n; j++) {
          // 计算 2x2 子矩阵 J^T A^T A J
          let alpha = 0, beta = 0, gamma = 0;

          for (let k = 0; k < m; k++) {
            alpha += U[k][i] * U[k][i];
            beta += U[k][j] * U[k][j];
            gamma += U[k][i] * U[k][j];
          }

          maxError = Math.max(maxError, Math.abs(gamma) / Math.sqrt(alpha * beta + EPSILON));

          if (Math.abs(gamma) < EPSILON) continue;

          // 计算旋转角度
          let zeta = (beta - alpha) / (2 * gamma);
          let t = Math.sign(zeta) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
          let c = 1 / Math.sqrt(1 + t * t);
          let s = c * t;

          // 更新 U
          for (let k = 0; k < m; k++) {
            let t1 = U[k][i];
            let t2 = U[k][j];
            U[k][i] = c * t1 - s * t2;
            U[k][j] = s * t1 + c * t2;
          }

          // 更新 V
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

    // 提取奇异值 S
    const S = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      let sum = 0;
      for (let k = 0; k < m; k++) sum += U[k][i] * U[k][i];
      S[i] = Math.sqrt(sum);
      // 归一化 U 的列
      if (S[i] > EPSILON) {
        for (let k = 0; k < m; k++) U[k][i] /= S[i];
      }
    }

    return { U, S, V };
  }

  /**
   * 矩阵乘法 C = A * B
   */
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

  /**
   * 矩阵转置
   */
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

  /**
   * 几何优化：根据方向计算杂化系数
   * 使用 SVD 求解正交 Procrustes 问题：寻找正交矩阵 Q 使得 ||A - Q|| 最小
   * Q = U * V^T
   */
  function computeHybridCoefficients(directions, orbitalParams) {
    const A = buildDirectionMatrix(directions, orbitalParams);

    // SVD 分解 A = U S V^T
    const { U, V } = jacobiSVD(A);

    // Q = U * V^T
    // 注意：jacobiSVD 返回的 V 是 V 本身，不是 V^T
    const VT = matTranspose(V);
    const Q = matMul(U, VT);

    return Q;
  }

  function getHybridCoefficients(numOrbitals, orbitalParams) {
    // 强制要求必须提供轨道参数，绝不允许 Fallback
    if (!orbitalParams || orbitalParams.length !== numOrbitals) {
      console.error("getHybridCoefficients: Missing or invalid orbitalParams!", numOrbitals, orbitalParams);
      throw new Error(`getHybridCoefficients requires orbitalParams matching numOrbitals (${numOrbitals})`);
    }

    // 【核心修复】使用约束方向，确保在轨道的可达角向空间内
    const directions = generateConstrainedDirections(orbitalParams);

    // 使用几何优化 (SVD Procrustes) 计算精确系数
    return computeHybridCoefficients(directions, orbitalParams);
  }

  /**
   * 计算单个杂化轨道的概率密度
   * 
   * @param {Array} paramsList - 原子轨道参数列表
   * @param {number} hybridIndex - 杂化轨道索引 (0, 1, 2, ...)
   * @param {number} r - 径向距离
   * @param {number} theta - 极角
   * @param {number} phi - 方位角
   * @param {string} atomType - 原子类型
   * @returns {number} - 概率密度
   */
  function singleHybridDensity3D(paramsList, hybridIndex, r, theta, phi, Z = 1, a0 = A0, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return 0;

    const numOrbitals = paramsList.length;
    const coeffMatrix = getHybridCoefficients(numOrbitals, paramsList);

    // 确保索引有效
    const idx = hybridIndex % numOrbitals;
    const coeffs = coeffMatrix[idx];

    let psiReal = 0;

    for (let i = 0; i < numOrbitals; i++) {
      const p = paramsList[i];
      const coeff = coeffs[i];

      const R = radialR(p.n, p.l, r, Z, a0, atomType);
      const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);

      psiReal += coeff * R * Y;
    }

    return psiReal * psiReal;
  }

  /**
   * 计算单个杂化轨道的波函数值（用于相位）
   */
  function singleHybridWavefunction(paramsList, hybridIndex, r, theta, phi, Z = 1, a0 = A0, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return 0;

    const numOrbitals = paramsList.length;
    const coeffMatrix = getHybridCoefficients(numOrbitals, paramsList);
    const idx = hybridIndex % numOrbitals;
    const coeffs = coeffMatrix[idx];

    let psiReal = 0;

    for (let i = 0; i < numOrbitals; i++) {
      const p = paramsList[i];
      const coeff = coeffs[i];

      const R = radialR(p.n, p.l, r, Z, a0, atomType);
      const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);

      psiReal += coeff * R * Y;
    }

    return psiReal;
  }

  /**
   * 计算杂化轨道的推荐采样边界
   * 取所有轨道中最大的 rMax
   */
  function hybridRecommendRmax(paramsList, a0 = A0) {
    if (!paramsList || paramsList.length === 0) return 40;

    let maxRmax = 0;
    for (const p of paramsList) {
      const rmax = recommendRmax(p.n, a0);
      if (rmax > maxRmax) maxRmax = rmax;
    }
    return maxRmax;
  }

  /**
   * 估算杂化轨道的最大概率密度（用于拒绝采样）
   * 通过数值搜索找到近似最大值
   */
  function hybridEstimateMaxDensity(paramsList, rMax, numSamples = 1000) {
    if (!paramsList || paramsList.length === 0) return 1;

    let maxDensity = 0;

    // 在球体内随机采样点来估计最大密度
    for (let i = 0; i < numSamples; i++) {
      // 使用更密集的内核采样（概率密度通常在内核区域最高）
      const r = Math.random() * rMax * Math.pow(Math.random(), 0.5);
      const cosTheta = 2 * Math.random() - 1;
      const theta = Math.acos(cosTheta);
      const phi = 2 * PI * Math.random();

      const density = hybridDensity3D(paramsList, r, theta, phi);
      if (density > maxDensity) {
        maxDensity = density;
      }
    }

    // 添加安全边际
    return maxDensity * 1.5;
  }

  function recommendRmax(n, a0 = A0) { return 15 * n * n * a0; }

  function radialGrid(n, l, rmax, num = 512, Z = 1, a0 = A0) {
    const rs = new Float32Array(num);
    const ps = new Float32Array(num);
    const dr = rmax / (num - 1);
    for (let i = 0; i < num; i++) {
      const r = i * dr; rs[i] = r; ps[i] = radialPDF(n, l, r, Z, a0);
    }
    return { r: rs, Pr: ps };
  }

  function histogramRadialFromSamples(rArray, bins = 160, rmax = null, normalize = true, smooth = true) {
    const N = rArray.length; if (N === 0) return { edges: [], counts: [] };
    // 【性能修复】使用循环替代Math.max(...array)，避免大数组栈溢出
    let maxr;
    if (rmax !== null) {
      maxr = rmax;
    } else {
      maxr = 0;
      for (let i = 0; i < N; i++) {
        if (rArray[i] > maxr) maxr = rArray[i];
      }
    }

    // 确保最小范围和bins数
    const effectiveMaxr = Math.max(maxr, 0.1);
    const effectiveBins = Math.max(bins, 10);

    const edges = new Float32Array(effectiveBins + 1);
    let counts = new Float32Array(effectiveBins);
    const dr = effectiveMaxr / effectiveBins;

    for (let i = 0; i <= effectiveBins; i++) edges[i] = i * dr;
    for (let i = 0; i < N; i++) {
      const r = rArray[i];
      if (r < 0 || r > effectiveMaxr) continue;
      const b = Math.min(effectiveBins - 1, Math.floor(r / dr));
      counts[b] += 1;
    }

    // 【物理优化】添加高斯平滑以减少采样毛刺
    // 使用自适应窗口大小：根据样本量和bin数量动态调整
    if (smooth && N > 100) {
      const smoothedCounts = new Float32Array(effectiveBins);

      // 计算每个bin的平均样本数
      const samplesPerBin = N / effectiveBins;

      // 自适应平滑窗口：
      // - 样本稀疏（<50/bin）时使用较大窗口（2-3）
      // - 样本中等（50-200/bin）时使用小窗口（1）
      // - 样本密集（>200/bin）时不平滑（0）
      let windowSize;
      if (samplesPerBin < 50) {
        windowSize = 2;
      } else if (samplesPerBin < 200) {
        windowSize = 1;
      } else {
        windowSize = 0;
      }

      if (windowSize > 0) {
        for (let i = 0; i < effectiveBins; i++) {
          let sum = 0;
          let weightSum = 0;
          for (let j = -windowSize; j <= windowSize; j++) {
            const idx = i + j;
            if (idx >= 0 && idx < effectiveBins) {
              // 高斯权重
              const sigma = windowSize * 0.6;
              const weight = Math.exp(-0.5 * (j / sigma) ** 2);
              sum += counts[idx] * weight;
              weightSum += weight;
            }
          }
          smoothedCounts[i] = sum / weightSum;
        }
        counts = smoothedCounts;
      }
    }

    if (normalize) {
      let area = 0;
      for (let i = 0; i < effectiveBins; i++) area += counts[i] * dr;
      if (area > 0) {
        for (let i = 0; i < effectiveBins; i++) counts[i] /= area;
      }
    }
    return { edges, counts, dr, rmax: effectiveMaxr };
  }

  function histogramThetaFromSamples(thetaArray, bins = 180, normalize = true) {
    const N = thetaArray.length; if (N === 0) return { edges: [], counts: [] };
    const edges = new Float32Array(bins + 1);
    const counts = new Float32Array(bins);
    const dth = Math.PI / bins;
    for (let i = 0; i <= bins; i++) edges[i] = i * dth;
    for (let i = 0; i < N; i++) {
      const t = thetaArray[i]; if (t < 0 || t > Math.PI) continue;
      const b = Math.min(bins - 1, Math.floor(t / dth)); counts[b] += 1;
    }
    if (normalize) { let area = 0; for (let i = 0; i < bins; i++) area += counts[i] * dth; if (area > 0) { for (let i = 0; i < bins; i++) counts[i] /= area; } }
    return { edges, counts, dθ: dth };
  }

  /**
   * 从采样数据生成 φ (方位角) 直方图
   * φ 范围是 [0, 2π]
   * @param {Array} phiArray - φ 值数组
   * @param {number} bins - 直方图的 bin 数量
   * @param {boolean} normalize - 是否归一化
   * @returns {Object} { edges, counts, dφ }
   */
  function histogramPhiFromSamples(phiArray, bins = 180, normalize = true) {
    const N = phiArray.length;
    if (N === 0) return { edges: [], counts: [] };

    const edges = new Float32Array(bins + 1);
    const counts = new Float32Array(bins);
    const dphi = TWO_PI / bins;

    // 生成边缘 [0, 2π]
    for (let i = 0; i <= bins; i++) {
      edges[i] = i * dphi;
    }

    // 计算直方图
    for (let i = 0; i < N; i++) {
      let p = phiArray[i];
      // 确保 phi 在 [0, 2π) 范围内
      if (p === undefined || isNaN(p)) continue;
      while (p < 0) p += TWO_PI;
      while (p >= TWO_PI) p -= TWO_PI;
      const b = Math.min(bins - 1, Math.floor(p / dphi));
      counts[b] += 1;
    }

    // 归一化
    if (normalize) {
      let area = 0;
      for (let i = 0; i < bins; i++) {
        area += counts[i] * dphi;
      }
      if (area > 0) {
        for (let i = 0; i < bins; i++) {
          counts[i] /= area;
        }
      }
    }

    return { edges, counts, dφ: dphi };
  }

  function orbitalKey(params) {
    return `${params.n}${params.l}${params.angKey.m}${params.angKey.t}`;
  }

  // 计算角向边缘概率密度 P(θ) = sin(θ) × ∫|Y|² dφ
  // 对于实球谐函数，无论 m=0 还是 m≠0，公式统一为 2π × N² × P_lm² × sin(θ)
  function angularPDF_Theta(l, m, theta) {
    const mm = Math.abs(m);
    const Plm = associatedLegendre(l, mm, Math.cos(theta));
    const N2 = ((2 * l + 1) / (4 * PI)) * (factorial(l - mm) / factorial(l + mm));
    return 2 * PI * N2 * Plm * Plm * Math.sin(theta);
  }

  /**
   * 计算方位角边缘概率密度 P(φ) = ∫|Y|² sin(θ) dθ
   * 对于实球谐函数：
   * - m=0 时：P(φ) = 1/(2π)，均匀分布
   * - m≠0 且 cos型：P(φ) = cos²(mφ)/π
   * - m≠0 且 sin型：P(φ) = sin²(mφ)/π
   * 
   * @param {number} m - 磁量子数
   * @param {string} t - 轨道类型：'c'=cos型, 's'=sin型, ''=m=0
   * @param {number} phi - 方位角 [0, 2π]
   * @returns {number} 概率密度
   */
  function angularPDF_Phi(m, t, phi) {
    const mm = Math.abs(m);

    if (mm === 0) {
      // m=0: 均匀分布
      return 1 / TWO_PI;
    }

    // m≠0: 依赖于轨道类型
    if (t === 'c' || t === 'cos') {
      // cos型：P(φ) = cos²(mφ)/π
      const cosMphi = Math.cos(mm * phi);
      return cosMphi * cosMphi / PI;
    } else {
      // sin型：P(φ) = sin²(mφ)/π
      const sinMphi = Math.sin(mm * phi);
      return sinMphi * sinMphi / PI;
    }
  }

  // ==================== 精确逆CDF采样 (Inverse CDF Sampling) ====================
  // 
  // 【物理准确性保证】
  // 利用解析式 P(r) = r² |R_nl(r)|² 构建精确的累积分布函数
  // 通过数值积分预计算 CDF，然后使用二分查找进行逆变换采样
  // 
  // 优势：
  // 1. 100% 接受率（无拒绝采样）
  // 2. 精确服从理论分布
  // 3. 没有提议分布匹配问题
  // ====================

  // CDF 表缓存
  const _cdfCache = {};

  /**
   * 构建径向分布的累积分布函数表
   * @param {number} n - 主量子数
   * @param {number} l - 角量子数
   * @param {number} numPoints - CDF表的点数（越多越精确）
   * @param {string} atomType - 原子类型
   * @returns {Object} { r: Float64Array, cdf: Float64Array, rMax: number }
   */
  function buildRadialCDF(n, l, numPoints = 2000, atomType = 'H') {
    const key = `${n}_${l}_${atomType}`;
    if (_cdfCache[key]) {
      return _cdfCache[key];
    }

    // 确定积分范围：到概率密度下降到峰值的 1e-8 为止
    const rMax = Math.max(4 * n * n * A0, 50 * A0);
    const dr = rMax / numPoints;

    const r = new Float64Array(numPoints + 1);
    const cdf = new Float64Array(numPoints + 1);

    // 使用梯形法则进行数值积分
    r[0] = 0;
    cdf[0] = 0;
    let cumulative = 0;

    for (let i = 1; i <= numPoints; i++) {
      r[i] = i * dr;
      const P_prev = radialPDF(n, l, r[i - 1], 1, 1, atomType);
      const P_curr = radialPDF(n, l, r[i], 1, 1, atomType);
      // 梯形法则：积分 ≈ (f(a) + f(b)) * h / 2
      cumulative += (P_prev + P_curr) * dr / 2;
      cdf[i] = cumulative;
    }

    // 归一化 CDF 到 [0, 1]
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
   * @param {number} n - 主量子数
   * @param {number} l - 角量子数
   * @param {string} atomType - 原子类型
   * @returns {number} 采样得到的 r 值
   */
  function sampleRadialExact(n, l, atomType = 'H') {
    const { r, cdf, numPoints } = buildRadialCDF(n, l, 2000, atomType);

    // 生成 [0, 1) 均匀随机数
    const u = Math.random();

    // 二分查找：找到 cdf[i] <= u < cdf[i+1] 的位置
    let lo = 0, hi = numPoints;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (cdf[mid] < u) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }

    // 线性插值得到更精确的 r 值
    const i = Math.max(0, lo - 1);
    const j = Math.min(numPoints, lo);

    if (i === j || cdf[j] === cdf[i]) {
      return r[i];
    }

    // 在区间 [r[i], r[j]] 内线性插值
    const t = (u - cdf[i]) / (cdf[j] - cdf[i]);
    return r[i] + t * (r[j] - r[i]);
  }

  // ==================== 重要性采样 (Importance Sampling) ====================
  // 
  // 【核心设计】多峰混合提议分布
  // 
  // 氢原子轨道的径向分布有 (n - l) 个峰值。为了准确采样所有峰，
  // 我们使用 **多峰混合 Gamma 分布** 作为提议分布：
  // 
  // q(r) = Σ w_k × Gamma(3, α_k)(r)
  // 
  // 其中每个成分覆盖一个概率峰，权重 w_k 近似于该峰的积分贡献。
  // ==================== 

  /**
   * 计算轨道的所有径向峰位置（数值方法）
   * 氢原子 (n,l) 轨道有 (n-l) 个径向峰
   */
  function findRadialPeaks(n, l, a0 = A0) {
    const peaks = [];
    const rmax = n * n * 3 * a0;  // 搜索范围
    const dr = 0.1 * a0;

    let prevVal = 0;
    let prevPrevVal = 0;

    for (let r = dr; r < rmax; r += dr) {
      const val = radialPDF(n, l, r);

      // 检测峰值：前一个点大于两侧
      if (r > 2 * dr && prevVal > prevPrevVal && prevVal > val) {
        peaks.push({ r: r - dr, pdf: prevVal });
      }

      prevPrevVal = prevVal;
      prevVal = val;
    }

    return peaks;
  }

  // 缓存峰位置以避免重复计算
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
   * 为每个峰创建一个 Gamma(3) 分布成分
   * 
   * 【物理优化】使用积分贡献估计权重，而非峰高度
   * 这样可以确保每个峰按其实际概率贡献获得采样机会
   */
  function getMultiPeakMixtureParams(n, l, a0 = A0) {
    const peaks = getCachedPeaks(n, l);

    if (peaks.length === 0) {
      // 回退：使用简单估计
      const r_peak = n * n * a0;
      return {
        components: [{ alpha: 2.0 / r_peak, weight: 1.0 }]
      };
    }

    // 【改进】计算每个峰的权重（基于估计的积分贡献）
    // 积分贡献 ≈ 峰高 × 峰宽，峰宽与 r_peak 成正比
    let totalWeight = 0;
    const components = [];

    for (const peak of peaks) {
      // 使用 pdf × r_peak 作为积分贡献的估计（更接近真实积分）
      // 这比单纯使用 sqrt(pdf) 更物理准确
      const estimatedIntegral = peak.pdf * peak.r;
      const weight = estimatedIntegral;
      totalWeight += weight;

      // Gamma(3) 分布在 r = 2/α 处达到峰值，所以 α = 2/r_peak
      const alpha = 2.0 / peak.r;
      components.push({ alpha, weight, r_peak: peak.r, pdf: peak.pdf });
    }

    // 归一化权重
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
      // Gamma(3, α) 的 PDF: α³/2 × r² × exp(-αr)
      const alpha = comp.alpha;
      const norm = alpha * alpha * alpha / 2.0;
      pdf += comp.weight * norm * r * r * Math.exp(-alpha * r);
    }

    return pdf;
  }

  /**
   * 从 Gamma(3) 分布采样（Erlang 方法）
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

    // 按权重随机选择一个成分
    const u = Math.random();
    let cumWeight = 0;

    for (const comp of params.components) {
      cumWeight += comp.weight;
      if (u < cumWeight) {
        return sampleGamma3(comp.alpha);
      }
    }

    // 回退（不应该到达这里）
    return sampleGamma3(params.components[params.components.length - 1].alpha);
  }

  /**
   * 从均匀球面分布采样角度 (θ, φ)
   */
  function sampleUniformSphere() {
    const phi = TWO_PI * Math.random();
    const cosTheta = 2 * Math.random() - 1;
    const theta = Math.acos(cosTheta);
    return { theta, phi };
  }

  // 缓存权重上界计算结果
  const _maxWeightCache = {};

  /**
   * 计算径向权重的理论上界
   * 【物理优化】使用数值搜索找到精确的最大权重
   * 这避免了保守估计导致的效率损失和激进估计导致的偏差
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
      // 在峰附近的范围内搜索
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

    // 添加 20% 安全边际，确保不会截断样本
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
   * 
   * 这确保最终分布精确等于 |ψ|² = R²(r) × |Y|²(θ,φ)
   */
  function importanceSample(n, l, angKey, samplingBoundary = Infinity, atomType = 'H') {
    // ==================== 第一步：径向采样（精确逆CDF） ====================
    // 直接从精确的 P(r) 分布采样，无需接受-拒绝
    let r = sampleRadialExact(n, l, atomType);

    // 边界检查
    if (r > samplingBoundary * 2) {
      return null;
    }

    // ==================== 第二步：角向采样 ====================
    const { theta, phi } = sampleUniformSphere();

    // 计算角向权重：w_angular = 4π |Y|²
    const Y2 = realYlm_abs2(angKey.l, angKey.m, angKey.t, theta, phi);
    const w_angular = 4 * PI * Y2;

    // 【物理修正】角向权重上界的精确计算
    // 对于实球谐函数：
    // - m=0: 最大值 = (2l+1)/(4π)，所以 4π|Y|² ≤ 2l+1
    // - m≠0: 由于 cos² 或 sin² 因子，最大值可能更大
    // 使用更保守的上界以确保物理准确性
    const maxAngularWeight = (angKey.m === 0) ? (2 * angKey.l + 1.2) : (2 * (2 * angKey.l + 1) + 0.5);

    // 角向接受-拒绝
    if (Math.random() * maxAngularWeight > w_angular) {
      return { x: 0, y: 0, z: 0, r, theta, phi, weight: w_angular, accepted: false };
    }

    // ==================== 采样成功 ====================
    const sinTheta = Math.sin(theta);
    const x = r * sinTheta * Math.cos(phi);
    const y = r * sinTheta * Math.sin(phi);
    const z = r * Math.cos(theta);

    // 【摩尔纹抑制】添加亚像素级的微小抖动
    const dither = 0.01;
    const dx = (Math.random() - 0.5) * dither;
    const dy = (Math.random() - 0.5) * dither;
    const dz = (Math.random() - 0.5) * dither;

    return { x: x + dx, y: y + dy, z: z + dz, r, theta, phi, weight: 1, accepted: true };
  }

  /**
   * 批量重要性采样
   */
  function batchImportanceSample(n, l, angKey, targetCount, samplingBoundary, maxAttempts) {
    const samples = [];
    let attempts = 0;

    while (samples.length < targetCount && attempts < maxAttempts) {
      attempts++;
      const result = importanceSample(n, l, angKey, samplingBoundary);
      if (result && result.accepted) {
        samples.push(result);
      }
    }

    return samples;
  }

  // 扩展轨道参数映射，添加 n=5 的所有轨道
  function orbitalParamsFromKey(key) {
    // Map UI key to {n,l,angKey:{l,m,t}}
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
      // seven real f
      case '4f_z3': return R(4, 3, 0, 'c');
      case '4f_xz2': return R(4, 3, 1, 'c');
      case '4f_yz2': return R(4, 3, 1, 's');
      case '4f_z(x2-y2)': return R(4, 3, 2, 'c');
      case '4f_xyz': return R(4, 3, 2, 's');
      case '4f_x(x2-3y2)': return R(4, 3, 3, 'c');
      case '4f_y(3x2-y2)': return R(4, 3, 3, 's');
      // n=5 轨道
      // 5s
      case '5s': return R(5, 0, 0, 'c');
      // 5p 轨道
      case '5pz': return R(5, 1, 0, 'c');
      case '5px': return R(5, 1, 1, 'c');
      case '5py': return R(5, 1, 1, 's');
      // 5d 轨道
      case '5d_z2': return R(5, 2, 0, 'c');
      case '5d_xz': return R(5, 2, 1, 'c');
      case '5d_yz': return R(5, 2, 1, 's');
      case '5d_xy': return R(5, 2, 2, 's');
      case '5d_x2-y2': return R(5, 2, 2, 'c');
      // 5f 轨道
      case '5f_z3': return R(5, 3, 0, 'c');
      case '5f_xz2': return R(5, 3, 1, 'c');
      case '5f_yz2': return R(5, 3, 1, 's');
      case '5f_z(x2-y2)': return R(5, 3, 2, 'c');
      case '5f_xyz': return R(5, 3, 2, 's');
      case '5f_x(x2-3y2)': return R(5, 3, 3, 'c');
      case '5f_y(3x2-y2)': return R(5, 3, 3, 's');
      // 5g 轨道 (l=4, m=-4,-3,-2,-1,0,1,2,3,4 共9个实轨道)
      case '5g_z4': return R(5, 4, 0, 'c');           // m=0
      case '5g_z3x': return R(5, 4, 1, 'c');          // m=1, cos
      case '5g_z3y': return R(5, 4, 1, 's');          // m=1, sin
      case '5g_z2xy': return R(5, 4, 2, 's');         // m=2, sin
      case '5g_z2(x2-y2)': return R(5, 4, 2, 'c');    // m=2, cos
      case '5g_zx(x2-3y2)': return R(5, 4, 3, 'c');   // m=3, cos
      case '5g_zy(3x2-y2)': return R(5, 4, 3, 's');   // m=3, sin
      case '5g_xy(x2-y2)': return R(5, 4, 4, 's');    // m=4, sin
      case '5g_x4-6x2y2+y4': return R(5, 4, 4, 'c');  // m=4, cos
      default: return R(1, 0, 0, 'c');
    }
  }

  // ==================== 杂化轨道高效采样 ====================
  // 
  // 【改进策略】
  // 使用混合提议分布：从各参与轨道的径向分布中采样
  // 然后使用重要性权重进行接受-拒绝
  // 
  // 这比纯拒绝采样效率高得多，因为提议分布自然覆盖了所有参与轨道的径向范围
  // ====================

  /**
   * 杂化轨道的精确采样（基于逆CDF混合）
   * 
   * 【物理原理】
   * 杂化轨道密度 |Σ c_i ψ_i|² 不是简单的各轨道密度之和
   * 但我们可以用各轨道的混合分布作为高效的提议分布
   * 
   * 【采样策略】
   * 1. 径向：从参与轨道中随机选择一个，使用其精确逆CDF采样
   * 2. 角向：均匀球面采样
   * 3. 接受-拒绝：按杂化密度权重决定接受
   * 
   * @param {string} atomType - 原子类型
   * @returns {Object|null} - 采样结果 { x, y, z, r, theta, phi, psi, accepted }
   */
  function hybridImportanceSample(paramsList, samplingBoundary, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return null;

    const numOrbitals = paramsList.length;
    const defaultCoeff = 1.0 / Math.sqrt(numOrbitals);

    // ==================== 第一步：从混合提议分布采样径向距离 ====================
    // 随机选择一个参与轨道，使用其精确逆CDF采样
    const orbitalIndex = Math.floor(Math.random() * numOrbitals);
    const selectedOrbital = paramsList[orbitalIndex];
    const r = sampleRadialExact(selectedOrbital.n, selectedOrbital.l);

    // 边界检查
    if (r > samplingBoundary * 2 || r < 1e-10) {
      return { x: 0, y: 0, z: 0, r, theta: 0, phi: 0, psi: 0, accepted: false };
    }

    // ==================== 第二步：均匀球面采样角度 ====================
    const { theta, phi } = sampleUniformSphere();

    // ==================== 第三步：计算杂化波函数和密度 ====================
    let psi = 0;
    for (const p of paramsList) {
      const coeff = p.coefficient !== undefined ? p.coefficient : defaultCoeff;
      const R = radialR(p.n, p.l, r, 1, 1, atomType);
      const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
      psi += coeff * R * Y;
    }
    const hybridDensity = psi * psi;

    // ==================== 第四步：计算提议分布密度 ====================
    // 提议分布 = (1/N) × Σ [P_i(r) × 1/(4π)]
    // 其中 P_i(r) = r² |R_i(r)|²
    let proposalDensity = 0;
    for (const p of paramsList) {
      const R = radialR(p.n, p.l, r, 1, 1, atomType);
      proposalDensity += r * r * R * R / (4 * PI);
    }
    proposalDensity /= numOrbitals;

    // ==================== 第五步：重要性权重接受-拒绝 ====================
    if (proposalDensity < 1e-300) {
      return { x: 0, y: 0, z: 0, r, theta, phi, psi, accepted: false };
    }

    // 权重 = 目标密度 / 提议密度
    const weight = hybridDensity / proposalDensity;

    // 估计最大权重（动态调整）
    // 由于杂化轨道的干涉效应，实际最大权重通常小于 4π
    // 但在节点附近可能会有极端值
    const maxWeight = 4 * PI * 1.5; // 保守估计

    // 接受-拒绝
    if (Math.random() * maxWeight > weight) {
      return { x: 0, y: 0, z: 0, r, theta, phi, psi, accepted: false };
    }

    // ==================== 采样成功 ====================
    const sinTheta = Math.sin(theta);
    const x = r * sinTheta * Math.cos(phi);
    const y = r * sinTheta * Math.sin(phi);
    const z = r * Math.cos(theta);

    // 添加亚像素抖动
    const dither = 0.01;
    const dx = (Math.random() - 0.5) * dither;
    const dy = (Math.random() - 0.5) * dither;
    const dz = (Math.random() - 0.5) * dither;

    return {
      x: x + dx, y: y + dy, z: z + dz,
      r, theta, phi,
      psi,  // 保留波函数值，用于相位着色
      accepted: true
    };
  }

  /**
   * 构建杂化轨道的径向CDF（用于最高效的采样）
   * 
   * 【物理说明】
   * 杂化轨道的径向概率密度不是简单的分离变量形式
   * P_hybrid(r) = r² × ∫∫ |Σ c_i R_i(r) Y_i(θ,φ)|² sin(θ) dθ dφ
   * 
   * 但由于球谐函数的正交性，交叉项在角向积分后消失
   * 所以 P_hybrid(r) = r² × Σ |c_i|² × |R_i(r)|²
   * 
   * @param {string} atomType - 原子类型
   * @returns {Object} { r, cdf, rMax, ... }
   */
  function buildHybridRadialCDF(paramsList, numPoints = 2000, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return null;

    // 生成缓存key
    const key = paramsList.map(p => `${p.n}_${p.l}`).join('|') + '_' + atomType;
    if (_cdfCache['hybrid_' + key]) {
      return _cdfCache['hybrid_' + key];
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

      // 计算 P(r) = r² × Σ |c_i|² × |R_i(r)|²
      let P_prev = 0, P_curr = 0;
      const r_prev = r[i - 1];
      const r_curr = r[i];

      for (const p of paramsList) {
        const coeff = p.coefficient !== undefined ? p.coefficient : defaultCoeff;
        const c2 = coeff * coeff;

        const R_prev = radialR(p.n, p.l, r_prev, 1, 1, atomType);
        const R_curr = radialR(p.n, p.l, r_curr, 1, 1, atomType);

        P_prev += c2 * r_prev * r_prev * R_prev * R_prev;
        P_curr += c2 * r_curr * r_curr * R_curr * R_curr;
      }

      // 梯形法则
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
    _cdfCache['hybrid_' + key] = result;
    return result;
  }

  /**
   * 杂化轨道的高精度采样（基于杂化CDF + 角向接受-拒绝）
   * 
   * 【最高效方法】
   * 1. 使用杂化轨道自己的径向CDF精确采样r
   * 2. 均匀球面采样角度
   * 3. 按角向权重接受-拒绝
   * 
   * 这个方法在理论上最接近"精确采样"，效率最高
   */
  function hybridPreciseSample(paramsList, samplingBoundary, atomType = 'H') {
    if (!paramsList || paramsList.length === 0) return null;

    const numOrbitals = paramsList.length;
    const defaultCoeff = 1.0 / Math.sqrt(numOrbitals);

    // 构建/获取杂化CDF
    const cdfData = buildHybridRadialCDF(paramsList, 2000, atomType);
    if (!cdfData) return null;

    // ==================== 第一步：从杂化径向CDF精确采样 ====================
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

    // ==================== 第二步：均匀球面采样角度 ====================
    const { theta, phi } = sampleUniformSphere();

    // ==================== 第三步：计算角向权重 ====================
    // 计算杂化波函数
    let psi = 0;
    let radialSum2 = 0; // Σ |c_i R_i(r)|²

    for (const p of paramsList) {
      const coeff = p.coefficient !== undefined ? p.coefficient : defaultCoeff;
      const R = radialR(p.n, p.l, r, 1, A0, atomType);
      const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
      psi += coeff * R * Y;
      radialSum2 += coeff * coeff * R * R;
    }

    // 角向权重 = |Σ c_i R_i Y_i|² / (Σ |c_i R_i|² × 1/(4π))
    // 这是因为径向采样已经按 Σ|c_i R_i|² 加权了
    // 均匀球面采样按 1/(4π) 分布
    const angularFactor = psi * psi;
    const expectedFactor = radialSum2 / (4 * PI);

    if (expectedFactor < 1e-300) {
      return { x: 0, y: 0, z: 0, r, theta, phi, psi, accepted: false };
    }

    const angularWeight = angularFactor / expectedFactor;

    // 【关键修复】角向权重的上界应与轨道数量成正比
    // 当多个轨道的球谐函数相干叠加时，|Σ c_i Y_i|² 可能远大于单轨道情况
    // 旧值 8π 导致高密度区域样本被错误拒绝，造成采样偏差
    const maxAngularWeight = 4 * PI * (numOrbitals + 2);

    // 角向接受-拒绝
    if (Math.random() * maxAngularWeight > angularWeight) {
      return { x: 0, y: 0, z: 0, r, theta, phi, psi, accepted: false };
    }

    // ==================== 采样成功 ====================
    const sinTheta = Math.sin(theta);
    const x = r * sinTheta * Math.cos(phi);
    const y = r * sinTheta * Math.sin(phi);
    const z = r * Math.cos(theta);

    // 添加亚像素抖动
    const dither = 0.01;
    const dx = (Math.random() - 0.5) * dither;
    const dy = (Math.random() - 0.5) * dither;
    const dz = (Math.random() - 0.5) * dither;

    return {
      x: x + dx, y: y + dy, z: z + dz,
      r, theta, phi,
      psi,
      accepted: true
    };
  }

  /**
   * 计算杂化轨道的主轴方向
   * 
   * 【物理原理】
   * 杂化轨道的主轴是角向分布最大值的方向
   * 例如 sp 杂化轨道的主轴在 p 轨道的方向
   * 
   * 【算法】
   * 1. 在单位球面上搜索角向强度 |Ψ(θ,φ)|² 的最大值方向
   * 2. 返回该方向的单位向量
   * 
   * @param {Array} paramsList - 轨道参数列表
   * @returns {Object} - {x, y, z} 主轴方向单位向量，以及 {theta, phi} 球坐标
   */
  function findHybridPrincipalAxis(paramsList) {
    if (!paramsList || paramsList.length === 0) {
      return { x: 0, y: 0, z: 1, theta: 0, phi: 0 }; // 默认 z 轴
    }

    // 如果只有一个 s 轨道（球对称），返回 z 轴
    if (paramsList.length === 1 && paramsList[0].l === 0) {
      return { x: 0, y: 0, z: 1, theta: 0, phi: 0 };
    }

    const defaultCoeff = 1.0 / Math.sqrt(paramsList.length);

    // 计算给定方向上的角向强度（只计算角向部分）
    function angularIntensity(theta, phi) {
      let psiSum = 0;
      for (const p of paramsList) {
        const coeff = p.coefficient !== undefined ? p.coefficient : defaultCoeff;
        const Y = realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
        psiSum += coeff * Y;
      }
      return psiSum * psiSum;
    }

    // 第一阶段：粗搜索 - 在整个球面上搜索
    let maxIntensity = 0;
    let bestTheta = 0;
    let bestPhi = 0;

    const thetaSteps = 36;
    const phiSteps = 72;

    for (let i = 0; i <= thetaSteps; i++) {
      const theta = (i / thetaSteps) * PI;
      for (let j = 0; j < phiSteps; j++) {
        const phi = (j / phiSteps) * TWO_PI;
        const intensity = angularIntensity(theta, phi);
        if (intensity > maxIntensity) {
          maxIntensity = intensity;
          bestTheta = theta;
          bestPhi = phi;
        }
      }
    }

    // 第二阶段：细化搜索 - 在最佳点附近精细搜索
    const dTheta = PI / thetaSteps;
    const dPhi = TWO_PI / phiSteps;

    for (let iter = 0; iter < 3; iter++) {
      const searchRange = dTheta / Math.pow(2, iter);
      let localMax = maxIntensity;
      let localBestTheta = bestTheta;
      let localBestPhi = bestPhi;

      for (let i = -5; i <= 5; i++) {
        const theta = Math.max(0, Math.min(PI, bestTheta + i * searchRange / 5));
        for (let j = -5; j <= 5; j++) {
          const phi = (bestPhi + j * searchRange / 5 + TWO_PI) % TWO_PI;
          const intensity = angularIntensity(theta, phi);
          if (intensity > localMax) {
            localMax = intensity;
            localBestTheta = theta;
            localBestPhi = phi;
          }
        }
      }
      maxIntensity = localMax;
      bestTheta = localBestTheta;
      bestPhi = localBestPhi;
    }

    // 将球坐标转换为笛卡尔坐标
    const x = Math.sin(bestTheta) * Math.cos(bestPhi);
    const y = Math.sin(bestTheta) * Math.sin(bestPhi);
    const z = Math.cos(bestTheta);

    return { x, y, z, theta: bestTheta, phi: bestPhi };
  }

  /**
   * 计算从 z 轴旋转到目标轴的四元数
   * 
   * @param {Object} targetAxis - {x, y, z} 目标轴方向单位向量
   * @returns {Object} - {x, y, z, w} 四元数
   */
  function getRotationQuaternionToAxis(targetAxis) {
    // 源轴：z 轴 (0, 0, 1)
    const sx = 0, sy = 0, sz = 1;
    // 目标轴
    const tx = targetAxis.x, ty = targetAxis.y, tz = targetAxis.z;

    // 计算旋转轴（叉乘）
    const cx = sy * tz - sz * ty;
    const cy = sz * tx - sx * tz;
    const cz = sx * ty - sy * tx;

    // 计算旋转角（点乘）
    const dot = sx * tx + sy * ty + sz * tz;

    // 如果几乎平行
    if (dot > 0.9999) {
      return { x: 0, y: 0, z: 0, w: 1 }; // 单位四元数
    }

    // 如果几乎反平行
    if (dot < -0.9999) {
      // 绕 x 轴旋转 180°
      return { x: 1, y: 0, z: 0, w: 0 };
    }

    // 标准化旋转轴
    const crossLen = Math.sqrt(cx * cx + cy * cy + cz * cz);
    const nx = cx / crossLen;
    const ny = cy / crossLen;
    const nz = cz / crossLen;

    // 计算半角
    const angle = Math.acos(dot);
    const halfAngle = angle / 2;
    const sinHalfAngle = Math.sin(halfAngle);
    const cosHalfAngle = Math.cos(halfAngle);

    return {
      x: nx * sinHalfAngle,
      y: ny * sinHalfAngle,
      z: nz * sinHalfAngle,
      w: cosHalfAngle
    };
  }

  /**
   * 获取轨道的最高次旋转对称轴方向
   * 
   * 基于理论计算，不依赖采样结果
   * 
   * 【对称性分析】
   * - s 轨道 (l=0)：球对称 O(3)，返回 z 轴
   * - p 轨道 (l=1)：
   *   - pz (m=0)：C∞v，σ∞ 轴在 z 方向
   *   - px (m=1, cos)：C∞v，σ∞ 轴在 x 方向
   *   - py (m=1, sin)：C∞v，σ∞ 轴在 y 方向
   * - d 轨道 (l=2)：
   *   - dz² (m=0)：D∞h，主轴在 z 方向
   *   - dxz (m=1, cos)：C2v，C2 轴在 y 方向（垂直于 xz 平面）
   *   - dyz (m=1, sin)：C2v，C2 轴在 x 方向（垂直于 yz 平面）
   *   - dx²-y² (m=2, cos)：D4h，C4 轴在 z 方向
   *   - dxy (m=2, sin)：D4h，C4 轴在 z 方向
   * - f 轨道 (l=3)：
   *   - fz³ (m=0)：主轴在 z 方向
   *   - 其他：根据对称性确定
   * 
   * @param {Object} angKey - 角向键 {l, m, t}
   * @returns {Object} - {x, y, z} 单位向量，最高次对称轴方向
   */
  function getOrbitalSymmetryAxis(angKey) {
    if (!angKey) return { x: 0, y: 0, z: 1 };

    const l = angKey.l;
    const m = angKey.m;
    const t = angKey.t; // 'c' for cos, 's' for sin

    // s 轨道：球对称，默认 z 轴
    if (l === 0) {
      return { x: 0, y: 0, z: 1 };
    }

    // p 轨道
    if (l === 1) {
      if (m === 0) {
        // pz：σ∞ 在 z 轴
        return { x: 0, y: 0, z: 1 };
      } else if (m === 1) {
        if (t === 'c') {
          // px：σ∞ 在 x 轴
          return { x: 1, y: 0, z: 0 };
        } else {
          // py：σ∞ 在 y 轴
          return { x: 0, y: 1, z: 0 };
        }
      }
    }

    // d 轨道
    if (l === 2) {
      if (m === 0) {
        // dz²：D∞h，主轴在 z
        return { x: 0, y: 0, z: 1 };
      } else if (m === 1) {
        if (t === 'c') {
          // dxz：C2v，C2 轴垂直于 xz 平面，即 y 方向
          return { x: 0, y: 1, z: 0 };
        } else {
          // dyz：C2v，C2 轴垂直于 yz 平面，即 x 方向
          return { x: 1, y: 0, z: 0 };
        }
      } else if (m === 2) {
        // dx²-y² 或 dxy：D4h，C4 轴在 z
        return { x: 0, y: 0, z: 1 };
      }
    }

    // f 轨道
    if (l === 3) {
      if (m === 0) {
        // fz³：主轴在 z
        return { x: 0, y: 0, z: 1 };
      } else if (m === 1) {
        // fxz², fyz²
        if (t === 'c') {
          return { x: 0, y: 1, z: 0 }; // 垂直于 xz 平面
        } else {
          return { x: 1, y: 0, z: 0 }; // 垂直于 yz 平面
        }
      } else if (m === 2) {
        // fz(x²-y²), fxyz
        return { x: 0, y: 0, z: 1 };
      } else if (m === 3) {
        // fx(x²-3y²), fy(3x²-y²)
        return { x: 0, y: 0, z: 1 };
      }
    }

    // 默认返回 z 轴
    return { x: 0, y: 0, z: 1 };
  }

  /**
   * 寻找杂化轨道的主轴方向（最大概率密度方向）
   * 用于对齐网格极点，避免极点处的纹理扭曲
   * 
   * 【改进】优先使用解析解（群论/几何构造），回退到数值搜索
   */
  function findHybridPrincipalAxis(paramsList, hybridIndex) {
    console.log('[DEBUG findHybridPrincipalAxis] hybridIndex:', hybridIndex, 'paramsList:', paramsList);
    if (!paramsList || paramsList.length === 0) return { x: 0, y: 0, z: 1 };

    // 移除所有硬编码的解析解检查
    // 直接使用下方的通用数值搜索，以保证对任意组合（如 2s + 2px）的正确性
    // 用户反馈表明解析解假设的标准方向（如 z 轴）往往与用户的实际选择（如 x 轴）不符

    // 【关键修复】对轨道进行排序，确保与 getHybridCoefficients 假设的顺序一致
    // getHybridCoefficients 期望的顺序是：s, p (px, py, pz), d (dz2, dx2-y2, ...)
    // 如果不排序，2s+2px 可能被错误地按 [2px, 2s] 传入，导致系数错配
    const sortedParams = sortOrbitalsForHybridization(paramsList);

    const numOrbitals = sortedParams.length;
    const coeffMatrix = getHybridCoefficients(numOrbitals, sortedParams);
    const coeffs = coeffMatrix[hybridIndex % numOrbitals];

    // 【关键修复】使用全3D密度计算，而非仅角向强度
    // 之前的 angularIntensity 忽略了径向波函数 R(r) 的符号
    // 例如 2s 轨道在远端是负的，而 2px 在 +x 是负的（由于 Condon-Shortley 相位）
    // 简单的角向叠加会导致错误的干涉预测（如预测在 -x 增强，实际在 +x 增强）
    // 使用 r = n^2 (玻尔半径单位) 作为特征半径进行采样
    const n = sortedParams[0].n;
    const r_characteristic = n * n;

    // 计算给定方向上的3D密度
    function getDensityAt(theta, phi) {
      return singleHybridDensity3D(sortedParams, hybridIndex, r_characteristic, theta, phi);
    }

    // 第一阶段：粗搜索 - 在整个球面上搜索
    let maxIntensity = -1;
    let bestTheta = 0;
    let bestPhi = 0;

    const thetaSteps = 36;
    const phiSteps = 72;

    for (let i = 0; i <= thetaSteps; i++) {
      const theta = (i / thetaSteps) * Math.PI;
      for (let j = 0; j < phiSteps; j++) {
        const phi = (j / phiSteps) * 2 * Math.PI;
        const intensity = getDensityAt(theta, phi);
        if (intensity > maxIntensity) {
          maxIntensity = intensity;
          bestTheta = theta;
          bestPhi = phi;
        }
      }
    }

    // 【修复】继续执行第二阶段细化搜索，而不是提前返回
    // 第二阶段：细化搜索 - 在最佳点附近精细搜索
    const dTheta = Math.PI / thetaSteps;

    for (let iter = 0; iter < 3; iter++) {
      const searchRange = dTheta / Math.pow(2, iter);
      let localMax = maxIntensity;
      let localBestTheta = bestTheta;
      let localBestPhi = bestPhi;

      for (let i = -5; i <= 5; i++) {
        const theta = Math.max(0, Math.min(Math.PI, bestTheta + i * searchRange / 5));
        for (let j = -5; j <= 5; j++) {
          const phi = (bestPhi + j * searchRange / 5 + 2 * Math.PI) % (2 * Math.PI);
          const intensity = getDensityAt(theta, phi);
          if (intensity > localMax) {
            localMax = intensity;
            localBestTheta = theta;
            localBestPhi = phi;
          }
        }
      }
      maxIntensity = localMax;
      bestTheta = localBestTheta;
      bestPhi = localBestPhi;
    }

    // 转换为直角坐标
    const sinTheta = Math.sin(bestTheta);
    let axis = {
      x: sinTheta * Math.cos(bestPhi),
      y: sinTheta * Math.sin(bestPhi),
      z: Math.cos(bestTheta)
    };

    // 【无需反转】
    // 因为我们直接搜索的是概率密度 |Ψ|² 的最大值
    // 所以找到的方向一定是指向电子云密度最大的"瓣"
    // 不需要再根据波函数符号进行反转（密度总是正的）

    console.log('[DEBUG findHybridPrincipalAxis] Final axis:', axis, 'maxDensity:', maxIntensity);
    return axis;
  }

  /**
   * 计算从 z 轴旋转到目标轴的四元数
   * 
   * @param {Object} targetAxis - {x, y, z} 目标轴方向单位向量
   * @returns {Object} - {x, y, z, w} 四元数
   */
  /**
   * Group Theory / Composition Analysis for "All" Mode
   * Determines the highest symmetry axis based on the set of constituent orbitals.
   */
  function getSetSymmetryAxis(paramsList) {
    if (!paramsList || paramsList.length === 0) return { x: 0, y: 0, z: 1 };

    let hasPx = false;
    let hasPy = false;
    let hasPz = false;

    for (const p of paramsList) {
      if (p.l === 1) { // p orbital
        if (p.angKey.m === 0) hasPz = true;
        else if (p.angKey.m === 1 && p.angKey.t === 'c') hasPx = true;
        else if (p.angKey.m === 1 && p.angKey.t === 's') hasPy = true;
      }
    }

    const pCount = (hasPx ? 1 : 0) + (hasPy ? 1 : 0) + (hasPz ? 1 : 0);

    if (pCount === 1) {
      if (hasPx) return { x: 1, y: 0, z: 0 };
      if (hasPy) return { x: 0, y: 1, z: 0 };
      if (hasPz) return { x: 0, y: 0, z: 1 };
    }

    if (pCount === 2) {
      if (hasPx && hasPy) return { x: 0, y: 0, z: 1 }; // xy plane -> z
      if (hasPx && hasPz) return { x: 0, y: 1, z: 0 }; // xz plane -> y
      if (hasPy && hasPz) return { x: 1, y: 0, z: 0 }; // yz plane -> x
    }

    return { x: 0, y: 0, z: 1 };
  }

  /**
   * Analytical Lobe Axis for "Single" Hybrid Mode
   * Calculates the direction vector of the specific hybrid lobe.
   */
  function getHybridLobeAxis(paramsList, hybridIndex) {
    if (!paramsList || paramsList.length === 0) return { x: 0, y: 0, z: 1 };

    const sortedParams = sortOrbitalsForHybridization(paramsList);
    const numOrbitals = sortedParams.length;
    const coeffMatrix = getHybridCoefficients(numOrbitals, sortedParams);
    const coeffs = coeffMatrix[hybridIndex % numOrbitals];

    let vecX = 0, vecY = 0, vecZ = 0;

    for (let i = 0; i < numOrbitals; i++) {
      const p = sortedParams[i];
      const c = coeffs[i];

      if (p.l === 1) {
        // p 轨道贡献：realYlm_value 已经处理了化学惯例
        // 正系数 -> 指向正轴方向
        if (p.angKey.m === 0) vecZ += c;  // pz
        else if (p.angKey.m === 1 && p.angKey.t === 'c') vecX += c;  // px
        else if (p.angKey.m === 1 && p.angKey.t === 's') vecY += c;  // py
      }
    }

    const len = Math.sqrt(vecX * vecX + vecY * vecY + vecZ * vecZ);
    if (len < 1e-6) return { x: 0, y: 0, z: 1 };

    return { x: vecX / len, y: vecY / len, z: vecZ / len };
  }

  function getRotationToAxis(targetAxis) {
    // 源轴：z 轴 (0, 0, 1)
    const tx = targetAxis.x, ty = targetAxis.y, tz = targetAxis.z;

    // 如果已经对齐 z 轴，返回单位四元数
    if (Math.abs(tz - 1) < 1e-6) {
      return { x: 0, y: 0, z: 0, w: 1 };
    }

    // 如果反向（指向 -z），绕 x 轴旋转 180°
    if (Math.abs(tz + 1) < 1e-6) {
      return { x: 1, y: 0, z: 0, w: 0 };
    }

    // 计算旋转轴（z × target）
    const cx = -ty;  // 0*tz - 1*ty
    const cy = tx;   // 1*tx - 0*tz
    const cz = 0;    // 0*ty - 0*tx

    // 标准化旋转轴
    const crossLen = Math.sqrt(cx * cx + cy * cy);
    const nx = cx / crossLen;
    const ny = cy / crossLen;

    // 计算半角
    const dot = tz; // (0,0,1) · (tx,ty,tz)
    const angle = Math.acos(Math.max(-1, Math.min(1, dot)));
    const halfAngle = angle / 2;
    const sinHalfAngle = Math.sin(halfAngle);
    const cosHalfAngle = Math.cos(halfAngle);

    return {
      x: nx * sinHalfAngle,
      y: ny * sinHalfAngle,
      z: 0,
      w: cosHalfAngle
    };
  }

  /**
   * 估算轨道的 95% 等值面半径
   * 用于渲染前自动缩放相机和调整点大小
   * 
   * @param {string} atomType - 原子类型（'H', 'C', 'Kr' 等）
   * @param {string} orbitalKey - 轨道键值（'1s', '2px', '3d_xy' 等）
   * @returns {number} - 估算的 95% 等值面半径（玻尔半径单位）
   */
  function estimateOrbitalRadius95(atomType, orbitalKey) {
    // 从轨道键值解析 n 和 l
    const params = orbitalParamsFromKey(orbitalKey);
    if (!params) return 15; // 默认值

    const n = params.n;
    const l = params.l;

    // 氢原子：使用解析解公式 15 × n²
    if (!atomType || atomType === 'H') {
      return 15 * n * n;
    }

    // 非氢原子：基于 STO 数据估算
    const orbitalId = getOrbitalKey(n, l); // 例如 '2s', '3p'
    const SlaterBasis = window.SlaterBasis;

    if (!SlaterBasis || !SlaterBasis[atomType] || !SlaterBasis[atomType].orbitals) {
      // 无 STO 数据，回退到氢原子公式
      return 15 * n * n;
    }

    const basis = SlaterBasis[atomType].orbitals[orbitalId];
    if (!basis || basis.length === 0) {
      return 15 * n * n;
    }

    // 【优化】使用有效最小 Zeta（最慢衰减项）来估算半径
    // 平均 Zeta 会被内层高 Zeta 项拉高，导致低估外层轨道大小
    let minZeta = Infinity;
    let maxN_atMinZeta = 1;

    // 1. 找到对轨道外沿有实质贡献的最小 Zeta
    // 忽略系数极小 (< 0.05 maxCoeff) 的长尾项，防止数值噪声干扰
    let maxCoeff = 0;
    for (const term of basis) {
      if (Math.abs(term.coeff) > maxCoeff) maxCoeff = Math.abs(term.coeff);
    }
    const threshold = maxCoeff * 0.05;

    let found = false;
    for (const term of basis) {
      // 只考虑有实质贡献的项
      if (Math.abs(term.coeff) > threshold) {
        if (term.zeta < minZeta) {
          minZeta = term.zeta;
          maxN_atMinZeta = term.nStar;
          found = true;
        }
      }
    }

    // 兜底：如果没找到（不应该发生），使用第一项
    if (!found) {
      minZeta = basis[0].zeta;
      maxN_atMinZeta = basis[0].nStar;
    }

    // STO 峰值半径 ≈ n* / zeta
    // 95% 等值面通常延伸到峰值的 3-4 倍远处（取决于 exp(-zeta*r) 的衰减）
    // 对于漫散轨道（小 zeta），需要更大的倍数
    const rPeak = maxN_atMinZeta / minZeta;
    const factor = 4.0;
    const r95 = rPeak * factor;

    // 返回估算值，但设置合理的最小值
    return Math.max(2, r95);
  }

  window.Hydrogen = {
    A0,
    radialR,
    radialWavefunction: radialR, // 别名：用于兼容性
    radialPDF,
    density3D_real,
    realYlm_abs2,
    realYlm_value,
    radialGrid,
    histogramRadialFromSamples,
    histogramThetaFromSamples,
    histogramPhiFromSamples,
    recommendRmax,
    orbitalParamsFromKey,
    orbitalKey,
    angularPDF_Theta,
    angularPDF_Phi,
    // 精确逆CDF采样
    buildRadialCDF,
    sampleRadialExact,
    // 有效核电荷（Slater规则）
    getEffectiveZ,
    // 重要性采样相关函数（保留兼容性）
    findRadialPeaks,
    getCachedPeaks,
    getMultiPeakMixtureParams,
    radialProposalPDF,
    sampleRadialProposal,
    sampleUniformSphere,
    getMaxRadialWeight,
    importanceSample,
    batchImportanceSample,
    // 杂化轨道相关函数
    hybridDensity3D,
    hybridRadialPDF, // 杂化轨道径向PDF（理论曲线用）
    allHybridOrbitalsDensity3D, // 所有杂化轨道总密度
    hybridRecommendRmax,
    hybridEstimateMaxDensity,
    hybridImportanceSample,
    buildHybridRadialCDF,
    hybridPreciseSample,
    // 单个杂化轨道
    getHybridCoefficients,
    singleHybridDensity3D,
    singleHybridWavefunction,
    sortOrbitalsForHybridization,
    findHybridPrincipalAxis,
    getSetSymmetryAxis,
    getHybridLobeAxis,
    // 轨道对称轴
    getOrbitalSymmetryAxis,
    getRotationToAxis,
    slaterRadialR, // Export STO function
    // 轨道半径估算（用于自动缩放）
    estimateOrbitalRadius95,

    // === 势能积分曲线相关 ===

    /**
     * 计算理论势能积分曲线
     * E(r) = -Z * ∫[0->r] x * R(x)^2 dx
     * 
     * @param {number} n - 主量子数
     * @param {number} l - 角量子数
     * @param {number} Z - 核电荷数
     * @param {string} atomType - 原子类型
     * @param {number} rMax - 最大半径
     * @param {number} steps - 积分步数
     * @returns {Object} { r: [], E: [] }
     */
    calculateCumulativePotential: function (n, l, Z, atomType, rMax, steps = 500) {
      const dr = rMax / steps;
      const rValues = new Float32Array(steps);
      const eValues = new Float32Array(steps); // 累积势能

      let integral = 0;
      // 使用梯形法则积分
      // dE/dr = (-Z/r) × P(r) = (-Z/r) × r²R² = -Z × r × R²

      let prevValue = 0;

      for (let i = 0; i < steps; i++) {
        const r = (i + 1) * dr;
        const R = radialR(n, l, r, 1, A0, atomType);
        const currentValue = -Z * r * R * R;

        // 梯形面积
        const area = 0.5 * (prevValue + currentValue) * dr;
        integral += area;

        rValues[i] = r;
        eValues[i] = integral;
        prevValue = currentValue;
      }

      return { r: rValues, E: eValues };
    },

    /**
     * 将概率密度直方图转换为势能积分曲线
     * 用于处理采样数据
     * 
     * @param {Array} counts - 直方图计数
     * @param {Array} edges - 直方图边界
     * @param {number} totalSamples - 总样本数
     * @param {number} Z - 核电荷数
     * @returns {Array} 累积能量值数组 (与counts长度相同)
     */
    transformHistogramToPotential: function (counts, edges, totalSamples, Z) {
      if (totalSamples === 0) return new Float32Array(counts.length);

      const eValues = new Float32Array(counts.length);
      let cumulativeE = 0;

      for (let i = 0; i < counts.length; i++) {
        const count = counts[i];
        if (count > 0) {
          // Bin center r
          const r = 0.5 * (edges[i] + edges[i + 1]);
          // 概率 P_i = count / totalSamples
          const prob = count / totalSamples;
          // 该 bin 的平均势能贡献贡献: P_i * (-Z / r)
          // 注意：这里我们用 bin 中心的势能近似整个 bin 的平均势能
          const contribution = prob * (-Z / r);
          cumulativeE += contribution;
        }
        eValues[i] = cumulativeE;
      }
      return eValues;
    }
  };
})();
