# 电子云可视化 - 算法技术文档

本文档详细说明了本项目中的核心物理算法和数学实现细节，涵盖多电子原子计算、杂化轨道处理及采样系统。

---

## 目录

1. [多电子原子物理引擎 (STO)](#1-多电子原子物理引擎-sto)
2. [相位校正算法](#2-相位校正算法)
3. [杂化轨道叠加计算](#3-杂化轨道叠加计算)
4. [半径估算方法](#4-半径估算方法)
5. [采样系统](#5-采样系统)

---

## 1. 多电子原子物理引擎 (STO)

**核心文件**: `physics.js`, `slater_basis.js`

对于非氢原子（He 到 Xe），项目采用 **Slater-Type Orbitals (STO)** 展开法来逼近 Hartree-Fock 极限。

### 1.1 数学形式
径向波函数 $R_{nl}(r)$ 表示为 Slater 基函数的线性组合：

$$R_{nl}(r) = \sum_{i=1}^M c_i \cdot \phi_i(r; n^*_i, \zeta_i)$$

其中基函数 $\phi_i$ 定义为：

$$\phi_i(r) = N_i \cdot r^{n^*_i - 1} \cdot e^{-\zeta_i r}$$

- $n^*_i$：有效主量子数
- $\zeta_i$：轨道指数
- $c_i$：展开系数
- $N_i$：归一化常数

### 1.2 数据来源
使用 **Clementi & Roetti (1974)** 的 Double-Zeta 及扩展基组数据。

---

## 2. 相位校正算法

**位置**: `physics.js` -> `radialR`

为解决 Clementi-Roetti 数据与物理学通用 **Condon-Shortley** 约定之间的符号不一致（如 Li-2s），引入了**渐近相位校正**。

### 2.1 算法逻辑
1. **计算目标相位**：根据氢原子解析解，轨道在无穷远处 ($r \to \infty$) 的符号应为：
   $$Sign_{target} = (-1)^{n-l-1}$$

2. **分析 STO 行为**：找到 STO 展开中 $\zeta$ 最小（衰减最慢）的主导项。
   
3. **校正**：
   若主导项符号与目标相位不符，则对计算结果取反：
   $$R_{final}(r) = -1 \cdot R_{STO}(r)$$

此方法确保了不同原子间波函数拓扑结构的一致性。

---

## 3. 杂化轨道叠加计算

**位置**: `physics.js` -> `hybridRadialPDF`, `sampling-worker.js` -> `buildHybridRadialCDF`

杂化轨道是原子轨道的线性组合：$\Psi = \sum c_i \psi_i$。在计算概率密度 $|\Psi|^2$ 时，必须包含**交叉项 (Cross-terms)**。

### 3.1 密度公式
概率密度函数 (PDF) 为：

$$P(\mathbf{r}) = r^2 \left( \sum_i |c_i R_i Y_i|^2 + \sum_{i \neq j} c_i c_j R_i R_j Y_i Y_j^* \right)$$

### 3.2 积分处理
在构建径向采样分布 (CDF) 时，通过对全立体角积分消除正交项。
- **正交组合** ($l$ 或 $m$ 不同)：交叉项积分为零。
- **非正交组合** (如 $1s + 2s$)：保留交叉项。

代码中实现了双重循环积分，以正确处理所有情况。

---

## 4. 半径估算方法

**位置**: `physics.js` -> `estimateOrbitalRadius95`

对于多指数 STO 轨道，采用**最小指数优先 (Minimum Zeta)** 策略来确定渲染边界。

1. 筛选基组中系数显著的项。
2. 找到最小的 $\zeta_{min}$。
3. 估算 95% 半径：
   $$r_{95} \approx \frac{n^*}{\zeta_{min}} \times 4.0$$

该方法能更准确地反映重原子较为弥散的价层电子云大小。

---

## 5. 采样系统

### 5.1 架构
- **多线程**：使用 Web Worker 进行后台采样。
- **会话管理**：通过 Session ID 区分不同采样任务，防止数据冲突。

### 5.2 滚动更新 (Rolling Update)
在比照模式下，采用滚动更新机制：
- 每帧随机替换约 1% 的旧点。
- 保持视觉连续性，避免画面闪烁。

---

## 6. 势能可视化计算

**核心文件**: `physics.js`, `data_panel.js`, `orbital.js`

为了展示电子-核库仑相互作用的径向分布，项目新增了势能积分和势能密度的可视化功能。

### 6.1 物理公式

电子-核库仑势能的径向分布：

$$\frac{dE}{dr} = -\frac{Z}{r} \cdot P(r)$$

其中 $P(r) = r^2 |R(r)|^2$ 是径向概率密度 (RDF)。

累积势能积分：

$$E(r) = \int_0^r \frac{dE}{dr'} dr' = -Z \int_0^r \frac{P(r')}{r'} dr'$$

等价形式：

$$E(r) = -Z \int_0^r r' \cdot |R(r')|^2 dr'$$

### 6.2 理论曲线计算

**函数**: `calculateCumulativePotential(n, l, Z, atomType, rMax, steps)`

使用梯形法则进行数值积分：

```javascript
for (let i = 0; i < steps; i++) {
    const r = (i + 1) * dr;
    const R = radialR(n, l, r, Z, A0, atomType);
    const currentValue = -Z * r * R * R;
    const area = 0.5 * (prevValue + currentValue) * dr;
    integral += area;
}
```

### 6.3 采样数据转换

**函数**: `transformHistogramToPotential(counts, edges, scaleFactor, Z)`

将归一化的概率密度直方图转换为累积势能：

1. 每个 bin 的概率：$P_i = \text{counts}[i] \times dr$
2. 每个 bin 的势能贡献：$\Delta E_i = P_i \times (-Z/r_i)$
3. 累积求和：$E_i = \sum_{j=0}^{i} \Delta E_j$

### 6.4 对数尺度表示

**dE/d(log r)**：乘以 $r$ 因子以适应对数横轴：

$$\frac{dE}{d(\log r)} = r \cdot \frac{dE}{dr} = -Z \cdot P(r)$$

这种表示方式使得内层电子（小 r）和外层电子（大 r）的势能贡献在图表上更加均衡可见。

---

*文档日期: 2025-12-18*
