# 核心算法与物理口径说明（重写版，保留全部信息）

更新时间：2025-12-29

本文档的定位是：用**心平气和、可复现、容易核对**的方式说明本项目在“浏览器端实时可视化”这一工程目标下，采用了哪些物理近似、哪些数学定义、以及这些定义在 UI 中对应哪些图表。

请注意：这里的“算法”并不是最重要的主题；更重要的是：**哪些量是被计算的、如何计算、哪些量容易被误读、以及数据/公式来自哪里**。为满足“内容不能删减”的要求，本文末尾保留了当前版本全文作为附录。

---

## 0. 核心约定与阅读方式

### 0.1 单位制

本项目所有计算默认使用**原子单位制 (Atomic Units, a.u.)**：
$$ \hbar = m_e = e = \frac{1}{4\pi\epsilon_0} = 1 $$
长度单位为玻尔半径 $a_0$，能量单位为 Hartree $E_h$。

### 0.2 一个重要“强调”（避免误读）

UI 中出现的 $V_{ee}(r)$、$Z_{eff}(r)$、$dE/dR$、$t_{loc}(r)$ 等曲线，均应理解为在如下实现框架下的**模型内自洽诊断量**：

- 中心场/构型平均（球平均的径向口径）
- 轨道径向来自 STO 数据表（Koga 系列）
- 电子-电子项使用（表格化的）有效势 $V_{ee}^{(i)}(r)$

它们用于解释“在这个实现框架里，哪些半径区间贡献更大”，但**不等同于严格可观测量**，也不等同于多体能量密度分解的唯一方案。

---

## 1. 物理模型边界与局限性（先讲清楚再看图）

本项目不是全功能量化计算软件（例如 Gaussian、ORCA、VASP 等），而是面向**实时可视化与交互诊断**的物理引擎。这里把“能做到什么/不能做到什么”放在开头，是为了让后面的公式与曲线更不容易被误解。

### 1.1 哈密顿量与近似层级

我们工作的出发点是非相对论、定核、定态薛定谔哈密顿量：
$$ H = -\sum_i \frac{1}{2}\nabla_i^2 - \sum_i \frac{Z}{r_i} + \sum_{i<j} \frac{1}{|\mathbf{r}_i - \mathbf{r}_j|} $$

在实现上，关键近似是：

- **非相对论近似**：不求解 Dirac-HF，不包含自旋-轨道耦合导致的精细结构分裂（如 $2p_{1/2}$ 与 $2p_{3/2}$）。
- **单电子轨道口径（HF 轨道）**：展示的是 Hartree-Fock 意义下的单电子轨道形状；因此不包含电子关联能（correlation）。
- **工程边界（缓存覆盖范围）**：与 $V_{ee}$/$Z_{eff}$/能量分解相关的查表数据仅覆盖 H–Kr（见 [vee_cache.js](vee_cache.js) 的文件头注释，以及生成脚本 [ _experiments/generate_production_cache.py ](_experiments/generate_production_cache.py)）。

### 1.2 冻结核近似（基态参数用于其它态）

[slater_basis.js](slater_basis.js) 中的基组参数针对**基态原子**优化。

- 当用户选择激发态（例如 C 的 $5d$）时，仍沿用基态屏蔽/基组参数。
- 这会引入明显误差：激发态电子感受到的有效核电荷 $Z_{eff}$ 可能与基态不同。此时可视化更接近“拓扑/定性正确”，半径与能量相关曲线不应被过度解读。

---

## 2. 波函数的数学定义（项目里真正用到的口径）

本项目的采样与可视化围绕单电子轨道波函数 $\psi(\mathbf{r})$ 展开。

### 2.1 分离变量与坐标

采用球坐标：
$$ r = \sqrt{x^2+y^2+z^2},\quad \theta = \arccos\left(\frac{z}{r}\right),\quad \phi = \operatorname{atan2}(y,x) $$

并写成径向与角向的乘积（或线性组合）：
$$ \psi(\mathbf{r}) = \sum_i R_i(r)\,Y_i(\theta,\phi) $$
其中 $i$ 可以表示单个轨道，也可以表示杂化轨道的线性组合通道。

### 2.2 复球谐与实球谐（化学惯例）

实现位于 [physics-core.js](physics-core.js) 的 `Ylm_complex`/`realYlm_value`。

复球谐的标准形式（含 Condon–Shortley 相位）为：
$$ Y_l^m(\theta,\phi) = (-1)^m\,\sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}\,P_l^m(\cos\theta)\,e^{im\phi} $$

实球谐通过 $Y_l^{m}$ 与 $Y_l^{-m}$ 的线性组合得到，并区分“cos 型/ sin 型”（代码中用 `type` 表示）：
$$ Y_{l,m}^{(c)} \propto \sqrt{2}\,\Re(Y_l^{|m|}),\qquad Y_{l,m}^{(s)} \propto \sqrt{2}\,\Im(Y_l^{|m|}) $$

强调：采用该组合是为了与常见量子化学/教材里的 $p_x,p_y$ 等实轨道约定保持一致，避免符号混乱。

### 2.3 氢原子径向函数（解析式 + 数值稳定性）

氢样原子径向波函数的常见写法是：
$$ R_{nl}(r) = \mathcal{N}_{nl}\,e^{-\rho/2}\,\rho^l\,L_{n-l-1}^{2l+1}(\rho),\quad \rho=\frac{2Zr}{n a_0} $$

本项目实现中强调“数值稳定性”而不是追求符号上的闭式展开：

- 广义拉盖尔多项式使用递推而非直接封闭形式，以避免 $n$ 较大时的溢出。
- 三项递推关系（项目文内保留原式）：
$$ (n+1)L_{n+1}^{(\alpha)}(x) = (2n+1+\alpha-x)L_n^{(\alpha)}(x) - (n+\alpha)L_{n-1}^{(\alpha)}(x) $$

### 2.4 多电子原子：Koga-STO 数据表径向函数

对于多电子原子，本项目不在浏览器端做 Roothaan-HF 自洽迭代，而使用 Koga 系列的 STO 展开数据（见 [slater_basis.js](slater_basis.js)）。

实现位于 [physics-core.js](physics-core.js) 的 `slaterRadialR`：

1) 每一项 STO 基函数采用
$$ \chi_j(r) = N_j\,r^{n_j^*-1}e^{-\zeta_j r} $$
2) 归一化常数在实现中采用
$$ N_j = \frac{(2\zeta_j)^{n_j^*+1/2}}{\sqrt{(2n_j^*)!}} $$
3) 径向函数按线性组合
$$ R(r)=\sum_{j=1}^{M} c_j\,\chi_j(r) $$
其中 $\{c_j,n_j^*,\zeta_j\}$ 由数据表给出。

---

## 3. $Z_{eff}(r)$ 与 $V_{ee}^{(i)}(r)$：定义、插值与边界条件

### 3.1 定义（从总势能反推等效电荷）

在“类氢势”的重写形式下定义：
$$ V(r) = -\frac{Z_{eff}(r)}{r} $$
并取
$$ V(r) = -\frac{Z}{r} + V_{ee}(r) $$
得到
$$ Z_{eff}(r) = Z - r\,V_{ee}(r) $$

该实现入口为 [physics.js](physics.js) 的 `calculateZeff`。

### 3.2 $V_{ee}^{(i)}(r)$ 的生成与意义边界

$V_{ee}^{(i)}(r)$ 是用于可视化的一种“中心场有效势”（直接库仑项 + 交换项的球平均/径向表达），以查表形式存储在 `VeeCache`（见 [vee_cache.js](vee_cache.js)）。

关于交换项角系数：生成脚本 [ _experiments/generate_production_cache.py ](_experiments/generate_production_cache.py) 中的 `EXCHANGE_COEFFS` 作用相当于原子结构理论中 Slater–Condon/Racah 角动量代数（可由 Gaunt 系数/Wigner $3j$ 等推导）的“角部分系数表”。本项目只是在工程上复用标准角动量结果来生成缓存数据，不宣称提出新的交换近似或新系数。

### 3.3 缓存插值

运行时对离散网格 $\{r_k\}$ 的 $V_{ee}(r_k)$ 采用二分查找 + 线性插值：
$$ V_{ee}(r) \approx (1-t)V_{ee}(r_k) + tV_{ee}(r_{k+1}),\quad t=\frac{r-r_k}{r_{k+1}-r_k} $$

### 3.4 远场边界条件（需要显式满足）

对中性原子 $N=Z$，当 $r\to\infty$ 且电子处在原子外部时，应有
$$ V_{ee}^{(i)}(r) \to \frac{N-1}{r},\qquad Z_{eff}(r)=Z-rV_{ee}(r)\to 1 $$

这一约束在生成端与运行端都被显式实现：

- 生成端在缓存尾部强制 $V_{ee}(r)=(N-1)/r$（避免数值尾部残余）。
- 运行端对 $r>r_{max}$ 做 $1/r$ 外推（避免常数夹断带来的非物理长程尾）。

---

## 4. 能量类图表：严格定义 + 误读提醒

本节的目标不是“争论物理对不对”，而是把项目里图表的**数学定义**写清楚，从而让读者能判断：

- 某条曲线是严格定义的量（在模型内）；
- 还是一个工程上方便的诊断量；
- 以及哪些解释是超出了该模型的。

### 4.1 径向概率密度

对一个球平均口径的轨道 $i$，径向概率密度定义为：
$$ P_i(r)=r^2|R_i(r)|^2,\qquad \int_0^\infty P_i(r)\,dr=1 $$

### 4.2 “累计贡献曲线”与“密度曲线”的关系

对任意局域径向算符 $U(r)$（例如 $-Z/r$ 或 $V_{ee}(r)$），期望值
$$ \langle U\rangle_i = \int_0^{\infty} P_i(r)\,U(r)\,dr $$

项目绘制的累计曲线在数学上是
$$ \langle U \rangle_i(<R)=\int_0^R P_i(r)\,U(r)\,dr $$
因此密度曲线（导数）严格为
$$ \frac{d}{dR}\langle U \rangle_i(<R)=P_i(R)\,U(R) $$

对应实现入口是 [physics.js](physics.js) 的 `calculateCumulativeOrbitalEnergy`。

### 4.3 局域动能 / 局域能量（代码曾称“重构”）

如果轨道满足一个局域势的一电子本征方程
$$ \left[-\frac12\nabla^2+V_{eff}(r)\right]\phi_i=\epsilon_i\phi_i $$
则有点态恒等式
$$ t_{loc,i}(r) \equiv \frac{-\frac12\nabla^2\phi_i}{\phi_i}=\epsilon_i-V_{eff}(r) $$

项目使用 $V_{eff}(r)=V_{nuc}(r)+V_{ee}^{(i)}(r)$ 构造 $t_{loc}$，并绘制 $t_{loc}(r)P(r)$ 等径向密度曲线。

强调：在 HF 中交换算符严格非局域；把交换项表成一个可乘的 $V_{ee}^{(i)}(r)$ 属于中心场/构型平均框架下的表达与近似。因此这些能量密度曲线属于“模型内部自洽的诊断量”，不应被当作严格可观测量。

---

## 5. 采样算法（点云生成）

### 5.1 径向：逆 CDF 采样

目标分布 $P(r)=r^2|R(r)|^2$。

1) 构造 CDF
$$ F(r)=\int_0^r P(t)\,dt $$
2) 生成 $\xi\sim\mathcal{U}(0,1)$，求 $F(r)=\xi$ 得到样本 $r$。

实现见 [physics.js](physics.js) 的 `buildRadialCDF`，并使用二分查找 + 线性插值实现逆变换。

### 5.2 角向：拒绝采样

在单位球面上均匀采样方向 $(\theta,\phi)$，以
$$ a = \frac{|Y(\theta,\phi)|^2}{M_{max}} $$
作为接受概率（$M_{max}$ 为上界估计），生成 $u\sim\mathcal{U}(0,1)$：若 $u<a$ 则接受，否则拒绝重采。

### 5.3 抖动（用于削弱离散 CDF 带来的环纹）

为削弱 CDF 离散化导致的同心圆伪影，对最终坐标加入小幅均匀抖动：
$$ x'=x+\delta_x,\ y'=y+\delta_y,\ z'=z+\delta_z,\quad \delta_*\sim\mathcal{U}(-\Delta,\Delta) $$
其中 $\Delta\approx 0.005\,a_0$（实现见 [sampling.js](sampling.js)）。

---

## 6. 95% 轮廓与等值面（Marching Cubes + Worker）

本项目存在两类“95%”相关量，容易混淆，这里把计算口径分开写：

### 6.1 $r_{95}$（半径分位数）

若对采样点的半径数组 $\{r_i\}$ 排序，则
$$ r_{95}=\mathrm{quantile}_{0.95}(\{r_i\}) $$
对应实现见 [visualization.js](visualization.js) 的 `calculate95PercentileRadius`。

### 6.2 等值面阈值（以 $|\psi|$ 分位选取）

网格等值面使用 Marching Cubes，阈值从采样点上的 $|\psi(\mathbf{r}_i)|$ 取分位数：
$$ \psi_{iso}=\mathrm{quantile}_{0.95}(\{|\psi(\mathbf{r}_i)|\}) $$
然后提取
$$ |\psi(\mathbf{r})| = \psi_{iso} $$

这是一种工程定义：它与“包含 95% 概率质量的等密度面”并不严格等价，但在交互式可视化里可提供稳定、可重复的阈值选择。

### 6.3 Worker 并行与连通域分离

等值面计算放入 Web Worker（见 [isosurface-worker.js](isosurface-worker.js)），主线程与 Worker 使用 Transferable 的 `ArrayBuffer` 传输顶点数组以降低拷贝开销。

Marching Cubes 输出三角形后进行连通域分离（BFS/并查集皆可），并按每个连通子网格中心点的 $\psi$ 符号进行上色：

- $\psi>0$ 与 $\psi<0$ 分别使用两种颜色
- 该策略对 $p/d/f$ 等多波瓣形状较鲁棒，无需硬编码几何规则

---

## 7. 杂化轨道自适应对齐（四元数优化）

为减少随机几何构型（Thomson 点等）与各向异性基组之间的“极点不对齐”，引入四元数旋转优化器。

目标函数写作
$$ J(q)=\sum_k \log\big(\sigma_k(A(q))\big) $$
其中 $\sigma_k$ 为投影矩阵的奇异值。实现与四元数/小型线代工具见 [physics-core.js](physics-core.js)。

实现细节（保持原信息不删减）：

- 随机重启爬山：20 次重启、每次 50 步迭代
- 计算量通常 <5ms

---

## 8. 公式总表（便于对照实现）

这一节把项目中“确实在计算里用到”的公式以索引形式汇总，便于你逐项对照代码。

1) 原子单位制：$\hbar=m_e=e=1/(4\pi\epsilon_0)=1$

2) 球坐标：$r=\sqrt{x^2+y^2+z^2}$，$\theta=\arccos(z/r)$，$\phi=\operatorname{atan2}(y,x)$

3) 复球谐：
$$ Y_l^m=(-1)^m\sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}P_l^m(\cos\theta)e^{im\phi} $$

4) 实球谐（cos/sin 线性组合）：
$$ Y_{l,m}^{(c)}\propto\sqrt{2}\Re(Y_l^{|m|}),\quad Y_{l,m}^{(s)}\propto\sqrt{2}\Im(Y_l^{|m|}) $$

5) 氢样径向形式：$R_{nl}(r)=\mathcal{N}e^{-\rho/2}\rho^lL_{n-l-1}^{2l+1}(\rho)$，$\rho=2Zr/(na_0)$

6) 拉盖尔三项递推：
$$ (n+1)L_{n+1}^{(\alpha)}(x)=(2n+1+\alpha-x)L_n^{(\alpha)}(x)-(n+\alpha)L_{n-1}^{(\alpha)}(x) $$

7) STO 基函数与归一化：
$$ \chi_j(r)=N_j r^{n_j^*-1}e^{-\zeta_j r},\quad N_j=\frac{(2\zeta_j)^{n_j^*+1/2}}{\sqrt{(2n_j^*)!}} $$

8) STO 线性组合：$R(r)=\sum_j c_j\chi_j(r)$

9) 径向概率密度：$P(r)=r^2|R(r)|^2$

10) CDF：$F(r)=\int_0^r P(t)dt$，逆变换采样：$F(r)=\xi,\ \xi\sim\mathcal{U}(0,1)$

11) $Z_{eff}$ 定义：$Z_{eff}(r)=Z-rV_{ee}(r)$

12) 线性插值：$V(r)\approx(1-t)V(r_k)+tV(r_{k+1})$

13) 远场边界：$V_{ee}(r)\to (N-1)/r$，$Z_{eff}(r)\to 1$

14) 累计贡献：$\langle U\rangle(<R)=\int_0^R P(r)U(r)dr$，密度：$d\langle U\rangle/dR=P(R)U(R)$

15) 局域动能恒等式（局域势本征方程假设下）：$t_{loc}=\epsilon-V_{eff}$

16) 抖动：$\delta\sim\mathcal{U}(-\Delta,\Delta)$

17) 95% 半径分位：$r_{95}=\mathrm{quantile}_{0.95}(\{r_i\})$

18) 等值面阈值分位：$\psi_{iso}=\mathrm{quantile}_{0.95}(\{|\psi(\mathbf{r}_i)|\})$，提取 $|\psi(\mathbf{r})|=\psi_{iso}$

19) 杂化对齐目标：$J(q)=\sum_k\log(\sigma_k(A(q)))$

---

## 9. 数据与公式来源清单（统一放在文末）

### 9.1 仓库内数据来源

- Koga 系列 STO 数据表：见 [slater_basis.js](slater_basis.js) 文件头注释（Koga 1999/2000 数据集，以及 `Au_R` 相关数据）。
- 电子-电子有效势查表（H–Kr）：见 [vee_cache.js](vee_cache.js) 与生成脚本 [ _experiments/generate_production_cache.py ](_experiments/generate_production_cache.py)。
- 交换项角系数表：见 [ _experiments/generate_production_cache.py ](_experiments/generate_production_cache.py) 中 `EXCHANGE_COEFFS`（角动量代数系数的工程化预制表）。

### 9.2 公式来源（标准教材/经典算法）

以下引用用于说明“采用的是标准定义”，不用于复制任何原文：

- 原子单位制、氢原子解析解、球谐函数与关联勒让德函数：常见量子力学教材（例如 Griffiths《量子力学导论》、Bransden & Joachain《Physics of Atoms and Molecules》等）。
- Slater 型轨道（STO）与归一化常数：量子化学教材（例如 Szabo & Ostlund《Modern Quantum Chemistry》）及 STO 基组文献常用定义。
- Slater–Condon/Racah 角动量代数、Gaunt 系数/Wigner $3j$ 与交换项角部分：原子结构理论与角动量代数参考书（例如 Edmonds《Angular Momentum in Quantum Mechanics》）。
- Marching Cubes：Lorensen & Cline (1987) 的经典等值面算法。

---

## 附录 A：当前版本全文保留（v3.0 物理深度解析版，原文不删减）

以下内容为重写前的文档全文保留，以满足“保留当前内容、不可删减”的要求。

# 核心算法技术文档 (v3.0 物理深度解析版)

本文档旨在提供对本项目物理引擎的**物理层级**解析。本文不仅描述算法实现，更侧重于揭示底层的物理近似、数学工具的选用依据以及系统边界。

**核心约定**：
*   所有计算均在**原子单位制 (Atomic Units, a.u.)** 下进行：$\hbar = m_e = e = 1 / 4\pi\epsilon_0 = 1$。
*   长度单位为玻尔半径 ($a_0$)，能量单位为 Hartree ($E_h$)。

**重要定位（产品与物理边界）**：
*   本项目的目标是**可视化与交互式诊断**，不是提出新的理论近似、更不是替代量化计算软件。
*   文档中出现的 $V_{ee}(r)$、$Z_{eff}(r)$、$dE/dR$、$t_{loc}(r)$ 等曲线，均应理解为在“中心场/构型平均 + STO 轨道 +（局域化/表格化的）有效势”这一实现框架下的**模型内自洽诊断量**；它们不等同于严格可观测量，也不等同于多体能量密度的唯一分解。

---

## 1. 物理模型边界与局限性 (Critical Physics Limitations)

本项目并非全功能的量子化学计算软件（如 Gaussian 或 VASP），而是一个侧重于**实时可视化**的物理引擎。必须明确以下核心近似：

### 1.1 哈密顿量近似
我们求解的是**非相对论性、定核、定态薛定谔方程**：
$$ H = -\sum_i \frac{1}{2}\nabla_i^2 - \sum_i \frac{Z}{r_i} + \sum_{i<j} \frac{1}{|r_i - r_j|} $$

*   **非相对论近似 (Non-relativistic Limit)**：
    *   轨道径向部分来自 [slater_basis.js](slater_basis.js) 中预置的 STO 展开数据（文件头标注为 *Koga (1999)/(2000)* 系列，并包含一个 `Au_R` 的特例数据集）。
    *   本引擎总体上仍处于**非相对论**框架：不求解 Dirac-HF，也不包含自旋-轨道耦合导致的精细结构分裂（如 $2p_{1/2}$ 与 $2p_{3/2}$）。
    *   **工程边界**：电子-电子有效势查表 [vee_cache.js](vee_cache.js) 当前仅覆盖 **H–Kr**（见文件头注释与生成脚本），因此与 $V_{ee}$/$Z_{eff}$/能量分解相关的曲线在更重元素上不具备同等物理含义。
*   **单电子近似 (Single Particle Approximation)**：
    *   我们展示的是 **Hartree-Fock (HF)** 意义下的单电子轨道。这意味着忽略了**电子关联能 (Electron Correlation Energy)**。因此，电子云形状是“平均场”下的结果，无法体现多体波函数的纠缠特性。

### 1.2 冻结核近似 (Frozen Core Approximation)
`slater_basis.js` 中的基组参数是针对**基态 (Ground State)** 原子优化的。
*   当用户可视化激发态（例如碳原子的 $5d$ 轨道）时，我们使用的是基态原子的内层电子屏蔽参数。
*   **误差分析**：这会引入显著误差，因为激发态电子感受到的有效核电荷 $Z_{eff}$ 与基态完全不同。可视化的高激发态仅具有定性的拓扑正确性，能量和半径可能有较大偏差。

---

## 2. 氢原子：解析解的数值稳定性

文件: `physics.js`

对于 $Z=1$ 系统，我们求解：
$$ -\frac{1}{2}\nabla^2\psi - \frac{1}{r}\psi = E\psi $$

### 2.1 广义拉盖尔多项式 (Generalized Laguerre Polynomials)
径向波函数 $R_{nl}(r) \propto e^{-\rho/2} \rho^l L_{n-l-1}^{2l+1}(\rho)$。

**物理惯例警告**：
数学界（如 Mathematica）与物理界（量子力学教材）对拉盖尔多项式的定义可能相差一个阶乘因子 $(n+l)!$。
*   本项目采用了**物理学家惯例**，确保归一化条件 $\int_0^\infty R_{nl}^2 r^2 dr = 1$ 严格成立。
*   **数值陷阱**：直接计算 $L_n^\alpha(x)$ 的封闭形式会在 $n>20$ 时溢出。我们实现了 **三项递推 (Three-Term Recurrence)** 算法：
    $$ (n+1)L_{n+1}^{(\alpha)}(x) = (2n+1+\alpha-x)L_n^{(\alpha)}(x) - (n+\alpha)L_{n-1}^{(\alpha)}(x) $$
    该算法在 $n=1\sim 100$ 范围内具有优异的数值稳定性。

### 2.2 实球谐函数 (Real Spherical Harmonics)
为了可视化轨道波瓣方向，我们在内部采用复球谐函数并线性组合得到实球谐函数。
**Condon-Shortley 相位**：
代码显式引入了 $(-1)^m$ 相位因子，复球谐的规范形式为：
$$ Y_l^m \propto (-1)^m \sqrt{\frac{(2l+1)}{4\pi} \frac{(l-m)!}{(l+m)!}} P_l^m(\cos\theta) e^{im\phi} $$
实球谐通过 $Y_l^{m}$ 与 $Y_l^{-m}$ 的线性组合得到，这保证了基组展开时的正负号与常见量子化学软件/教材的约定一致，避免了 $p_x, p_y$ 混合时的符号混乱。


---

## 3. 多电子原子：Koga-STO 轨道数据与中心场诊断量（VeeCache 仅 H–Kr）

文件: `slater_basis.js`, `physics.js`

对于多电子原子，解析解不存在。本项目不在浏览器端求解 Roothaan-HF 迭代，而是使用预置的高精度 STO 展开来近似表示 HF 轨道径向形状与轨道能量，并在此基础上构造用于可视化的“中心场有效势”诊断量。

### 3.1 基组选择：Koga (1999)
代码库内置了 [slater_basis.js](slater_basis.js) 中的 STO 展开数据（文件头标注为 Koga (1999)/(2000) 系列）。本项目在文档层面不试图“重新推导/提出”这些基组或系数，仅将其作为既有数据源用于可视化。
*   **基函数类型**：Slater-Type Orbitals (STOs)。
    $$ \chi(r) = N r^{n-1} e^{-\zeta r} $$
*   **非相对论性**：该类数据集一般基于非相对论哈密顿量优化；对于重元素，相对论效应会改变内层轨道与能级（这里不做系统性修正）。
*   **展开式**：每个原子轨道表示为多个 Slater 函数的线性组合 (LCAO)：
    $$ R_{nl}(r) = \sum_{j=1}^{M} c_j \chi_j(r; n_j, \zeta_j) $$
    其中 $M,\{c_j,n_j,\zeta_j\}$ 由数据表给定。

### 3.2 构型平均 (Configuration Averaging)
对于开壳层原子（Open-shell atoms，如 C, O），Hartree-Fock 轨道依赖于具体的 $L,S$ 耦合态。
*   **近似**：本项目使用了**构型平均 (Average of Configuration)** 的轨道参数。即我们不区分 $2p_x, 2p_y, 2p_z$ 在晶体场下的微小差异，认为它们是简并的。

---

## 4. 有效核电荷 $Z_{eff}$：从 $V_{ee}$ 反推

文件: `physics.js` -> `calculateZeff`

教科书通常只给出 Slater Rules 的常数 $Z_{eff}$。本项目实现了**空间分辨 (Spatially Resolved)** 的 $Z_{eff}(r)$。

### 4.1 理论定义
基于重构核势的思想，可以定义一个等效的核电荷 $Z_{eff}(r)$，使得该点的总势能 $V(r)$ 看起来像是一个类氢势：
$$ V(r) = -\frac{Z_{eff}(r)}{r} $$
代入总势能 $V(r) = -\frac{Z}{r} + V_{ee}(r)$：
$$ -\frac{Z_{eff}(r)}{r} = -\frac{Z}{r} + V_{ee}(r) \implies Z_{eff}(r) = Z - r \cdot V_{ee}(r) $$

### 4.2 电子-电子排斥势 $V_{ee}(r)$
这是用于可视化的、与轨道相关的“中心场有效势” $V_{ee}^{(i)}(r)$：它在生成端以 STO 解析积分形式构造，包含直接库仑项（Hartree）以及按角动量代数预制系数加权的交换项（Fock 的局域化/球平均表达）。
*   **预计算 (Pre-computation)**：由于计算 $V_{ee}$ 需要对所有被占轨道进行双重积分，极为耗时。我们采用了 Look-up Table (LUT) 策略，存储在 `VeeCache` 中。
*   **插值**：运行时通过二分查找 + 线性插值获取 $V_{ee}$。
*   **物理意义**：当 $r \to 0$，电子穿透到核附近，$r V_{ee} \to 0$，故 $Z_{eff} \to Z$（裸核）；当 $r \to \infty$，电子在该壳层之外，$Z_{eff} \to Z - (N-1)$（完全屏蔽）。

**关于“系数从哪来”**：生成脚本 [_experiments/generate_production_cache.py](_experiments/generate_production_cache.py) 内的 `EXCHANGE_COEFFS` 是一个 $(l_{target}, l_{source}) \mapsto \{(k,\mathrm{coeff})\}$ 的预制表，用于将交换项写成有限个径向核 $Y_k(r)$ 的线性组合。其角色等价于原子结构理论中 Slater–Condon/Racah 角动量代数（Gaunt/Wigner $3j$ 等）给出的“角部分系数”。本项目只是在工程上复用这类标准结果来生成查表数据，并不宣称提出新的交换近似或新系数。

**实现范围提醒**：当前 [vee_cache.js](vee_cache.js) 文件头明确标注仅包含 **H–Kr** 的预计算数据；当所选原子/轨道不存在对应缓存时，[physics.js](physics.js) 会回退为 $V_{ee}=0$ 或 $Z_{eff}=Z$ 的退化行为，因此此时相关曲线不应被作物理解释。

### 4.3 自旋统计缩放因子（构型平均）

在生成 [vee_cache.js](vee_cache.js) 的过程中（见 [_experiments/generate_production_cache.py](_experiments/generate_production_cache.py)），为了避免把闭壳层错误当作“全平行自旋”从而导致交换项过大，我们对交换项引入一个基于洪特规则的**自旋统计缩放因子 $\chi$**。

这不是新理论：它是“构型平均/自旋统计权重”的一种工程化实现。其有效性依赖于“采用高自旋基态占据分配”的假设；当用户关注具体 $LS$ 项、激发态或离子态时，该近似不一定成立。

*   **物理背景**：交换作用仅发生在自旋平行的电子之间。对于一个占据数为 $N$、容量为 $C$ 的亚层，自旋平行的统计概率取决于自旋组态 $N_{\uparrow}, N_{\downarrow}$。
*   **修正公式**：
    $$ V_{exch}^{corrected} = \chi \cdot V_{exch}^{geo} $$
    其中几何交换势 $V_{exch}^{geo}$ 对应“全自旋平行”的情况。自旋因子 $\chi$ 定义为：
    $$ \chi = \frac{N_{\uparrow}^2 + N_{\downarrow}^2}{N^2} $$
*   **典型案例验证**：
    *   **He ($1s^2$)**: $N=2, N_{\uparrow}=1, N_{\downarrow}=1 \implies \chi = 0.5$。用于避免闭壳层被错误当作“全平行自旋”而导致交换项过大。
    *   **O ($2p^4$)**: $N=4, N_{\uparrow}=3, N_{\downarrow}=1 \implies \chi = 0.625$。较“开壳层因子=1.0”的粗略近似更细致。



---

## 5. 能量图表计算：径向累计期望值与局域动能（local kinetic energy）

这一节回答两个“物理严谨性”问题：

1) 势能曲线 $E(R)$ / 势能密度 $dE/dR$ 在本项目里到底是什么量？

2) 代码里曾被叫作“动能重构/逆向重构”的部分是什么？（结论：它是**由一电子本征方程定义的局域动能**，不是新物理。）

### 5.1 势能曲线的严格定义（不是“能量随半径变化”）

文件: `physics.js` -> `calculateCumulativeOrbitalEnergy`

对任意一个（球平均的）轨道 $i$，其径向概率密度定义为：
$$ P_i(r) = r^2 |R_i(r)|^2, \qquad \int_0^{\infty} P_i(r)\,dr = 1 $$

对任何**局域的一电子径向算符** $U(r)$（例如 $-Z/r$ 或某个局域有效势），其期望值可写成：
$$ \langle U \rangle_i = \int_0^{\infty} P_i(r)\,U(r)\,dr $$

本项目绘制的“能量曲线”在数学上是一个**径向累计期望值/累计贡献函数**：
$$ \langle U \rangle_i(<R) = \int_0^{R} P_i(r)\,U(r)\,dr $$
因此其导数是严格的：
$$ \frac{d}{dR}\langle U \rangle_i(<R) = P_i(R)\,U(R) $$

在实现中：

*   $U_{nuc}(r)=-Z/r$：该项使用不完全伽马函数对 STO 组合做解析积分，并返回核势贡献密度（见 [physics.js](physics.js) 中 `dVnucDr`）。
*   $U_{ee}^{(i)}(r)=V_{ee}^{(i)}(r)$：该项来自 `VeeCache`，并以 $P_i(r)\,V_{ee}^{(i)}(r)$ 的形式进入 $dE/dR$（见 `dEdr` 的构造）。

重要：这里的 $V_{ee}^{(i)}(r)$ 是**与轨道有关的中心场有效势**（直接库仑项 + 交换项的径向表达/近似），用于可视化一电子方程意义下的 $V_{eff}$。它不是多体体系中严格意义的“电子-电子势能密度”。

**关于“potential”图名的澄清**：当前 UI 中 `potential`/`dEdr` 使用的是上述 $U(r)$ 的累计与密度；其中非氢原子分支在 [physics.js](physics.js) 内**有意省略**对 $T$ 的直接数值计算（避免 cusp 处二阶导数噪声），因此 `E(R)` 默认反映的是“有效势贡献”的累计，而不是严格的总能量。

### 5.2 远场边界条件（必须显式满足）

对中性原子（总电子数 $N=Z$），当 $r\to\infty$ 且电子处在原子外部时：

*   电子-电子项的长程主导是库仑（Hartree）项，交换作用由于轨道重叠迅速衰减，为短程。
*   因此必须有
    $$ V_{ee}^{(i)}(r) \xrightarrow[r\to\infty]{} \frac{N-1}{r}, \qquad Z_{eff}(r)=Z-rV_{ee}(r)\xrightarrow[r\to\infty]{} 1 $$

本项目在两处显式保证该条件：

*   **生成端**：`_experiments/generate_production_cache.py` 在缓存尾部强制 $V_{ee}(r)=(N-1)/r$（最后一段采样点），防止开壳层 $p/d$ 的交换“构型平均系数”在数值尾部留下非物理残余。
*   **运行端**：`physics.js` 对 $r > r_{max}$ 使用 $1/r$ 外推（由最后一个缓存点估计系数），避免“常数夹断”产生非物理长程尾巴。

这不是“拍脑袋创新”，而是对已知物理边界条件的实现约束。

### 5.3 局域动能（为何代码里曾叫“重构”）

文件: `physics.js` -> `calculateLocalKineticEnergy`（旧名 `calculateReconstructedKineticEnergy`）

如果某轨道满足一个**局域势的一电子本征方程**：
$$ \left[-\frac12\nabla^2 + V_{eff}(r)\right]\phi_i = \epsilon_i\phi_i $$
则有点态恒等式（把等式除以 $\phi_i$）：
$$ t_{loc,i}(r) \equiv \frac{-\frac12\nabla^2\phi_i}{\phi_i} = \epsilon_i - V_{eff}(r) $$

这就是量子化学/DFT 语境中常见的 **local kinetic energy / local energy** 概念（注意：动能密度的定义并不唯一，但这个 $\epsilon - V_{eff}$ 在“局域势本征方程”假设下是自洽的）。

本项目用 $V_{eff}(r)=V_{nuc}(r)+V_{ee}^{(i)}(r)$ 构造 $t_{loc}$，并绘制两条径向密度：

*   $\epsilon_i\,P_i(r)$：轨道能量密度（对 $r$ 积分得到 $\epsilon_i$）。
*   $t_{loc,i}(r)\,P_i(r)$：局域动能密度（对 $r$ 积分得到一个“模型内”的动能期望）。

为什么这么做：它避免对 STO 波函数做二阶数值微分（Cusp 处会非常不稳定），并且在同一近似框架下保持 $t_{loc}(r)+V_{eff}(r)=\epsilon$ 的一致性。

### 5.4 严格性边界（哪些是模型假设，不要误读）

*   Hartree-Fock 的交换算符严格来说是**非局域**的；把它表示成一个对每个轨道都可乘的 $V_{ee}^{(i)}(r)$ 属于中心场/构型平均框架下的表达与近似。
*   因此 $t_{loc}(r)$ 与这些能量密度曲线属于“**该模型内部自洽的诊断量**”，不能直接等同于实验可观测量或严格的多体能量密度分解。
*   你要拿它做“物理结论”，应当把它限定为：在 Koga STO + 构型平均 +（局域化的）交换处理这套近似内的解释。

---

## 6. 核心采样算法 (Sampling Algorithms)

这是本项目实现“百万级点云”实时生成的关键。

### 6.1 径向：数值逆 CDF (Inverse CDF)
对于径向分布函数 $P(r) = r^2 |R(r)|^2$，我们避免使用效率较低的 Metropolis-Hastings。
*   **实现**：
    1.  预计算 CDF: $F(r) = \int_0^r P(t) dt$（默认 2000 个采样点，梯形积分；见 [physics.js](physics.js) 中 `buildRadialCDF`）。
    2.  逆变换：生成均匀随机数 $\xi \in [0,1]$，求解 $F(r) = \xi$。
*   **精度说明**：由于 CDF 是单调递增的，二分查找 (Binary Search) 保证了 $O(\log N)$ 的时间复杂度。由于 $P(r)$ 是平滑函数，线性插值在可视化尺度上通常足够。
*   **无拒绝 (Rejection-Free)**：该步骤无额外拒绝环节，但仍属于数值近似采样。

### 6.2 角向：局部拒绝采样 (Local Rejection)
角向分布 $|Y_{lm}(\theta, \phi)|^2$ 比较复杂，且难以构建 2D 逆 CDF。
*   **策略**：我们在单位球面上均匀布点，但按概率拒绝。
*   **优化**：由于 $Y_{lm}$ 的最大值 $M_{max}$ 已知（或容易估算），接受率通常在 30%-50% 之间，远高于全局 3D 采样的 <1%。

### 6.3 摩尔纹 (Moiré Patterns) 与抖动 (Dithering)
*   **现象**：由于逆 CDF 表的离散性（默认 2000 点），在生成数百万个粒子时，径向会出现肉眼可见的同心圆环（离散化伪影）。
*   **解决方案**：在最终坐标上叠加小幅随机抖动（Dithering），以削弱离散采样造成的视觉纹理。当前实现为对 $(x,y,z)$ 各分量加入均匀抖动 `±0.005 a0`（即范围宽度约 $0.01 a_0$；见 [sampling.js](sampling.js)）。

---

## 7. 等值面提取：Marching Cubes 算法

文件: `marching_cubes.js`, `isosurface-worker.js`

### 7.1 并行计算架构
为了不阻塞 UI 线程（保持拖拽旋转流畅），等值面计算完全卸载到 Web Worker。
*   **Transferable Objects**：主线程与 Worker 之间传输大量顶点数据时，使用 `ArrayBuffer` 的所有权转移（Transfer），实现零拷贝 (Zero-copy) 通信。

### 7.2 连通域分离 (Connected Component Analysis)
这是正确渲染 $p, d, f$ 轨道相位的关键。
*   **问题**：$2p_z$ 轨道由正瓣（上）和负瓣（下）组成。Marching Cubes 输出的是一堆无序三角形。
*   **算法**：
    1.  建立三角形邻接图。
    2.  执行广度优先搜索 (BFS) 或并查集 (Union-Find) 将网格分离为独立的连通子网格。
    3.  对每个子网格，计算其几何中心的波函数值 $\psi$。
    4.  若 $\psi > 0$ 染蓝色，$\psi < 0$ 染红色。
*   **鲁棒性**：该算法能自动处理任意复杂的轨道形状（如 $d_{z^2}$ 的环形瓣），无需硬编码几何规则。

---

## 8. 杂化轨道自适应对齐 (Adaptive Hybridization Alignment)

为了解决随机生成的几何构型 (Thomson Points) 与各向异性基组（如 $d_{z^2}$）的不对齐问题，我们引入了**四元数旋转优化器**。

### 8.1 核心原理
算法寻找一个最佳旋转 $q$，使得基组在几何构型上的投影矩阵 $A$ 的体积最大化：
$$ J(q) = \sum \log(\sigma_k(A(q))) $$
其中 $\sigma_k$ 是投影矩阵的奇异值。最大化这些奇异值的对数和等价于最大化投影体积（det值），确保几何构型与基组最匹配。

### 8.2 优势
- **sp3d (TBP)**: 自动将轴向对齐到 Z 轴，使 $d_{z^2}$ 正确参与轴向成键。
- **sp3d2 (Oct)**: 自动对齐 XYZ 轴，促成轨道成分的清晰分离。
- **通用性**: 支持任意 $s, p, d$ 轨道的组合，无需硬编码。
- **鲁棒性**: 能够检测并处理“线性相关”或“不兼容”的杂化请求（此时体积 $\det \to 0$）。

### 8.3 实现细节
- **优化算法**: Random Restart Hill Climbing (20 次重启，每次 50 步迭代)。
- **数学库**: 内置微型四元数库 (`physics-core.js`)。
- **性能**: 计算量极小，通常在 <5ms 内收敛。

---
*文档生成规范：基于 v3.0 物理引擎代码审计*
