# 核心算法与物理口径说明

本文件是“信任工程（Trust Engineering）”的一部分：它只描述**当前仓库代码**在浏览器端实际计算与展示的量，并给出对应的数学定义、近似边界、以及仓库内的数据/公式来源，便于独立复核。

项目的核心定位与设计原则参照 [publication_prep/paper_cn_v3.md](publication_prep/paper_cn_v3.md)（尤其是“信任工程”小节）：用户在屏幕上看到的一切都应能被追溯到明确的公式、明确的数据表、或明确的代码路径。

---

## 1. 约定、定位与可验证性

### 1.1 单位制与符号

本项目默认使用原子单位制（Atomic Units, a.u.）：

$$ \hbar = m_e = e = \frac{1}{4\pi\epsilon_0} = 1 $$

长度单位为玻尔半径 $a_0$，能量单位为 Hartree $E_h$。

### 1.2 核心定位（paperCNV3 口径）

本项目作为“原子轨道虚拟实验室”运行，其核心定位是：

> **一个基于「Koga & Thakkar 多指数 STO 近似」的虚拟实验室，用户在屏幕上看到的一切，都是该近似框架下的数值结果。**

在当前实现中：

- 对氢原子：使用薛定谔方程解析解（径向函数 + 球谐函数）。
- 对多电子原子（H–Kr）：使用 Koga 系列多指数 STO 数据表来近似非相对论 Hartree–Fock 极限下的轨道形状与轨道本征值（见 [slater_basis.js](slater_basis.js)）。

### 1.3 信任工程：可验证性护栏（理论曲线 vs 采样）

UI 同时展示（或可生成）两类信息：

- **理论曲线**：由波函数定义直接给出，例如径向概率密度 $P(r)=r^2|R(r)|^2$。
- **采样直方图**：由蒙特卡洛点云统计得到，并随样本数增加逐步收敛到理论曲线（大数定律层面的可验证性）。

这条“理论—采样”一致性链路是当前工程质量护栏：当图形看起来“奇怪”时，用户可通过是否收敛于理论曲线来判断它是物理效应（在该近似框架内）还是实现问题。

### 1.4 覆盖范围与回退行为（透明说明）

与 $V_{ee}$/$Z_{eff}$/能量分解相关的查表数据目前仅覆盖 H–Kr（见 [vee_cache.js](vee_cache.js) 及生成脚本 [_experiments/generate_production_cache.py](_experiments/generate_production_cache.py)）。

当所选原子/轨道不存在对应缓存时，[physics.js](physics.js) 会回退到退化行为（例如 $V_{ee}=0$ 或 $Z_{eff}=Z$ 的简化），此时相关曲线不应被作物理解释。

---

## 2. 物理模型边界与局限性（先写清楚“没做什么”）

### 2.1 哈密顿量层级与近似

出发点是非相对论、定核、定态薛定谔框架：

$$ H = -\sum_i \frac{1}{2}\nabla_i^2 - \sum_i \frac{Z}{r_i} + \sum_{i<j} \frac{1}{|\mathbf{r}_i - \mathbf{r}_j|} $$

明确的近似边界：

- 非相对论：不求解 Dirac-HF，不含自旋-轨道耦合导致的精细结构分裂（如 $2p_{1/2}$ 与 $2p_{3/2}$）。
- 单电子轨道口径（HF 轨道）：忽略电子关联能（correlation）；电子云形状是平均场结果。

### 2.2 冻结核/基态参数外推（激发态的定性使用）

[slater_basis.js](slater_basis.js) 的基组参数针对基态原子优化。

- 当用户可视化激发态（例如 C 的 $5d$）时，仍沿用基态的内层屏蔽与径向参数。
- 这会引入显著误差：激发态电子感受到的 $Z_{eff}$ 与基态不同。此时可视化更接近“拓扑/定性正确”，能量与半径可能偏差很大。

---

## 3. 波函数与概率密度：项目里采用的严格数学定义

### 3.1 坐标与分离变量

采用球坐标：

$$ r = \sqrt{x^2+y^2+z^2},\quad \theta = \arccos\left(\frac{z}{r}\right),\quad \phi = \operatorname{atan2}(y,x) $$

单电子轨道常写作径向与角向的乘积（或线性组合）：

$$ \psi(\mathbf{r}) = R(r)\,Y(\theta,\phi) $$

### 3.2 径向概率密度（用于采样与能量累计）

对球平均口径的轨道 $i$，径向概率密度定义为：

$$ P_i(r) = r^2 |R_i(r)|^2,\qquad \int_0^{\infty} P_i(r)\,dr = 1 $$

---

## 4. 球谐函数：复球谐与实球谐（化学惯例一致性）

实现主要在 [physics-core.js](physics-core.js) 的 `Ylm_complex` 与 `realYlm_value`。

### 4.1 复球谐（含 Condon–Shortley 相位）

$$ Y_l^m(\theta,\phi) = (-1)^m\,\sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}\,P_l^m(\cos\theta)\,e^{im\phi} $$

### 4.2 实球谐（cos/sin 线性组合）

实球谐通过 $Y_l^{m}$ 与 $Y_l^{-m}$ 线性组合得到，并区分 cos 型与 sin 型：

$$ Y_{l,m}^{(c)} \propto \sqrt{2}\,\Re(Y_l^{|m|}),\qquad Y_{l,m}^{(s)} \propto \sqrt{2}\,\Im(Y_l^{|m|}) $$

采用该组合的意图是：在 $p_x,p_y$ 等实轨道上与常见量子化学软件/教材的约定一致，避免符号混乱。

---

## 5. 径向函数：氢解析式与 Koga-STO 表

### 5.1 氢样原子的解析形式（强调数值稳定实现）

氢样径向函数通常写为：

$$ R_{nl}(r) = \mathcal{N}_{nl}\,e^{-\rho/2}\,\rho^l\,L_{n-l-1}^{2l+1}(\rho),\qquad \rho=\frac{2Zr}{n a_0} $$

为避免高阶溢出，广义拉盖尔多项式使用三项递推：

$$ (n+1)L_{n+1}^{(\alpha)}(x) = (2n+1+\alpha-x)L_n^{(\alpha)}(x) - (n+\alpha)L_{n-1}^{(\alpha)}(x) $$

### 5.2 多电子原子：Koga-STO 展开

项目不在浏览器端求解 Roothaan-HF 自洽迭代，而使用预置 STO 展开近似 HF 轨道径向形状与轨道能量（数据见 [slater_basis.js](slater_basis.js)，文件头注明 Koga (1999)/(2000) 系列）。

STO 基函数：

$$ \chi(r) = N\,r^{n-1} e^{-\zeta r} $$

在实现（见 [physics-core.js](physics-core.js) 的 `slaterRadialR`）中，STO 项采用 $n^*$ 与系数：

$$ \chi_j(r) = N_j\,r^{n_j^*-1} e^{-\zeta_j r} $$

并使用归一化常数：

$$ N_j = \frac{(2\zeta_j)^{n_j^*+1/2}}{\sqrt{(2n_j^*)!}} $$

轨道径向函数是线性组合：

$$ R_{nl}(r) = \sum_{j=1}^{M} c_j\,\chi_j(r) $$

### 5.3 构型平均（开壳层简并口径）

对开壳层原子（如 C, O），HF 轨道严格依赖具体 $L,S$ 耦合态。项目采用构型平均（Average of Configuration）的口径：不区分 $2p_x,2p_y,2p_z$ 在晶体场下的微小差异，按简并处理。

---

## 6. $V_{ee}(r)$ 与 $Z_{eff}(r)$：定义、生成、插值与边界条件

### 6.1 从总势能反推 $Z_{eff}(r)$

定义一个等效核电荷 $Z_{eff}(r)$ 使总势能呈现类氢形式：

$$ V(r) = -\frac{Z_{eff}(r)}{r} $$

而项目取

$$ V(r) = -\frac{Z}{r} + V_{ee}(r) $$

因此

$$ Z_{eff}(r) = Z - r\,V_{ee}(r) $$

实现入口为 [physics.js](physics.js) 的 `calculateZeff`。

### 6.2 $V_{ee}^{(i)}(r)$ 的意义边界

$V_{ee}^{(i)}(r)$ 是用于可视化的、与轨道相关的中心场有效势：包含直接库仑项（Hartree）以及交换项的径向表达（以工程化方式球平均/局域化并查表）。

强调：HF 的交换算符严格非局域。将其表示成对每个轨道都可乘的 $V_{ee}^{(i)}(r)$ 属于中心场/构型平均框架下的表达与近似。

### 6.3 交换项角系数：来源说明（避免误解为“自创系数”）

生成脚本 [_experiments/generate_production_cache.py](_experiments/generate_production_cache.py) 内的 `EXCHANGE_COEFFS` 是一个 $(l_{target}, l_{source}) \mapsto \{(k,\mathrm{coeff})\}$ 的预制表，用于将交换项写成有限个径向核 $Y_k(r)$ 的线性组合。

它的角色等价于原子结构理论中 Slater–Condon/Racah 角动量代数（例如 Gaunt 系数、Wigner $3j$ 等）给出的角部分系数。本项目在工程上复用这类标准结果来生成缓存数据，不宣称提出新的交换近似或新系数。

### 6.4 运行时插值

若缓存提供在离散网格 $\{r_k\}$ 上的 $V_{ee}(r_k)$，运行时使用二分查找定位区间并线性插值：

$$ V_{ee}(r) \approx (1-t)V_{ee}(r_k) + tV_{ee}(r_{k+1}),\quad t=\frac{r-r_k}{r_{k+1}-r_k} $$

### 6.5 远场边界条件（必须显式满足）

对中性原子 $N=Z$，当 $r\to\infty$ 且电子在原子外部时，应满足：

$$ V_{ee}^{(i)}(r) \xrightarrow[r\to\infty]{} \frac{N-1}{r},\qquad Z_{eff}(r) = Z - rV_{ee}(r) \xrightarrow[r\to\infty]{} 1 $$

项目在两处显式保证该条件：

- 生成端：在缓存尾部强制 $V_{ee}(r)=(N-1)/r$，避免数值尾部残余。
- 运行端：对 $r>r_{max}$ 使用 $1/r$ 外推，避免常数夹断产生非物理长程尾巴。

### 6.6 自旋统计缩放因子 $\chi$（构型平均下的工程化处理）

为避免把闭壳层错误当作全平行自旋导致交换项过大，在生成 [vee_cache.js](vee_cache.js) 的过程中引入自旋统计缩放因子 $\chi$（见 [_experiments/generate_production_cache.py](_experiments/generate_production_cache.py)）。

交换项修正写作：

$$ V_{exch}^{corrected} = \chi\,V_{exch}^{geo} $$

其中 $V_{exch}^{geo}$ 对应“全自旋平行”的几何交换势。缩放因子定义为：

$$ \chi = \frac{N_{\uparrow}^2 + N_{\downarrow}^2}{N^2} $$

典型例子：

- He ($1s^2$)：$N=2, N_{\uparrow}=1, N_{\downarrow}=1 \Rightarrow \chi=0.5$
- O ($2p^4$)：$N=4, N_{\uparrow}=3, N_{\downarrow}=1 \Rightarrow \chi=0.625$


---

## 7. 能量相关图表：定义与代码映射

本节只描述当前代码给 UI 的三个能量相关输出：`potential`、`dEdr`、`localEnergy`。

### 7.1 `dEdr`（势能密度）与 `potential`（其径向累计）

能量相关曲线由 [physics.js](physics.js) 的 `calculateCumulativeOrbitalEnergy` 提供。当前实现中：

- 返回值里的 `T` 恒为 `null`（不输出通过 $\nabla^2$ 直接计算的动能积分项）。
- 返回值里的 `dEdr` 对应“势能密度”（更准确地说，是 $dV/dr$ 的离散化）。

在数学口径上，`dEdr` 可理解为

$$ \frac{dV}{dr} = \frac{dV_{nuc}}{dr} + P(r)\,V_{ee}(r) $$

其中

$$ P(r)=r^2|R(r)|^2,\qquad V_{nuc}(r)=-\frac{Z}{r} $$

而 $V_{ee}(r)$ 来自 [vee_cache.js](vee_cache.js) 的 `VeeCache`（轨道相关查表 + 线性插值 + $1/r$ 远场外推）。

`calculateCumulativeOrbitalEnergy` 在氢原子与非氢原子路径下返回值的含义略有不同（这是代码层面的事实，应按返回值解释 UI）：

- 对氢原子：`E(R)=\epsilon\,F(R)`（$F(R)$ 为径向 CDF），`dEdr` 返回的是核势密度 `dVnucDr`（由于 $V_{ee}=0$，它也等于 $dV/dr$），并同时返回核势累计 `V(R)`（键名为 `V`）。
- 对非氢原子：`dEdr` 是势能密度（核势密度 + $P(r)V_{ee}(r)$），`E(R)` 由对 `dEdr` 做梯形积分得到；同时返回核势累计 `V(R)`（键名为 `V`）。

### 7.2 `localEnergy`：两条径向能量密度曲线

`localEnergy` 图由 [orbital.js](orbital.js) 预计算并交给 [data_panel.js](data_panel.js) 渲染。它包含两条“径向能量密度”曲线：

- `epsDensity(r) = \epsilon\,P(r)`
- `Tdensity(r) = t_{loc}(r)\,P(r)`

其中 $\epsilon$ 来自 [slater_basis.js](slater_basis.js) 的 `SlaterBasis[atomType].energies[orbitalKey]`（多选/叠加时取选中轨道的平均值）。

`t_{loc}(r)` 由 [physics.js](physics.js) 的 `calculateReconstructedKineticEnergy` / `calculateLocalKineticEnergy` 计算，其定义是

$$ t_{loc}(r)=\epsilon - V_{nuc}(r) - V_{ee}(r) $$

并使用

$$ V_{nuc}(r)=-\frac{Z}{r} $$

以及从 `VeeCache` 插值得到的 $V_{ee}(r)$。

这一定义对应一电子局域势本征方程的恒等式形式

$$ \left[-\frac12\nabla^2 + V_{eff}(r)\right]\phi = \epsilon\phi\quad \Rightarrow\quad \frac{-\frac12\nabla^2\phi}{\phi}=\epsilon - V_{eff}(r) $$

在当前实现中取 $V_{eff}(r)=V_{nuc}(r)+V_{ee}(r)$。该量用于可视化/诊断：它由“轨道本征值 $\epsilon$ + 局域化的有效势 $V_{eff}$”共同确定。

---

## 8. 采样算法（点云生成）

### 8.1 径向：逆 CDF（无拒绝采样）

目标分布 $P(r)=r^2|R(r)|^2$。构造

$$ F(r)=\int_0^r P(t)\,dt $$

生成 $\xi\in[0,1]$ 并求解 $F(r)=\xi$ 得到样本 $r$。

实现见 [physics.js](physics.js) 的 `buildRadialCDF`（默认 2000 点，梯形积分），用二分查找实现 $O(\log N)$ 反演，并线性插值。

### 8.2 角向：局部拒绝采样

角向分布为 $|Y_{lm}(\theta,\phi)|^2$。在单位球面上均匀采样方向并按概率拒绝。

由于 $|Y_{lm}|^2$ 的最大值可估，接受率通常在 30%–50% 左右，远高于全局 3D 拒绝采样。

### 8.3 摩尔纹与抖动

由于逆 CDF 表离散，生成大量粒子时可能出现同心圆伪影。项目在最终坐标上叠加小幅随机抖动以削弱视觉纹理。

抖动实现：对 $(x,y,z)$ 各分量加入均匀抖动 $\pm 0.005\,a_0$（实现见 [sampling.js](sampling.js)）。

---

## 9. 等值面提取与显示：Marching Cubes + Worker + 连通域分离

相关实现见 [marching_cubes.js](marching_cubes.js)、[isosurface-worker.js](isosurface-worker.js)、以及 [visualization.js](visualization.js) 中的调用封装。

### 9.1 并行计算架构（Worker）

等值面计算放入 Web Worker，主线程与 Worker 之间使用 Transferable 的 `ArrayBuffer` 转移所有权以降低拷贝开销。

### 9.2 连通域分离与相位着色

Marching Cubes 输出无序三角形后：

1. 建立三角形邻接图
2. BFS 或并查集分离连通子网格
3. 计算每个子网格几何中心的 $\psi$ 符号
4. $\psi>0$ 与 $\psi<0$ 使用不同颜色

该方法能处理复杂轨道形状（如 $d_{z^2}$ 的环形瓣），无需硬编码几何规则。

---

## 10. 杂化轨道自适应对齐：四元数旋转优化

为解决几何构型（Thomson 点）与各向异性基组不对齐，引入四元数旋转优化。

目标函数：

$$ J(q) = \sum \log\left(\sigma_k(A(q))\right) $$

其中 $\sigma_k$ 为投影矩阵奇异值。实现使用 [physics-core.js](physics-core.js) 内的四元数/线代工具。

实现参数：

- Random Restart Hill Climbing（20 次重启，每次 50 步迭代）
- 计算量通常 <5ms

---

## 11. 公式总表（项目中实际出现/使用的计算公式）

下面以索引形式汇总本文涉及的全部关键公式，便于核对实现：

1) 原子单位制：$\hbar = m_e = e = 1/(4\pi\epsilon_0) = 1$

2) 球坐标：$r = \sqrt{x^2+y^2+z^2}$，$\theta = \arccos(z/r)$，$\phi = \operatorname{atan2}(y,x)$

3) 复球谐（含 Condon–Shortley）：

$$ Y_l^m = (-1)^m\sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}P_l^m(\cos\theta)e^{im\phi} $$

4) 实球谐线性组合：

$$ Y_{l,m}^{(c)} \propto \sqrt{2}\Re(Y_l^{|m|}),\quad Y_{l,m}^{(s)} \propto \sqrt{2}\Im(Y_l^{|m|}) $$

5) 氢样径向形式：

$$ R_{nl}(r)=\mathcal{N}e^{-\rho/2}\rho^lL_{n-l-1}^{2l+1}(\rho),\quad \rho=\frac{2Zr}{na_0} $$

6) 拉盖尔三项递推：

$$ (n+1)L_{n+1}^{(\alpha)}(x)=(2n+1+\alpha-x)L_n^{(\alpha)}(x)-(n+\alpha)L_{n-1}^{(\alpha)}(x) $$

7) STO 基函数：$\chi(r)=Nr^{n-1}e^{-\zeta r}$

8) STO（含 $n^*$）与归一化：

$$ \chi_j(r)=N_j r^{n_j^*-1}e^{-\zeta_j r},\quad N_j = \frac{(2\zeta_j)^{n_j^*+1/2}}{\sqrt{(2n_j^*)!}} $$

9) STO 线性组合：$R(r)=\sum_j c_j\chi_j(r)$

10) 径向概率密度：$P(r)=r^2|R(r)|^2$，$\int_0^\infty P(r)dr=1$

11) CDF：$F(r)=\int_0^r P(t)dt$；逆变换：$F(r)=\xi,\ \xi\sim\mathcal{U}(0,1)$

12) $Z_{eff}$ 定义：$Z_{eff}(r)=Z-rV_{ee}(r)$

13) 线性插值：$V(r)\approx(1-t)V(r_k)+tV(r_{k+1})$

14) 远场边界：$V_{ee}(r)\to(N-1)/r$，$Z_{eff}(r)\to 1$

15) 势能密度（代码中的 `dEdr`）：

$$ \frac{dV}{dr} = \frac{dV_{nuc}}{dr} + P(r)\,V_{ee}(r) $$

16) 局域动能（代码中的 `calculateReconstructedKineticEnergy` / `calculateLocalKineticEnergy`）：

$$ t_{loc}(r)=\epsilon - V_{nuc}(r) - V_{ee}(r),\quad V_{nuc}(r)=-\frac{Z}{r} $$

17) `localEnergy` 图的两条径向能量密度：

$$ \epsilon\,P(r),\quad t_{loc}(r)\,P(r) $$

18) 累计贡献与导数（一般形式）：

$$ \langle U\rangle(<R)=\int_0^R P(r)U(r)dr,\quad \frac{d}{dR}\langle U\rangle(<R)=P(R)U(R) $$

19) 自旋统计缩放：

$$ V_{exch}^{corrected}=\chi V_{exch}^{geo},\quad \chi = \frac{N_{\uparrow}^2 + N_{\downarrow}^2}{N^2} $$

20) 杂化对齐目标函数：

$$ J(q)=\sum \log(\sigma_k(A(q))) $$

---

## 12. 数据/公式来源清单（统一放在文末）

### 12.1 仓库内数据来源

- 设计原则与“信任工程”口径：见 [publication_prep/paper_cn_v3.md](publication_prep/paper_cn_v3.md)。
- Koga 系列 STO 数据表：见 [slater_basis.js](slater_basis.js) 文件头注释（Koga 1999/2000 数据集，以及 `Au_R` 相关数据）。
- 电子-电子有效势查表（H–Kr）：见 [vee_cache.js](vee_cache.js) 与生成脚本 [_experiments/generate_production_cache.py](_experiments/generate_production_cache.py)。
- 交换项角系数表：见 [_experiments/generate_production_cache.py](_experiments/generate_production_cache.py) 的 `EXCHANGE_COEFFS`。

### 12.2 公式来源（标准教材/经典算法，供口径核对）

以下列为“定义与算法的常见出处”，用于说明采用的是标准定义：

- 原子单位制、氢原子解析解、球谐函数与关联勒让德函数：常见量子力学教材（例如 Griffiths《量子力学导论》、Bransden & Joachain《Physics of Atoms and Molecules》等）。
- Slater 型轨道（STO）与归一化：量子化学教材（例如 Szabo & Ostlund《Modern Quantum Chemistry》）及 STO 基组文献常用定义。
- Slater–Condon/Racah 角动量代数、Gaunt 系数、Wigner $3j$ 等：角动量代数/原子结构理论参考书（例如 Edmonds《Angular Momentum in Quantum Mechanics》）。
- Marching Cubes：Lorensen & Cline (1987)。
