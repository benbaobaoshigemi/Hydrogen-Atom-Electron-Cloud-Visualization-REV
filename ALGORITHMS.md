# 氢原子电子云可视化 - 算法工程文档

本文档详细说明项目中每个核心算法的设计原理、数学基础和工程实现。

---

## 目录

1. [量子力学计算引擎 (physics.js)](#1-量子力学计算引擎)
2. [采样算法系统 (sampling.js)](#2-采样算法系统)
3. [等值面生成算法 (marching_cubes.js)](#3-等值面生成算法)
4. [可视化渲染系统 (visualization.js)](#4-可视化渲染系统)

---

## 1. 量子力学计算引擎

**文件**: `physics.js` (1840行)

### 1.1 数学基础函数

#### 阶乘与二项式系数
```javascript
// 预计算阶乘表 (0! 到 64!)，避免运行时递归
const FACT = [1, 1, 2, 6, 24, ...]; // 预计算到 64!
function factorial(n) { return FACT[n] ?? Infinity; }
function binomialInt(n, k) { return n! / (k! × (n-k)!); }
```

#### 广义 Laguerre 多项式

$$L_k^{(\alpha)}(x) = \sum_{i=0}^{k} (-1)^i \binom{k+\alpha}{k-i} \frac{x^i}{i!}$$

用于径向波函数 $R_{nl}(r)$ 的计算。

#### 连带 Legendre 函数

$$P_l^m(x) = (-1)^m (1-x^2)^{m/2} \frac{d^m}{dx^m} P_l(x)$$

实现采用 **Condon-Shortley 相位约定**，确保与化学教科书一致。

---

### 1.2 氢原子波函数

#### 径向波函数 $R_{nl}(r)$

$$R_{nl}(r) = N_{nl} \cdot \rho^l \cdot e^{-\rho/2} \cdot L_{n-l-1}^{2l+1}(\rho)$$

其中：
- $\rho = \frac{2Zr}{na_0}$（无量纲径向坐标）
- $N_{nl}$ 为归一化常数
- 归一化条件：$\int_0^\infty |R_{nl}|^2 r^2 dr = 1$

```javascript
function radialR(n, l, r, Z = 1, a0 = 1) {
    const rho = (2 * Z * r) / (n * a0);
    const k = n - l - 1;
    const pref = Math.pow(2*Z/(n*a0), 1.5) * 
                 Math.sqrt(factorial(n-l-1) / (2*n*factorial(n+l)));
    return pref * Math.exp(-rho/2) * Math.pow(rho, l) * 
           generalizedLaguerre(k, 2*l+1, rho);
}
```

#### 实球谐函数 $Y_l^m(\theta, \phi)$

复球谐函数的实线性组合，符合化学惯例：

| m 值 | 类型 | 定义 |
|------|------|------|
| m = 0 | - | $Y_l^0$ (纯实) |
| m > 0 | cos | $\sqrt{2} \cdot \text{Re}(Y_l^m)$ |
| m > 0 | sin | $\sqrt{2} \cdot \text{Im}(Y_l^m)$ |

例如：$p_x \sim Y_1^1(\cos)$，$p_y \sim Y_1^1(\sin)$，$p_z \sim Y_1^0$

---

### 1.3 杂化轨道理论

#### 杂化系数矩阵

杂化轨道是原子轨道的正交线性组合：$\psi_{hybrid} = \sum_i c_i \psi_i$

支持的杂化类型及系数矩阵：

| 类型 | 轨道数 | 几何构型 |
|------|--------|----------|
| sp   | 2      | 直线型 (180°) |
| sp²  | 3      | 三角平面 (120°) |
| sp³  | 4      | 四面体 (109.5°) |
| sp³d | 5      | 三角双锥 |
| sp³d² | 6     | 八面体 |

```javascript
// sp³ 杂化系数矩阵示例
const sp3Matrix = [
    [0.5,  0.5,  0.5,  0.5],   // 指向 (+,+,+)
    [0.5,  0.5, -0.5, -0.5],   // 指向 (+,+,-,-)
    [0.5, -0.5,  0.5, -0.5],   // 指向 (+,-,+,-)
    [0.5, -0.5, -0.5,  0.5]    // 指向 (+,-,-,+)
];
```

#### 杂化轨道密度计算

```javascript
function hybridDensity3D(paramsList, r, theta, phi) {
    let psi = 0;
    for (const p of paramsList) {
        psi += p.coefficient * radialR(p.n, p.l, r) * 
               realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
    }
    return psi * psi;  // |Ψ|²
}
```

---

### 1.4 主轴搜索算法

为杂化轨道找到最大概率密度方向（用于网格极点对齐）：

**算法流程**：
1. **粗搜索**：36×72 网格遍历整个球面
2. **细化搜索**：在最佳点附近迭代精细搜索
3. **解析优先**：对标准杂化类型使用解析解（如 sp³ 的四面体方向）

```javascript
function findHybridPrincipalAxis(paramsList, hybridIndex) {
    // 粗搜索
    for (theta = 0; theta <= π; step = π/36) {
        for (phi = 0; phi < 2π; step = 2π/72) {
            intensity = getDensityAt(theta, phi);
            if (intensity > maxIntensity) updateBest();
        }
    }
    // 细化搜索 (3轮迭代)
    for (iter = 0; iter < 3; iter++) {
        localSearch(searchRange / 2^iter);
    }
    return cartesianCoordinates(bestTheta, bestPhi);
}
```

---

## 2. 采样算法系统

**文件**: `sampling.js` (1456行)

### 2.1 采样方法概述

| 方法 | 适用场景 | 接受率 | 物理精确性 |
|------|----------|--------|------------|
| 逆 CDF 采样 | 单轨道径向 | 100% | 精确 |
| 重要性采样 | 单轨道完整 | ~60-80% | 精确 |
| 杂化精确采样 | 杂化轨道 | ~40-60% | 精确 |
| 拒绝采样 | 回退方案 | ~5-20% | 精确 |

---

### 2.2 逆 CDF 采样（最高效）

**原理**：对径向概率密度 $P(r) = r^2 |R_{nl}(r)|^2$ 构建累积分布函数：

$$F(r) = \int_0^r P(r') dr'$$

然后使用均匀随机数 $u \in [0,1]$，通过二分查找求 $F^{-1}(u)$。

```javascript
function buildRadialCDF(n, l, numPoints = 2000) {
    const rMax = 4 * n² * a0;
    const dr = rMax / numPoints;
    
    // 梯形法则数值积分
    for (i = 1; i <= numPoints; i++) {
        cumulative += (P_prev + P_curr) * dr / 2;
        cdf[i] = cumulative;
    }
    
    // 归一化到 [0, 1]
    for (i = 0; i <= numPoints; i++) {
        cdf[i] /= totalProb;
    }
    return { r, cdf, rMax };
}

function sampleRadialExact(n, l) {
    const u = Math.random();
    // 二分查找 + 线性插值
    return interpolatedRadius;
}
```

**优势**：100% 接受率，无拒绝采样开销

---

### 2.3 重要性采样

**策略**：分离变量采样
1. **径向**：使用精确逆 CDF 采样
2. **角向**：均匀球面采样 + 权重接受-拒绝

```javascript
function importanceSample(n, l, angKey, samplingBoundary) {
    // 第一步：径向采样（精确）
    const r = sampleRadialExact(n, l);
    
    // 第二步：角向采样
    const { theta, phi } = sampleUniformSphere();
    
    // 第三步：角向权重
    const Y2 = realYlm_abs2(angKey.l, angKey.m, angKey.t, theta, phi);
    const w_angular = 4π × Y2;
    
    // 接受-拒绝
    if (Math.random() × maxAngularWeight > w_angular) {
        return { accepted: false };
    }
    
    // 转换为笛卡尔坐标 + 亚像素抖动
    return { x, y, z, r, theta, phi, accepted: true };
}
```

**亚像素抖动**：添加 ±0.005 的随机偏移，抑制摩尔纹

---

### 2.4 杂化轨道采样

**挑战**：杂化波函数 $\psi = \sum c_i R_i Y_i$ 的密度分布复杂

**解决方案**：混合提议分布 + 权重校正

```javascript
function hybridPreciseSample(paramsList, samplingBoundary) {
    // 构建杂化径向 CDF
    // P_hybrid(r) = r² × Σ |c_i|² × |R_i(r)|²
    const cdfData = buildHybridRadialCDF(paramsList);
    
    // 从杂化 CDF 采样 r
    const r = sampleFromCDF(cdfData);
    
    // 均匀球面采样角度
    const { theta, phi } = sampleUniformSphere();
    
    // 计算杂化波函数
    let psi = 0;
    for (p of paramsList) {
        psi += p.coefficient × R(r) × Y(theta, phi);
    }
    
    // 角向权重校正
    const angularWeight = |psi|² / expectedDensity;
    
    // 接受-拒绝
    if (Math.random() × maxWeight > angularWeight) {
        return { accepted: false };
    }
    
    return { x, y, z, psi, accepted: true };
}
```

---

### 2.5 Web Worker 并行采样

**架构**：Worker 池模式

```
┌─────────────────────────────────────────────────────┐
│                    主线程                            │
│  ┌───────────┐  ┌───────────┐  ┌───────────────┐   │
│  │ UI 交互   │  │ 渲染循环   │  │ 结果合并      │   │
│  └───────────┘  └───────────┘  └───────────────┘   │
└────────────────────────┬────────────────────────────┘
                         │ postMessage
         ┌───────────────┼───────────────┐
         ▼               ▼               ▼
    ┌─────────┐    ┌─────────┐    ┌─────────┐
    │ Worker 1│    │ Worker 2│    │ Worker N│
    └─────────┘    └─────────┘    └─────────┘
```

**配置**：
- 池大小：`navigator.hardwareConcurrency` (CPU 核心数)
- 任务分配：每个 Worker 处理 `totalAttempts / workerCount` 次采样
- 会话管理：`samplingSessionId` 防止旧结果污染

---

### 2.6 滚动生成算法

**场景**：采样完成后持续更新部分点

**算法**：
1. 随机选择现有点进行替换
2. 使用重要性采样生成新点
3. 高亮新生成的点（渐变动画）

```javascript
function performRollingUpdate() {
    const pointsToUpdate = Math.floor(pointCount / 100);  // 1%
    
    for (k = 0; k < pointsToUpdate; k++) {
        const targetIndex = randomInt(pointCount);
        
        // 检查：跳过隐藏轨道的点
        if (isHidden(targetIndex)) continue;
        
        // 从映射表移除旧点
        removeFromOrbitalMap(targetIndex);
        
        // 采样新点
        const newPoint = sampleNewPoint();
        
        // 更新位置和颜色
        updatePosition(targetIndex, newPoint);
        
        // 添加高亮动画
        registerHighlight(targetIndex);
    }
}
```

---

## 3. 等值面生成算法

**文件**: `marching_cubes.js` (411行)

### 3.1 Marching Cubes 原理

**输入**：3D 标量场 $f(x,y,z)$ 和等值 $\sigma$

**输出**：满足 $f = \sigma$ 的三角网格表面

**核心步骤**：
1. 将空间划分为立方体网格
2. 对每个立方体的 8 个顶点采样
3. 根据顶点在等值面内/外的状态，确定交叉边和三角形

### 3.2 查找表

```javascript
// 256 种顶点状态 → 对应的三角形边索引
const TRI_TABLE = [
    [],                          // 0: 全在外
    [0, 8, 3],                   // 1: 顶点0在内
    [0, 1, 9],                   // 2: 顶点1在内
    [1, 8, 3, 9, 8, 1],          // 3: 顶点0,1在内
    // ... 256 种情况
];

// 12 条边的顶点连接关系
const EDGE_TABLE = [
    [0, 1], [1, 2], [2, 3], [3, 0],  // 底面
    [4, 5], [5, 6], [6, 7], [7, 4],  // 顶面
    [0, 4], [1, 5], [2, 6], [3, 7]   // 竖边
];
```

### 3.3 等值面提取

```javascript
function marchingCubes(calcPsi, bounds, resolution, isovalue) {
    // 1. 构建 3D 采样网格
    const grid = new Float32Array((resolution+1)³);
    for (iz, iy, ix) {
        grid[idx] = calcPsi(x, y, z);
    }
    
    // 2. 遍历每个立方体
    for (iz, iy, ix) {
        const values = getCubeVertexValues();
        
        // 3. 提取正相位等值面 (ψ = +threshold)
        processCube(values, +isovalue, positiveTriangles);
        
        // 4. 提取负相位等值面 (ψ = -threshold)
        processCube(values, -isovalue, negativeTriangles);
    }
    
    return { positive: positiveTriangles, negative: negativeTriangles };
}
```

### 3.4 连通域分离

**目的**：将单个等值面分离为独立的波瓣

**算法**：并查集 (Union-Find)

```javascript
function separateComponents(triangles) {
    const n = triangles.length / 3;  // 三角形数量
    const parent = new Int32Array(n);
    
    // 初始化：每个三角形是独立的
    for (i = 0; i < n; i++) parent[i] = i;
    
    // 构建顶点 → 三角形映射
    const vertexToTris = new Map();
    for (tri of triangles) {
        for (vertex of tri) {
            const key = quantize(vertex);  // 空间量化
            vertexToTris.get(key).push(tri.index);
        }
    }
    
    // 合并共享顶点的三角形
    for (tris of vertexToTris.values()) {
        for (i = 1; i < tris.length; i++) {
            union(tris[0], tris[i]);
        }
    }
    
    // 按根节点分组
    return groupByRoot(triangles);
}
```

---

## 4. 可视化渲染系统

**文件**: `visualization.js` (1337行)

### 4.1 3D 角向分布网格

**目的**：显示球谐函数 $|Y_l^m(\theta, \phi)|^2$ 的形状

**实现**：基于 Icosahedron 细分球

```javascript
function createOrbitalMesh(group, params, baseRadius) {
    // 1. 创建高密度 Icosphere (detail=50 → ~156k 顶点)
    const geometry = new THREE.IcosahedronGeometry(baseRadius, 50);
    
    // 2. 计算对称轴，确定网格极点方向
    const symmetryAxis = getOrbitalSymmetryAxis(angKey);
    const quaternion = new THREE.Quaternion()
        .setFromUnitVectors(spherePole, symmetryAxis);
    
    // 3. 遍历顶点，计算角向强度并变形
    for (vertex of vertices) {
        // 转换到物理坐标系
        const physDir = vertex.clone().applyQuaternion(quaternion);
        const { theta, phi } = toSpherical(physDir);
        
        // 计算角向强度
        const intensity = calcAngularIntensity(theta, phi);
        
        // 径向变形
        const newRadius = baseRadius × (1 + √intensity × 1.5);
        vertex.multiplyScalar(newRadius / baseRadius);
    }
    
    // 4. 创建线框材质
    const material = new THREE.MeshBasicMaterial({
        wireframe: true,
        blending: THREE.AdditiveBlending,
        transparent: true
    });
}
```

### 4.2 95% 等值面轮廓

**两种实现方案**：

#### A. 点云高亮方案（性能优先）

```javascript
function enableContourHighlight() {
    // 1. 计算每个点的密度排名
    for (point of points) {
        density = calculateDensity(point);
        densityRanks[point.index] = density;
    }
    
    // 2. 按密度排序
    sortedIndices.sort((a, b) => density[b] - density[a]);
    
    // 3. 计算排名百分位
    for (rank = 0; rank < pointCount; rank++) {
        contourRanks[sortedIndices[rank]] = rank / pointCount;
    }
    
    // 4. 高亮 95% 附近的点
    for (point of points) {
        if (|contourRanks[point] - 0.95| < 0.01) {
            point.color = highlightColor;
        } else {
            point.color = dimmedColor;
        }
    }
}
```

#### B. Marching Cubes 方案（质量优先）

```javascript
function createContourMesh(group, baseRadius) {
    // 1. 计算 5% 分位的 |ψ| 值作为阈值
    const psiValues = points.map(p => |calcPsi(p)|).sort();
    const isovalue = psiValues[Math.floor(length × 0.05)];
    
    // 2. 运行 Marching Cubes
    const result = MarchingCubes.run(calcPsi, bounds, resolution=120, isovalue);
    
    // 3. 分离连通域，创建独立网格
    for (component of result.positive) {
        const mesh = createMesh(component, positiveColor);
        group.add(mesh);
    }
    for (component of result.negative) {
        const mesh = createMesh(component, negativeColor);
        group.add(mesh);
    }
}
```

### 4.3 相位着色

**颜色约定**（与采样模块统一）：

| 波函数符号 | 颜色 | RGB |
|------------|------|-----|
| ψ > 0 (正) | 蓝色 | (0, 0.5, 1) |
| ψ < 0 (负) | 红色 | (1, 0.2, 0.2) |
| ψ = 0 (中) | 白色 | (1, 1, 1) |

```javascript
function updatePointColors() {
    const phaseColors = constants.phaseColors;
    
    for (point of points) {
        const psi = calculatePsi(point);
        
        if (psi > 0) {
            color = phaseColors.positive;  // 蓝
        } else if (psi < 0) {
            color = phaseColors.negative;  // 红
        } else {
            color = phaseColors.neutral;   // 白
        }
    }
}
```

### 4.4 黄金螺旋等值面点生成

用于创建均匀分布的等值面高亮点：

```javascript
function createContourInterpolationPoints(contourPoints) {
    const numPoints = 30000;
    const goldenAngle = π × (3 - √5);  // ~137.5°
    
    for (i = 0; i < numPoints; i++) {
        // 黄金螺旋均匀分布
        const y = 1 - (i / (numPoints - 1)) × 2;  // [-1, 1]
        const theta = arccos(y);
        const phi = goldenAngle × i;
        
        // 二分搜索找等值面半径
        let rMin = 0.1, rMax = avgRadius × 2.5;
        for (iter = 0; iter < 20; iter++) {
            const rMid = (rMin + rMax) / 2;
            if (calcDensity(rMid, theta, phi) > threshold) {
                rMin = rMid;
            } else {
                rMax = rMid;
            }
        }
        const radius95 = (rMin + rMax) / 2;
        
        // 生成点
        addPoint(toCartesian(radius95, theta, phi));
    }
}
```

---

## 附录：性能优化

### 缓存策略

| 缓存对象 | 用途 | 失效条件 |
|----------|------|----------|
| `_cdfCache` | 径向 CDF 表 | 轨道变更 |
| `_peakCache` | 径向峰位置 | 轨道变更 |
| `_maxWeightCache` | 最大权重估计 | 轨道变更 |
| `coeffMatrix` | 杂化系数矩阵 | 轨道数变更 |

### 内存管理

- 使用 `Float32Array` / `Int16Array` 减少内存占用
- 点云几何体复用（修改顶点而非重建）
- 及时释放 Three.js 资源（geometry.dispose, material.dispose）

### 并行化

- Web Worker 池：CPU 核心数个 Worker
- 任务粒度：每个 Worker 处理 ~10k 次采样尝试
- 会话 ID 机制防止竞态条件

---

*文档版本: 2024-12-06*
