# 补充材料 (Supporting Information)

## 1. Clementi-Roetti (1974) 基组解析逻辑

本虚拟实验室使用的核心数据来源于 Clementi 和 Roetti (1974) 的原始文献。数据以 JSON 格式存储，解析逻辑如下：

```javascript
// 简单的 STO 径向波函数计算示例 (JavaScript)
function slaterRadialR(nStar, zeta, r) {
    // 归一化常数 N = (2*zeta)^(nStar + 0.5) / sqrt(factorial(2 * nStar))
    // 注意：nStar 在原始文献中通常为整数
    const normalization = Math.pow(2 * zeta, nStar + 0.5) / 
                         Math.sqrt(factorial(2 * nStar));
    return normalization * Math.pow(r, nStar - 1) * Math.exp(-zeta * r);
}

// 多指数展开 (Linear Combination)
function multiZetaRadialR(basisList, r) {
    let R = 0;
    for (const term of basisList) {
        R += term.coeff * slaterRadialR(term.nStar, term.zeta, r);
    }
    return R;
}
```

## 2. 几何优化杂化算法伪代码

### 算法 1：Thomson 问题梯度下降 (求方向)
1. **输入**：轨道数量 $N$
2. **初始化**：在球面上随机（或按黄金比例分布）生成 $N$ 个单位向量 $v_1, \dots, v_N$。
3. **循环**（直到收敛）：
    a. 计算库仑势能梯度 $G_i = \sum_{j \neq i} \frac{v_i - v_j}{|v_i - v_j|^3}$。
    b. 更新位置 $v_i \leftarrow v_i - \eta \cdot G_i$（$\eta$ 为学习率）。
    c. 投影回球面 $v_i \leftarrow v_i / |v_i|$。
4. **输出**：最优空间指向向量集 $V_{target} = \{v_1, \dots, v_N\}$。

### 算法 2：SVD 求解正交 Procrustes 问题 (求系数)
1. **输入**：目标方向集 $V_{target}$，基组原子轨道集 $\Psi_{basis}$。
2. **构建矩阵**：$A_{ij} = \Psi_j(v_i)$。
3. **SVD 分解**：$A = U S V^T$。
4. **正交转换矩阵**：$Q = U V^T$。
5. **输出**：系数矩阵 $Q_{ij}$，即为第 $i$ 个杂化轨道中第 $j$ 个基组轨道的系数。

## 3. 多平台适配与性能 (BYOD)

为支持“自带设备 (BYOD)”进入课堂，本工具进行了深度性能优化：
- **Web Workers**：采样计算完全卸载至后台线程，确保 iPad/手机浏览器不卡死。
- **浮点纹理/TypedArrays**：使用 JavaScript 高性能数组处理百万级采样点。
- **响应式设计**：面板采用玻璃态弹性布局，自动适配移动端触控。

## 4. 教学指南：15 分钟课堂活动 (教案示例)

**主题**：探索屏蔽效应与有效核电荷
**目标**：通过观察对比，理解为何同一周期的轨道尺寸显著减小。

| 时间 | 活动步骤 | 教师引导语 |
|---|---|---|
| 0-5 min | 选择 H-2p 观察 | “大家看，氢原子的激发态 2p 轨道，它的『云』分布在多大范围内？” |
| 5-10 min | 使用“比照模式”添加 C-2p 和 O-2p | “现在添加碳和氧。注意它们的颜色。为什么碳的 2p 比氢小，氧的又比碳小？” |
| 10-15 min | 观察势能图 (dE/dr) | “看下方的能谱图。我们可以看到，随着核电荷 Z 增加，势能坑变得更深，电子被更紧地束缚在核心附近。” |

**结论总结**：电子间的屏蔽虽然存在，但核电荷的增加（Z）是决定轨道收缩的主导因素。
