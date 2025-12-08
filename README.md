# 电子云可视化 (Electron Cloud Visualization)

一个基于 WebGL 的原子与杂化轨道电子云可视化工具。本项目在原氢原子模型基础上进行了扩展，集成了 Slater-Type Orbitals (STO) 和量子力学计算，支持模拟从氢到氙（He-Xe）的多电子原子以及杂化轨道系统。

采用蒙特卡洛采样算法生成电子概率密度分布的3D点云，并使用 Marching Cubes 算法生成等值面。

[License](https://img.shields.io/badge/License-MIT-green)

---

## 核心特性

### 物理计算
- **多电子原子支持**：内置 Clementi-Roetti (1974) STO 基组数据，支持 He 到 Xe 等原子的轨道模拟。
- **相位校正**：实现了渐近分析算法，校正 STO 数据与 Condon-Shortley 相位约定之间的差异。
- **杂化轨道**：支持 $sp, sp^2, sp^3, sp^3d, sp^3d^2$ 杂化，计算包含交叉项（Cross-terms）的完整波函数叠加。
- **标准化约定**：遵循物理学界的 Condon-Shortley 相位定义。

### 可视化功能
- **实时采样**：基于概率密度的重要性采样，生成电子云点云。
- **轨道轮廓**：生成 95% 概率密度的等值面网格（Orbital Linear Mesh）。
- **着色模式**：支持按概率密度或波函数相位（红/蓝）着色。
- **角向分布**：独立显示球谐函数 $Y_{lm}$ 的 3D 形状。

### 交互模式
- **比照模式 (Compare Mode)**：并排对比不同轨道（如观察不同原子的轨道大小差异），支持滚动更新。
- **杂化演示**：观察杂化轨道的形成过程。
- **数据图表**：实时绘制径向概率分布 (RDF) 和角向分布图，并提供理论曲线对照。
- **手势控制**：支持通过摄像头手势（握拳旋转、捏合缩放）控制视图。

---

## 快速开始

### 方式一：直接运行（Windows）

1. 双击运行目录下的 **`Hydrogen_Electron_Cloud.exe`**。
2. 程序会自动启动本地服务器并打开默认浏览器。

### 方式二：源码运行

1. 确保已安装 Python 3.8+。
2. 运行启动脚本：
   ```bash
   python start_server.py
   ```
3. 访问 `http://127.0.0.1:8000`。

---

## 操作说明

### 模式介绍
- **Single (单轨)**：显示单个原子轨道。
- **Multi (多选)**：叠加显示多个轨道。
- **Compare (比照)**：分屏对比多个轨道。
- **Hybrid (杂化)**：显示杂化轨道组合。

### 界面功能
- **左侧面板**：选择原子（H, He... Xe）、轨道参数（n, l, m）。
- **右侧面板**：调整可视化参数（点数、大小）和显示选项。
- **底部图表**：显示概率分布统计数据。

---

## 物理背景

### 1. 氢原子
单电子原子的解析解：
$$\psi_{n,l,m}(r,\theta,\phi) = R_{nl}(r) \cdot Y_l^m(\theta,\phi)$$

### 2. 多电子原子 (STO)
使用 Slater-Type Orbitals 近似：
$$R_{nl}(r) = \sum_i c_i \frac{(2\zeta_i)^{n^*_i+0.5}}{\sqrt{(2n^*_i)!}} r^{n^*_i-1} e^{-\zeta_i r}$$

### 3. 杂化轨道
包含干涉项的密度计算：
$$|\Psi_{hybrid}|^2 = \sum_i |c_i \psi_i|^2 + \sum_{i \neq j} c_i c_j \psi_i \psi_j^*$$

---

## 文件结构

- `physics.js`: 核心物理计算（波函数、STO、杂化）
- `sampling-worker.js`: 采样 Worker 线程
- `visualization.js`: Three.js 渲染与网格生成
- `slater_basis.js`: STO 基组数据
- `gesture_v2.js`: 手势识别

---

## License

MIT License © 2024
