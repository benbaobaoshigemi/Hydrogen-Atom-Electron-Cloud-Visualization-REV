"""
计算 [px, py, pz, dz²] 被当作 sp³ 处理时的结果
"""
import numpy as np

# sp³ 系数矩阵（设计用于 [s, px, py, pz]）
c4 = 0.5
sp3_coeffs = np.array([
    [c4, c4, c4, c4],     # ψ₁
    [c4, c4, -c4, -c4],   # ψ₂
    [c4, -c4, c4, -c4],   # ψ₃
    [c4, -c4, -c4, c4]    # ψ₄
])

print("=" * 60)
print("用户选择的轨道（忘记s）：[px, py, pz, dz²]")
print("代码排序后：[px, py, pz, dz²]（按 l 升序）")
print("代码认为这是 4 个轨道 → 套用 sp³ 系数矩阵")
print("=" * 60)

orbital_names = ["px", "py", "pz", "dz²"]

print("\n生成的4个'杂化轨道'：")
print("-" * 60)

for i in range(4):
    terms = []
    for j in range(4):
        coeff = sp3_coeffs[i, j]
        sign = "+" if coeff > 0 else "-"
        if j == 0:
            terms.append(f"{abs(coeff):.1f}·{orbital_names[j]}")
        else:
            terms.append(f" {sign} {abs(coeff):.1f}·{orbital_names[j]}")
    print(f"ψ_{i+1} = {''.join(terms)}")

print("\n" + "=" * 60)
print("关键分析：这些轨道为什么是'两个腰果交叉'的诡异形状？")
print("=" * 60)

print("""
【物理解释】

在正常 sp³ 杂化中（[s, px, py, pz]）：
  - s 轨道是球对称的，提供"各向同性基座"
  - [+,+,+,+] 系数使 s 的球对称与 px+py+pz 的(1,1,1)方向叠加
  - 结果：4个等价瓣指向四面体顶点

但在您的 [px, py, pz, dz²] 组合中：
  - 没有球对称的 s 轨道！
  - dz² 是"沿 z 轴的甜甜圈+两极"结构，本身就是不对称的
  
当 sp³ 系数作用于这个组合时：

ψ₁ = 0.5·px + 0.5·py + 0.5·pz + 0.5·dz²
     ↑________________↑         ↑
     三个p叠加指向(1,1,1)     dz²在z轴正向增强
     
  → 结果：一个指向(1,1,1)方向但被dz²拉长的不对称"腰果"

ψ₂ = 0.5·px + 0.5·py - 0.5·pz - 0.5·dz²
     ↑________________↑         ↑
     p部分指向(1,1,-1)      -dz²翻转z轴贡献
     
  → 结果：另一个指向(1,1,-1)但形状不同的变形瓣

ψ₃ 和 ψ₄ 类似，但由于 dz² 的 z 轴极性，
它们与 ψ₁、ψ₂ 不等价！
""")

print("\n" + "=" * 60)
print("核心发现")
print("=" * 60)
print("""
【为什么 sp³d 看起来对称，而 p³d 看起来诡异？】

sp³d（5轨道）使用的系数矩阵（第421-437行）：
  赤道面3个：只用 s, px, py（dz² 系数为 0）
  轴向2个：  只用 pz, dz²（s, px, py 系数为 0）
  
这是一个精心设计的分层结构，确保对称性！

但 p³d（4轨道误当sp³处理）使用 sp³ 矩阵：
  所有4个杂化轨道都同时混合 px, py, pz, dz²
  没有任何分层，dz² 的 z 轴极性破坏了四面体对称性

→ 结果：4个完全不等价的"腰果"形状
""")
