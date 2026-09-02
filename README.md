# 质量扩散对流（对流-扩散-反应）方程求解器

> A finite element solver for the convection–diffusion–reaction (CDR) equation with variable polynomial coefficients, built on [FEALPy](https://github.com/weihuayi/fealpy). Correctness is verified via the method of manufactured solutions (MMS) and convergence-order analysis.
>
> 基于 FEALPy 的变系数对流-扩散-反应方程有限元求解器，含制造解正确性验证、收敛阶分析与可视化展示。

## 项目简介

本项目求解一类「质量扩散对流」方程——即对流-扩散-反应方程（Convection–Diffusion–Reaction, CDR）。一维形式为

$$-(a(x)u')' + b(x)u' + c(x)u = f(x), \qquad x\in(0,1), \qquad u(0)=u_0,\ u(1)=u_1$$

二维推广：

$$-\nabla\cdot\big(a(\boldsymbol x)\nabla u\big) + \boldsymbol b(\boldsymbol x)\cdot\nabla u + c(\boldsymbol x)u = f(\boldsymbol x), \qquad \boldsymbol x\in\Omega\subset\mathbb{R}^2$$

方程左边三项的系数不再是常数 1，而是任意给定的**多项式函数**：

| 项 | 形式 | 物理含义 | 系数条件 |
|---|---|---|---|
| 扩散项 | $-(a(x)u')'$ | 扩散（热量、污染物向四周散开） | $a(x)>0$（椭圆性） |
| 对流项 | $b(x)u'$ | 对流（流场把物质往下游搬运） | 任意多项式 |
| 反应/质量项 | $c(x)u$ | 反应或质量变化（衰减、吸收） | $c(x)\ge 0$（适定性） |

取 $a=1,\ b=0,\ c=0$ 时方程退化回 Poisson 方程 $-\Delta u=f$——因此本项目是 Poisson 方程有限元求解器的自然推广，可视为「一般二阶椭圆型方程」的教学型通用求解器。

## 项目亮点

- **变系数通用求解**：三个系数均可取任意多项式，扩散/对流/反应三项可自由组合（含退化回 Poisson 的对照算例）
- **完整的数值正确性验证管线**：制造解（MMS）反推右端项 → 求解 → 误差表 → **收敛阶与理论阶比对**（P1 元期望 L2 阶 2 / H1 阶 1，P2 元期望 3 / 2）
- **符号推导自动化**：右端项 $f$ 由 SymPy 自动推导，杜绝手算错误
- **一维 + 二维双支持**：一维完成快速验证与手算对照，二维做展示算例
- **全流程工程化**：任务分析 → 调研报告 → 设计方案 → 实现 → 验证 → 可视化，文档齐备

## 技术栈

| 组件 | 用途 |
|---|---|
| Python 3.13 | 开发语言 |
| [FEALPy](https://github.com/weihuayi/fealpy) v3.4 | 有限元框架（网格、函数空间、积分器、边界条件） |
| SciPy | 稀疏线性方程组求解（`spsolve`） |
| SymPy | 制造解右端项与边界值的符号推导 |
| Matplotlib | 解、误差与收敛阶可视化 |

## 数值方法

Galerkin 有限元离散：方程两边乘检验函数 $v\in H_0^1$ 并分部积分，得弱形式

$$\int_\Omega\Big(a\,\nabla u\cdot\nabla v + (\boldsymbol b\cdot\nabla u)\,v + c\,uv\Big)\,dx = \int_\Omega f\,v\,dx$$

离散后总矩阵为三块之和 $\mathbf A = \mathbf K + \mathbf B + \mathbf M$：

- 扩散项 → 对称刚度矩阵 $\mathbf K$；
- 对流项 → **非对称**矩阵 $\mathbf B$（与 Poisson 的关键差异）；
- 反应项 → 对称质量矩阵 $\mathbf M$。

采用 Lagrange 元（P1/P2），系数在积分点局部求值组装（`process_coef_func`），Dirichlet 边界条件按行修正，稀疏直接法求解。

**验证方法学**：给定精确解 $u_{ex}$（取非多项式，如三角函数），由方程反算 $f=-(au_{ex}')'+bu_{ex}'+cu_{ex}$ 与边界值，再解数值方程比对；网格逐级加密（$h\to h/2$）计算收敛阶 $\log_2(e_h/e_{h/2})$，实测阶逼近理论阶即判定实现正确。

## 项目结构（目标结构，代码开发中）

```
.
├── README.md                # 项目说明（本文件）
├── docs/                    # 文档
│   ├── 任务了解.md           # 任务解读：方程背景、交付物、待确认问题
│   ├── 调研报告.md           # 同类求解器调研 + FEALPy 可行性核查 + 风险清单
│   └── 设计方案.md           # 弱形式推导、离散、组装、验证方案（待写）
├── src/                     # 求解器实现（待写）
│   └── cdr_solver.py        # CDR 求解器类：输入 a,b,c,f + 网格 + 元次数，输出解与误差
├── examples/                # 数值算例（待写）
│   ├── example1_polynomial_coefs.py
│   ├── example2_2d.py
│   └── example0_poisson_degenerate.py   # 退化 Poisson 对照
├── verification/            # 正确性验证（待写）
│   ├── manufactured_solution.py         # SymPy 推导 f 与边界值
│   └── convergence.py       # 收敛阶计算与 log-log 图
└── figures/                 # 解、误差、收敛阶图片输出
```

## 快速开始

> 代码开发中，完成后本节将提供可直接复现的命令。

预计使用方式（示意）：

```python
from cdr_solver import CdrSolver

# 系数：a(x)=1+x^2, b(x)=x, c(x)=1+x
solver = CdrSolver(a=lambda x: 1 + x**2, b=lambda x: x,
                   c=lambda x: 1 + x, p=1)
uh = solver.solve(f, gd, nx=40)          # FEALPy 组装 + DirichletBC + spsolve
l2, h1 = solver.compute_errors(u_exact)  # mesh.error
```

## 开发进度

- [x] 任务了解文档：方程背景、交付物清单、技术路线
- [x] 调研报告：同类求解器横向对比（FiPy / deal.II / FreeFEM / OpenFOAM / scikit-fem 等）+ FEALPy 源码级可行性核查
- [ ] 设计方案文档
- [ ] 一维求解器主体 + 制造解验证管线
- [ ] 收敛阶验证（P1/P2，L2/H1）
- [ ] 二维推广算例
- [ ] 可视化与文档整合

## 项目背景

本项目源于算海暑期培训（2026 夏）的课程任务，是前期 Poisson 方程有限元工作（[用 FEALPy 实现的 Poisson 方程算例](https://github.com/weihuayi/fealpy)）的推广：将常系数 Poisson 方程推广为变系数对流-扩散-反应方程的一般求解器，并按要求交付程序、设计方案、数值算例、数值正确性验证与可视化。

## 参考资料

- [FEALPy: Finite Element Analysis Library in Python](https://github.com/weihuayi/fealpy)
- deal.II 教程 step-6 / step-26（变系数组装与时间步进的标准做法）
- FiPy（扩散/对流/反应项自由组合的接口设计参考）
- Roache, *Code Verification by the Method of Manufactured Solutions*, ASME J. Fluids Eng. 2002（MMS 方法论）
- 详细调研见 `docs/调研报告.md`

## 作者

祁靖（作者简介与联系方式）
