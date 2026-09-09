# 一维变系数对流–扩散–反应（CDR）方程有限元求解器

> A 1D finite element solver for the convection–diffusion–reaction (CDR) equation with variable polynomial coefficients, built on [FEALPy](https://github.com/weihuayi/fealpy). Correctness is verified via the method of manufactured solutions (MMS) and convergence-order analysis.
>
> 基于 FEALPy 的一维变系数对流–扩散–反应方程有限元求解器：程序、设计方案、数值算例、数值正确性验证与可视化，完整可复现。

## 方程与项目简介

本项目求解「质量扩散对流」方程——即对流-扩散-反应方程（Convection–Diffusion–Reaction, CDR）。一维稳态形式为

$$-(a(x)u')' + b(x)u' + c(x)u = f(x), \qquad x\in(0,1), \qquad u(0)=u_0,\ u(1)=u_1$$

| 项 | 形式 | 物理含义 | 系数条件 |
|---|---|---|---|
| 扩散项 | $-(a(x)u')'$ | 扩散（热量、污染物向四周散开） | $a(x)>0$（椭圆性） |
| 对流项 | $b(x)u'$ | 对流（流场把物质往下游搬运） | 任意多项式 |
| 反应项 | $c(x)u$ | 反应或质量变化（衰减、吸收） | $c(x)\ge 0$ |

要点：

- 三个系数 a、b、c 均可取**任意多项式**（含常数），右端源项 **f 是独立的必需输入**，可以是任意函数（含 f=0 的齐次情形）；
- 展开后是 $-(au')'=-au''-a'u'$，扩散系数变化时 $-a'u'$ 项不可遗漏；
- 取 $a=1,\ b=0,\ c=0$ 时方程退化回 Poisson 方程 $-\Delta u=f$——本求解器是常系数 Poisson 有限元求解器的自然推广。

## 数值方法

取零端点检验函数 v，仅对扩散项分部积分得弱形式

$$\int_0^1 a u'v' + \int_0^1 b u'v + \int_0^1 c uv = \int_0^1 f v.$$

离散后总矩阵为三块之和 $\mathbf A=\mathbf K+\mathbf B+\mathbf M_c$、右端 $\mathbf F$ 来自 f：

- 扩散项 → 刚度矩阵 $\mathbf K$（对称）
- 对流项 → 矩阵 $\mathbf B$（**非对称**，与 Poisson 求解器的关键差异）
- 反应项 → 质量矩阵 $\mathbf M_c$（对称）

实现采用连续 Lagrange 元（P1/P2），用 FEALPy 积分器组装：`ScalarDiffusionIntegrator(coef=a)`、`ScalarConvectionIntegrator(coef=b)`、`ScalarMassIntegrator(coef=c)`、`ScalarSourceIntegrator(f)`，系数在积分点求值；端点 Dirichlet 用 `DirichletBC` 修正矩阵与右端（非零边界贡献移到右端）；稀疏线性方程组经 FEALPy `spsolve(..., solver='scipy')` 求解。

## 验证方法学（MMS + 收敛阶）

制造解验证管线：给定精确解 $u_*=1+x+\sin(\pi x)$ → SymPy 反推源项 $f=-(au_*')'+bu_*'+cu_*$ → FEALPy 求解 → `mesh.error` 计算 L2/H1 误差 → 网格逐级加密计算经验收敛阶 $\log_2(e_h/e_{h/2})$ → 与理论阶比对（P1 期望 L2≈2/H1≈1，P2 期望 ≈3/≈2）。

四组系数算例：`poisson`（退化对照）、`variable`、`reverse`、`cubic`；每组 P1/P2 × n=8/16/32/64，共 **32 组全部通过**，且另有手算单元矩阵对照、f 作用验证、两形式转换等价性、强对流限制案例检出等独立检查。最细网格 n=64 实测摘要：

| 算例 | p | L2 误差 | L2 阶 | H1 阶 |
| --- | --- | --- | --- | --- |
| poisson | 1 | 1.555290e-04 | 2.00 | 1.00 |
| poisson | 2 | 4.809369e-07 | 3.00 | 2.00 |
| variable | 1 | 1.239706e-04 | 2.00 | 1.00 |
| variable | 2 | 4.809394e-07 | 3.00 | 2.00 |

完整 32 组数据见[docs/一维求解器_任务报告.md](./docs/一维求解器_任务报告.md)与 `results/convergence.csv`。

## 项目结构

```
.
├── README.md
├── run.cmd                    # 一键复跑：环境检查 + 32 组验证 + 独立检查
├── requirements.txt           # 依赖版本锁定（fealpy/numpy/scipy/sympy/matplotlib）
├── .gitignore
├── docs/                      # 项目文档
│   ├── 一维求解器_入门说明.md   # 入门阅读：方程与符号 → 手算例子 → 代码逐段拆解
│   ├── 一维求解器_任务报告.md   # 设计方案、完整程序、32 组实测表、验证方法与图片
│   ├── 质量扩散对流求解器_任务了解.md  # 早期任务分析（方程背景、交付物、待确认清单）
│   └── 质量扩散对流求解器_调研报告.md  # 同类软件调研 + FEALPy 可行性核查
├── cdr_solver.py              # 示例层：主接口 solve_cdr + 制造解 + run()（32 组验证与绘图）
├── cdr_lfem_solver_1d.py      # 核心类 CdrLFEMSolver1D：linear_system / apply_bc / solve
├── example_source.py          # 最小示例：解 -u''=2，f 显式给出
├── verify_integrators.py      # 单元矩阵与手算值对照 + 边界自由度检查
├── verify_formulation.py      # f 作用、多项式精确解、形式转换、强对流限制案例
├── verify_integration.py      # 后端兼容与可集成性检查（NumPy/PyTorch，需自行装 PyTorch）
├── check_environment.py       # 依赖版本核对（与 requirements.txt 比对）
└── results/                   # 验证输出：convergence.csv、verification.png、JSON 等
```

说明：`.py` 脚本之间通过同目录模块导入互相依赖（如 `verify_*.py` 导入 `cdr_solver`），因此代码保持扁平置于根目录；文档统一归入 `docs/`，输出数据与图片在 `results/`。`docs/一维求解器_任务报告.md` 附录内嵌全部源码全文，与根目录文件一一对应。

## 快速开始

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\run.cmd
```

`run.cmd` 依次执行环境检查、32 组收敛验证（输出 `results/convergence.csv` 与 `results/verification.png`）、单元矩阵对照与独立方程检查。

最小示例（`example_source.py`，解 $-\ u''=2$、零端点）：

```python
from cdr_solver import solve_cdr, expression_function as fun

mesh, space, uh, residual, boundary_error = solve_cdr(
    a=fun(1), b=fun(0, vector=True), c=fun(0),
    f=fun(2), gd=fun(0), n=16, p=2)
print('Maximum solution value:', max(uh[:]))
```

主接口 `solve_cdr(a, b, c, f, gd, n=32, p=1)`：a、b、c、f、gd 可为常数或 Cartesian 可调用函数（SymPy 表达式用 `expression_function` 包装；一维下 b 需返回形状 `(..., 1)`，即 `vector=True`）。返回网格、函数空间、数值解 uh 及残差指标。

可复用核心（适合集成到其他 FEALPy 程序）：

```python
from fealpy.mesh import IntervalMesh
from fealpy.functionspace import LagrangeFESpace
from cdr_lfem_solver_1d import CdrLFEMSolver1D

mesh = IntervalMesh.from_interval_domain([0, 1], nx=32)
space = LagrangeFESpace(mesh, p=1)
model = CdrLFEMSolver1D(space, a, b, c, f, gd)   # a/b/c/f/gd 为常数或可调用
uh = model.solve()
```

## 已知边界与范围

- **已实现并验证**：一维稳态、全 Dirichlet、连续 P1/P2、系数光滑且扩散正定的问题；NumPy CPU 全支持，PyTorch CPU 支持 P1（P2 受当前 IntervalMesh 插值点接口限制，显式拒绝）。
- **尚未实现**：二维、时间项、混合边界、稳定化（如 SUPG）、逐单元守恒重构与严格非负性保证。
- **已知局限**：无稳定化 Galerkin 在强对流下会出现伪振荡。程序对网格 Peclet 数 $\mathrm{Pe}=\frac{|b|h}{2a}>1$ 给出警告；`results/form_verification.json` 记录了一个 a=0.001、b=1、Pe≈15.6 的案例出现非物理负值，`accepted_as_supported=false`——该案例被程序正确检出，不作为支持工况。

## 复现环境（2026-09-07 实测）

| 组件 | 版本 |
| --- | --- |
| Python | 3.13.2 |
| fealpy | 3.4.0 |
| numpy | 2.3.4 |
| scipy | 1.16.3 |
| sympy | 1.14.0 |
| matplotlib | 3.10.7 |

另以本地更新的 FEALPy 开发源码复跑全部检查通过（详见[docs/一维求解器_任务报告.md](./docs/一维求解器_任务报告.md)）。

## 参考资料

- [FEALPy: Finite Element Analysis Library in Python](https://github.com/weihuayi/fealpy)
- Roache, *Code Verification by the Method of Manufactured Solutions*, ASME J. Fluids Eng. 2002（MMS 方法论）
- deal.II 教程 step-6 / step-26（变系数组装与时间步进的标准做法）
- FiPy（扩散/对流/反应项自由组合的接口设计参考）
- 详细调研见[docs/质量扩散对流求解器_调研报告.md](./docs/质量扩散对流求解器_调研报告.md)
