# 高维鲁棒的能量稳定变分移动网格算法研究及其应用

[科研工作库](../README.md) · [研究网页](index.html) · [论文与学习记录](论文精读/README.md)

本科毕业论文研究仓库，以 FEALPy 为主要计算平台。当前研究底座已经覆盖二维有向几何、一维等分布、固定网格 Poisson 制造解和现有 GFMMPDE 最小基线。

## 研究定位（2026-10-23 更新）

**正式题目暂保持《高维鲁棒的能量稳定变分移动网格算法研究及其应用》。** 已明确的实际研究主线是：**以导师团队2026年的迹—对数变分移动网格论文为基础，面向三维四面体网格完成数值算法实现、验证与应用。** 基础论文为 Wang、Huang、Wei，[*A Trace–Logarithmic Variational Functional for Equidistribution and Alignment in Moving Mesh Adaptation*](https://arxiv.org/abs/2601.20235)。

- **三维是最终核心任务**：围绕迹—对数泛函、四面体几何离散、单元与全局导数、边界处理、求解器及三维数值算例展开。
- **二维是学习和核验基础**：用于理解原理、复现经典工作、检查公式和调试，不再把“二维做深、三维验证”作为课题定位。
- **网格能量与鲁棒性是研究性质和评价对象**：区分连续、半离散、全离散结论；检查固定度量阶段的能量演化与单元正定向，但不再将“固定度量下网格能量稳定优先”设为独立选题假设。网格能量下降不等于物理 PDE 能量稳定。
- **继承与新增贡献必须分清**：原论文部分理论已按一般维数建立；此前已有三维变分网格研究。最终成果重点在特定新泛函的三维可运行实现、数值验证及经证据支持的改进，不声称首次发明三维变分网格。
- **开放问题**：具体时间离散改进、非线性求解优化、复杂三维应用与进一步理论证明，均须由复现结果和导师反馈逐步确定，不能预设研究结论。

以上研究主线依据已沟通的导师任务范围；不代表导师已经逐项认可下面的实验或创新设想。

## 长期学习与研究路线

采用**基础理解 → 手算推导 → 小型可验证实验 → 经典论文复现 → 三维研究**的路径。时间是参考节奏，不是强制截止日期；进阶依据是可复述的理论、可核验的计算和原始证据。

| 参考阶段 | 核心文献与问题 | 主要学习 / 验收产物 |
| --- | --- | --- |
| 2026.10–12 · 几何基础 | Huang 2006《Mathematical Principles》；联系已读 Huang–Kamenski 2015 | 理解 Jacobian、SVD、度量椭球、等分布与对齐；完成二维手算及基础实验 |
| 2026.12–2027.02 · 变分与 MMPDE | Huang 2001《Variational Mesh Adaptation》；按需选读 Huang–Ren–Russell 1994 | 从 EA 条件理解泛函构造，推导简单梯度流并实现小型算例 |
| 2027.02–04 · 度量与误差 | Huang 2005《Metric Tensors》；Huang–Sun 2003；Huang 2005《Measuring Mesh Qualities》 | 构造 Hessian 度量，比较插值误差、等分布及对齐指标 |
| 2027.04–07 · 数值实现与基准 | Huang 2001《Practical Aspects》；Huang–Kamenski–Russell 2015 | 复现经典二维/三维基准中的可行部分，形成可重复程序、图表和误差/耗时记录 |
| 2027.07–09 · 非退化性 | Huang–Kamenski 2018 | 列明定理假设，区分半离散证明与实际全离散算法；设计困难网格压力测试 |
| 2027.10–12 · 新泛函与三维实现 | Wang–Huang–Wei 2026；FEALPy 对应代码 | 梳理已有三维支持与二维专用实现，核验四面体能量/梯度，推进完整三维路径 |
| 2028.01–04 · 深化与收束 | 三维问题驱动补充阅读 | 完成三维数值应用、鲁棒性与成本比较，并整理毕业论文 |

**当前阅读顺序：** 01 → 03 → 04（结合09、10）→ 05、07 → 16 → 导师2026年论文。已阅读的 Huang–Kamenski 2015 作为几何离散桥梁，按问题回看，不必从头重读。02 和 11 保留为背景与旁支选读，不要求近期通读。

**近期小实验（计划，非已完成记录）：** F01 一维等分布；F02 度量椭圆与特征分解；F03 Jacobian/SVD 与 EA 条件；之后再做简单泛函下降、Hessian 度量、质量指标与经典复现。与现有仓库已实现的实验脚本分开管理，不能仅凭列入计划标记为已复现。

**阶段里程碑：** 2027 年 6 月前可视进度整理一份用于夏令营的可选二维基准或小改进成果；不以此强迫三维工作提前。毕业阶段以三维任务为主线，成果和个人贡献均以实际测试及导师确认为准。

## 进度与证据原则

研究进度区分“已阅读 / 已独立推导 / 已编程核验 / 已复现 / 待验证”。**阅读过或程序能运行不等于算法已经被证明，也不等于创新成果。** 每次实验至少保存论文公式对应关系、代码提交、配置、环境、原始结果、失败情况和验收结论。本地交接可维护独立进度文件，以仓库中的可检验记录为准。

下面保留的实验、误差和环境数字来自仓库已有记录；本次 README 更新**没有重新运行数值实验、修改 FEALPy 或验收三维代码**。

## 变分移动网格可视化

```powershell
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual
```

生成 `results/e_variational_mesh.html`，直接用浏览器打开，离线播放扰动平滑、陡变带固定边界、陡变带切向边界三个场景。可点击顶点跟踪四种坐标角色，查看边界位移、能量、有向面积及对齐/等分布指标。

复用 FEALPy `MetricTensorAdaptive` 的含对数项能量、局部导数、装配、平衡和边界投影，以 SciPy BDF 演化冻结度量的计算坐标。新版默认参数 `--nx 16 --duration 0.5 --frames 61`。误差为采样插值 RMSE；公式及边界见 [实验说明](docs/variational_mesh_visual.md)。

## 移动网格上的 PDE 精度对照

```powershell
.\.venv\Scripts\python.exe -m experiments.f_poisson_mesh_comparison
```

生成 [PDE 交互对比](results/f_poisson_comparison.html)：三个稳态 Poisson 制造解、三档规模、均匀/扰动/平滑/自适应固定边界/自适应切向边界五种物理网格，比较 P1 解误差、质量和成本。默认共 45 次主对照求解，参数 `--levels 8 16 32 --duration 0.5`。

默认 289 自由度带状算例，均匀/固定边界自适应/切向自适应的 L² 误差分别约 0.06413、0.02546、0.01522。度量使用已知解结构，移动成本高于一次小型 PDE 求解。公式、全表与限制见 [实验 F 说明](docs/poisson_mesh_comparison.md)。

## 仓库结构

```text
移动网格/
├── docs/           # 研究主档、环境审查、源码地图、计划与交接
├── experiments/    # 可独立运行的编号实验
├── references/     # 文献索引与公开参考实现说明，不存放无权分发的原文
├── results/        # 实验输出说明；可重复生成的数据和图片默认不提交
├── src/            # 自有的可复用几何、网格与算法模块
├── tests/          # 数学性质和回归测试
├── requirements.txt
└── requirements-lock-windows-py313.txt
```

FEALPy 是外部依赖；本机源码已原位更新到官方 develop 提交 `69faa5f41e3c12633ed714520c3ec512eb1a784b`。自有数值驱动在本仓库维护，没有改写库内数值算法。可生成含当前源码快照的离线分发包：`python tools/package_demos.py`。环境见 `docs/environment_audit.md`。

## 环境准备

在本目录创建 Python 3.13 虚拟环境，安装本机 FEALPy 源码与依赖：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e C:\Users\10251\fealpy
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

FEALPy 的版本与接口要求见 [环境审查](docs/environment_audit.md)。其他机器将 FEALPy 路径替换为自己的源码目录。

## 最短运行

```powershell
cd D:\桌面\movingmesh\移动网格
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual
.\.venv\Scripts\python.exe -m experiments.f_poisson_mesh_comparison
```

运行结果写入 `results/`。环境、源码地图、实验口径和本轮交接分别见 `docs/environment_audit.md`、`docs/fealpy_code_map.md`、`docs/experiment_protocol.md` 和 `docs/session_handoff.md`。

## 当前验证状态

| 项目 | 状态 | 主要证据 |
| --- | --- | --- |
| 有向三角形几何 | 通过 | 翻转为负、退化为零、单位方形总面积为 1 |
| 一维等分布 E0 | 通过 | 40 段监控质量最大相对偏差约 `1.92e-9` |
| P1 Poisson 制造解 | 通过 | L2 约二阶、H1 半范数约一阶 |
| GFMMPDE 最小二维路径 | 通过 | 移动后单元保持正定向、总面积为 1 |

## 当前边界

- FEALPy 基线来自本机 `C:\Users\10251\fealpy` 的 3.4.0 源码；精确提交与 dirty 状态见环境审计。
- GFMMPDE 最小基线没有暴露可比较的离散网格能量历史；当前只报告有向面积、位移、插值误差和耗时。
- 未把示例运行、质量指标或能量曲线解释为新的算法贡献。
- 当前没有三维、Burgers 或全离散能量稳定结论。

## 参考来源

- FEALPy: <https://github.com/weihuayi/fealpy>
- Huang and Kamenski, *A geometric discretization and a simple implementation for variational mesh generation and adaptation*, 2015.
- Huang and Kamenski, *On the mesh nonsingularity of the moving mesh PDE method*, 2018.
- 导师团队2026年论文：Wang、Huang、Wei，[*A Trace–Logarithmic Variational Functional for Equidistribution and Alignment in Moving Mesh Adaptation*](https://arxiv.org/abs/2601.20235)。
- Huang 2006，*Mathematical Principles of Anisotropic Mesh Adaptation*；Huang 2001，*Variational Mesh Adaptation: Isotropy and Equidistribution*。
- Huang 2005，*Metric tensors for anisotropic mesh generation*；Huang、Kamenski、Russell 2015，*A comparative numerical study of meshing functionals for variational mesh adaptation*。
- 历史项目主档：docs/高维鲁棒能量稳定变分移动网格_研究主档.md（旧阶段假设以本页2026-10-23更新为准；本次提交未同步修订该文件）。
