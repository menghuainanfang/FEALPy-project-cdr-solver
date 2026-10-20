# 高维鲁棒的能量稳定变分移动网格算法研究及其应用

[科研工作库](../README.md) · [研究网页](index.html) · [论文与学习记录](论文精读/README.md)

本科毕业论文研究仓库，以 FEALPy 为主要计算平台。当前研究底座已经覆盖二维有向几何、一维等分布、固定网格 Poisson 制造解和现有 GFMMPDE 最小基线。

正式题目保持不变；“二维做深、三维验证”和“固定度量下网格能量稳定优先”仍是待导师确认的工作假设。

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
- 项目规划依据：`docs/高维鲁棒能量稳定变分移动网格_研究主档.md`。
