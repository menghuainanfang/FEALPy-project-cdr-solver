# 本轮交接

日期：2026-09-27

## 已完成

- 建立项目专用 venv，并从本机 FEALPy 源码 editable 安装 3.4.0。
- 建立 README、环境审计、源码地图、研究计划、实验协议和本交接页。
- 完成 A：FEALPy 二维三角网格与有向几何核验。
- 完成 B：主档 E0 一维等分布教学实验。
- 完成 C：固定网格 P1 Poisson 制造解四级收敛实验。
- 完成 D：现有 GFMMPDE 最小二维基线，记录正定向、面积、位移、插值误差和耗时。
- 几何单元测试 3 项全部通过。

## 关键结果

| 核验 | 状态 | 结果 |
| --- | --- | --- |
| A | 通过 | 2×2 方形网格 8 个单元，最小有向面积 0.125，总面积 1；手算三角形面积 1，翻转 -1，退化 0 |
| B | 通过 | 40 段最大监控质量相对偏差 `1.92e-9`；自适应/均匀 L2 插值误差比 0.9522 |
| C | 通过 | n=8→64，L2 阶 1.9745/1.9935/1.9984，H1 阶 0.9891/0.9973/0.9993；边界误差 0；最大相对残差 `1.37e-13` |
| D | 通过（最小代码路径） | 64 单元，移动后最小有向面积 0.009569，总面积 1，最大位移 0.09748；L2 插值误差 0.3296→0.3371；未提供离散能量历史 |

原始结果：`results/a_geometry.json`、`results/b_equidistribution.json`、`results/b_equidistribution_segments.csv`、`results/c_poisson_convergence.json`、`results/c_poisson_convergence.csv`、`results/d_mmesher_baseline.json`。节点图在 `results/b_equidistribution.png`。

## 复现命令

```powershell
cd D:\桌面\moving-mesh-thesis
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m experiments.a_geometry
.\.venv\Scripts\python.exe -m experiments.b_equidistribution_1d
.\.venv\Scripts\python.exe -m experiments.c_poisson_manufactured
.\.venv\Scripts\python.exe -m experiments.d_mmesher_baseline
```

## 已知问题

- FEALPy 源码仓库整体为 dirty，虽未见移动网格核心文件改动，正式基准仍应在导师指定分支或干净检出上复核。
- D 是缩小后的 4×4、3 步 GFMMPDE 代码路径核验；插值误差略增，不能作为适应性收益证据。
- GFMMPDE 当前接口没有离散网格能量历史；不能据此声称能量稳定。
- 尚未运行三维、Burgers、MetricTensorAdaptive/EAGAdaptive、方向差分或强各向异性扫描。

## 接下来三件事

1. 读清 `metrictensoradaptive.py` 中 `I_h` 的公式与冻结 M 约定，用单三角形手算和中心方向差分核对梯度。
2. 在二维小网格上运行 `return_info=True`，同时记录 `I_h`、最小有向面积、步长与停止原因，寻找真实失效而非先改算法。
3. 自己从 Poisson 强形式推到弱形式，并解释本轮 L2 二阶、H1 一阶和边界/代数残差分别验证了什么。

