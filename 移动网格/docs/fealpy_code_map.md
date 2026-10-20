# FEALPy 移动网格源码地图

基线：FEALPy 3.4.0，提交 `e033533bd8e24f796f9c5284685b20a98b4d11cb`，本机源码目录 `C:\Users\10251\fealpy`。

## 公共入口与调度

| 实际文件 | 类/函数 | 输入、输出与作用 |
| --- | --- | --- |
| `fealpy/mmesh/mmesher.py` | `MMesher` | 接收 mesh、离散解 `uh`、有限元空间、`beta` 与 `Config`；按 `active_method` 动态加载实现。`initialize()` 建立实例，`run()` 返回更新后的 mesh 与 `uh`。 |
| `fealpy/mmesh/config.py` | `Config` | 保存方法、监控函数、平滑、插值、伪时间、容差与最大迭代数。默认方法是 `Harmap`；本轮显式选择 `GFMMPDE`。 |
| `fealpy/mmesh/base.py` | `PREProcessor` | 读取 `mesh.meshdata['vertices']`，建立逻辑网格、边界节点次序、P1/参数空间、单元与自由度映射及几何尺度。 |

支持的主网格类型在 `mmesher.py` 的 `_U` 中列出：TriangleMesh、TetrahedronMesh、QuadrangleMesh、HexahedronMesh 及二维高阶三角形/四边形。GFMMPDE 工厂在 `fealpy/mmesh/gfmmpde.py` 中按拓扑维数选择 `GFMMPDE2d`、`GFMMPDE3d` 或二维高阶实现。

## 本轮 GFMMPDE 调用链

| 环节 | 代码位置 | 已确认行为 |
| --- | --- | --- |
| 监控量 | `fealpy/mmesh/monitor.py` | `Monitor` 依据配置选择监控函数；本轮使用默认 `arc_length`。监控量由当前离散解导出，不是固定的网格能量。 |
| 平滑/正则化 | `fealpy/mmesh/monitor.py` 与 `gfmmpde_2d.py` | 配置选择 `heatequ` 平滑；每轮还以 0.2 权重混合上一轮监控量。 |
| 网格方程装配 | `fealpy/mmesh/gfmmpde_2d.py` | `fast_matrix_assembly()` 在二维单纯形上装配扩散、对流、质量项；`vector_assembly()` 组装右端；`func_solver()` 施加 Dirichlet 边界并解分块系统。 |
| 边界移动 | `gfmmpde_2d.py:_process_boundary_group` | 顶点固定在多边形顶点，边界节点沿各边求一维重新分布；需要 `mesh.meshdata['vertices']`。 |
| 几何可行步 | `gfmmpde_2d.py:_get_physical_node` | 由候选位移生成局部二次约束系数，选正根并乘 `alpha` 缩放更新。本轮仍独立检查所有单元的有向面积。 |
| 伪时间/停止 | `gfmmpde_2d.py:mesh_redistributor` | 每轮监控量→平滑→线性系统→可行更新→解转移；以最大节点位移与 `tol` 比较，或达到 `maxit`。默认容差由逻辑网格尺度计算。 |
| 解转移 | `fealpy/mmesh/interpolater.py` | `Interpolater` 提供多种转移策略；GFMMPDE 每轮调用当前 `interpolate(node)` 更新离散解。本轮使用默认 `comass` 路径。 |
| 质量诊断 | `fealpy/mmesh/mesh_quality.py` | `MeshQuality` 可计算等分布、对齐和几何质量指标；这些指标不能自动命名为离散网格能量。 |

映射关系是“物理网格上的监控量驱动逻辑/物理网格重分布”。`gfmmpde_2d.py` 同时持有 `self.node` 与 `self.logic_mesh.node`，矩阵装配显式使用逻辑节点。后续推导必须沿当前实现确认未知量与映射方向，不能直接套用 Huang–Kamenski 2015 的另一套几何离散公式。

## 能量与非退化能力边界

- `GFMMPDE2d.mesh_redistributor()` 不返回离散网格能量历史，本轮没有找到可直接比较的 GFMMPDE 能量接口。
- `fealpy/mmesh/metrictensoradaptive.py`、`eagadaptivehuang.py` 和 `eagadaptivefb.py` 的 `mesh_redistributor(return_info=True)` 会记录 `I_h`、相邻步能量差与最小单元测度。这些方法应作为下一轮能量诊断候选，但尚未在本轮运行。
- GFMMPDE 的 `_get_physical_node` 含局部可行步长构造；源码并未在公共返回值中给出整条节点运动路径的非退化证明或全局无重叠诊断。
- FEALPy `entity_measure('cell')` 用于通常的单元测度，不能代替方向检查。本项目在 `src/geometry.py` 单独保留符号面积。
- 二维和三维 GFMMPDE 类均存在；“存在三维类”只说明代码路径，不构成三维数值验证。

## PDE 装配路径

`experiments/c_poisson_manufactured.py` 使用当前 API：TriangleMesh → LagrangeFESpace(P1) → BilinearForm/ScalarDiffusionIntegrator → LinearForm/ScalarSourceIntegrator → DirichletBC → FEALPy 稀疏矩阵 `to_scipy()` → SciPy `spsolve`。误差由 `mesh.error()` 对精确解/梯度与有限元函数/梯度进行积分。

## 现有能力、缺口和个人工作边界

FEALPy 已提供二维/三维移动网格框架、多个监控/平滑/插值策略、边界处理、PDE 示例和若干质量指标。本项目当前新增的是独立可运行的几何方向测试、等分布教学实验、Poisson 收敛基线和带有向面积诊断的最小 GFMMPDE 运行记录；它们是研究基础设施，不是新算法。

下一项真实缺口是：在一个明确暴露 `I_h` 的现有方法上，固定度量与能量定义，核对解析/实现梯度、全离散能量变化、正定向和步长策略之间的关系。创新性仍未知。

