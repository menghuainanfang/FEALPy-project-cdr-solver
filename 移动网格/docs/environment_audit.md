# 环境审查

## 当前本地入口

工作库固定为 `D:/桌面/movingmesh`，主线项目位于 `移动网格/`，现用解释器为 `D:/桌面/movingmesh/移动网格/.venv/Scripts/python.exe`，FEALPy 仍从 `C:/Users/10251/fealpy` 导入。

## 当前演示环境（2026-10-17）

本机 `C:/Users/10251/fealpy` 已在原目录快进到官方 develop 提交 `69faa5f41e3c12633ed714520c3ec512eb1a784b`（2026-10-14）。本地分支名仍为 feature；项目 `.venv` 的 editable 导入直接使用更新后的代码，包元数据仍显示 3.4.0，应以源码提交区分版本。

本地三份 tutorial 教学脚本和 AGENTS.md/CLAUDE.md 已保留，另有 `D:/桌面/fealpy_local_backup_2026-10-17` 文件备份与 Git stash。官方新版移除了旧 tutorial 目录，用户的三份脚本已恢复到该路径。没有发布或推送操作。

当前 E/F 直接使用新版 schema 网格和 FEALPy MetricTensorAdaptive 数值核心；19 项测试通过，E 默认三场景和 F 默认 45 次主对照完成。新版接口的 `mesh.error` 使用函数参数和 `cell_axis=True`。A—D 以下历史基线没有据本次 E/F 验证自动认定全部适配新版。

## 历史环境基线

审查日期：2026-09-27（Asia/Shanghai）

## 结论

项目专用环境已建立并通过依赖检查。实际执行路径为 `D:\桌面\moving-mesh-thesis\.venv\Scripts\python.exe`，FEALPy 从本机源码仓库 editable 导入，而不是从系统 site-packages 导入。

## 主机与工具

| 项目 | 实际值 |
| --- | --- |
| 操作系统 | Windows 11 家庭版中文版，10.0.26200 |
| 终端 | PowerShell 7.6.6 |
| CPU | Intel Core Ultra 7 155H，16 核/22 逻辑处理器 |
| 内存 | 31.47 GiB，总审查时约 15.96 GiB 可用 |
| GPU | NVIDIA RTX 4060 Laptop GPU；Intel Arc Graphics |
| 磁盘 | C: 约 44.57 GiB 可用；D: 约 157.87 GiB 可用 |
| Git | 2.53.0.windows.2 |
| Python | CPython 3.13.2 |
| Conda | 24.11.3，可调用；本项目未使用 |
| uv | 未安装 |

本轮使用 NumPy 后端、CPU、float64。GPU 存在，但 A—D 不需要 CUDA 或第二后端。

## FEALPy 基线

| 项目 | 实际值 |
| --- | --- |
| 版本 | 3.4.0 |
| 源码目录 | `C:\Users\10251\fealpy` |
| 分支 | `feature` |
| 提交 | `e033533bd8e24f796f9c5284685b20a98b4d11cb` |
| 项目环境导入位置 | `C:\Users\10251\fealpy\fealpy\__init__.py` |
| 移动网格入口 | `fealpy.mmesh.mmesher.MMesher` |

源码仓库不是干净状态：`tutorial/poisson.py`、`tutorial/poisson_example.py` 有修改，另有未跟踪的仓库说明和 `tutorial/poisson_convergence.py`。本轮没有拉取、覆盖或修改该仓库。移动网格核心目录未出现在已有改动中；研究项目仍把整个 FEALPy 基线标为 dirty，不能只用提交号概括其完整工作树状态。

系统 Python 原先也装有 FEALPy 3.4.0，导入位置为系统 site-packages；项目 venv 已用 editable 安装覆盖为上述源码路径。`pip check` 返回 `No broken requirements found`。

## 创建与安装命令

```powershell
C:\Users\10251\AppData\Local\Programs\Python\Python313\python.exe -m venv D:\桌面\moving-mesh-thesis\.venv
D:\桌面\moving-mesh-thesis\.venv\Scripts\python.exe -m pip install --upgrade pip
D:\桌面\moving-mesh-thesis\.venv\Scripts\python.exe -m pip install -e C:\Users\10251\fealpy
D:\桌面\moving-mesh-thesis\.venv\Scripts\python.exe -m pip install pytest
D:\桌面\moving-mesh-thesis\.venv\Scripts\python.exe -m pip check
```

直接研究依赖见 `requirements.txt`；本次 Windows/CPython 3.13 解析结果见 `requirements-lock-windows-py313.txt`。FEALPy 上游把 Gmsh 和 VTK 列为直接安装依赖，故它们出现在锁定结果中，但 A—D 未调用它们。该锁定文件不宣称跨平台通用。

## 最短验证

```powershell
cd D:\桌面\moving-mesh-thesis
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m experiments.a_geometry
.\.venv\Scripts\python.exe -m experiments.b_equidistribution_1d
.\.venv\Scripts\python.exe -m experiments.c_poisson_manufactured
.\.venv\Scripts\python.exe -m experiments.d_mmesher_baseline
```

## 遗留问题

- 尚未获得导师指定的 FEALPy 分支；当前提交只是可运行的暂定基线。
- FEALPy 源码仓库为 dirty。若后续准备正式复现实验，应由用户先决定保留现有教程改动还是另建干净参考检出。
- MATLAB、MMPDElab、WSL 和 CUDA 未检查，因为它们不是本轮 A—D 的必要条件。

