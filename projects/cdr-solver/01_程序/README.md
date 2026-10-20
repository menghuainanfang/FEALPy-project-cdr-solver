# 程序使用与分类

## 安装和第一次运行

推荐 Python 3.13；依赖要求 Python 3.11 或更高。已验证平台为 Windows、NumPy/CPU。Windows 双击 setup.cmd 安装，再双击 run.cmd 运行双分量演示。也可在本目录终端执行：

```text
python -m venv .venv
.venv\Scripts\python -m pip install -e .
.venv\Scripts\python example_minimal.py
.venv\Scripts\python example_vector.py --dim 2 --n 4 --p 1 --steps 4 --output output/demo
```

macOS/Linux 使用 python3 创建环境，把后续解释器替换成 .venv/bin/python。该操作方式可移植，但本次没有在这两种系统实机测试。首次下载依赖需要网络；网页展示不需要任何安装。不要只复制 index.html，需保留整个汇报文件夹。

安装后，在其他目录的 Python 程序中也能 `from cdr import CDRSolver`，不需要修改 sys.path。依赖版本由 pyproject.toml 与 requirements.txt 明确固定。FEALPy 本体由安装工具获取，不包含在本包中。

## 程序分类与来源

| 文件 | 类别 | 职责与来源 |
|---|---|---|
| cdr/__init__.py | 项目接口 | 导出项目四个公共名字 |
| cdr/pde.py | 项目模型 | CDRModel、SymbolicCDR、at_time，统一数据形状 |
| cdr/solver.py | 项目求解流程 | box_mesh、CDRSolver，调用 FEALPy 组装、边界与线性求解 |
| cdr/stabilization.py | 项目数值扩展 | SUPG 局部残差、质量与源项修正 |
| cdr/cases.py | 项目测试模型 | 制造解与边界层算例，不是 FEALPy 原生模型 |
| example_minimal.py | 入门示例 | 三节点手算对照 |
| example_vector.py | 展示示例 | 双分量、时间推进、VTU 输出 |
| verify_vector.py | 数值验证 | 空间/时间收敛、耦合、SUPG、输出 |
| verify_vector_interfaces.py | 接口验证 | 导入不污染库、后端边界、多分量、回调与无真解模型 |
| check_environment.py | 环境诊断 | 检查固定依赖版本与实际导入位置 |
| setup.cmd / run.cmd | Windows 入口 | 安装与运行，不参与数学求解 |
| pyproject.toml / requirements.txt | 安装配置 | 项目安装与依赖清单 |

旧版 cdr_solver.py、cdr_lfem_solver*.py、pde/、example_source.py、example_extension.py 及配套旧测试属于历史单分量实现，已从本交付包移除。旧绘图/报告生成脚本和重复说明也不作为当前接口。当前测试覆盖所有保留的求解模块；示例、环境入口另行冒烟验证。图像和历史数值数据保留在汇报文件夹。

## 完整验证

```text
python check_environment.py
python verify_vector.py --output output/verification
python -m pip install torch
python verify_vector_interfaces.py output/interfaces.json
```

接口验证包含 PyTorch 后端切换检查，所以额外需要 torch；实际求解不需要它。也可以安装 `pip install -e ".[test]"`。测试不是运行示例的前置要求。

verify_vector 的产物写入 output，不覆盖交付的成果数据。VTU/PVD 用 ParaView 打开；网站直接展示已有 PNG。程序可通过 --output 指定新目录。

## 接入本地 FEALPy 开发仓库

先在独立环境安装本项目的固定依赖，再执行 `python -m pip install --no-deps -e 路径/fealpy`，随后执行上述两套验证。这样切换的是 FEALPy 的安装来源，不是把 CDR 文件直接塞进官方目录。查看 `python -c "import fealpy; print(fealpy.__file__)"` 确认实际来源。

本项目沿用 FEALPy 的模型方法、积分器、函数空间、BlockForm 与 backend_manager 风格；独立包安装通过不等于已被 FEALPy 官方收录。上游收录仍需遵循其目录、测试与审查约定。

## 当前支持范围

1/2/3 维仿射单纯形网格；P1/P2 拉格朗日元；任意正分量数；分量共享正标量扩散 a 和容量 d；反应可为全矩阵；系数可随空间和时间变化。NumPy/CPU、全边界 Dirichlet、常步长后向欧拉或 BDF2，默认 SUPG。不提供通用非线性反应、各向异性扩散、混合边界、自适应步长或 GPU 求解。SUPG 不保证所有网格下严格非负；P2 的 VTU 输出是顶点采样，不完整保留二次自由度。
