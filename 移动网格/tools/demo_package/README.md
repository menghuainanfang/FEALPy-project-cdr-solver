# 移动网格演示包

## 直接看演示

解压整个压缩包，双击 `index.html`。Windows、macOS、Linux 均可用现代浏览器打开。不需要 Python、FEALPy、联网或本地服务器。

- 网格映射演示：对照固定/切向边界，播放计算坐标及物理网格重建，点击顶点跟踪坐标，查看位移、能量、面积和质量指标。
- PDE 对比：切换三个 Poisson 制造解、三档网格、五种物理网格以及近似解/误差显示。

网页已内嵌本次实际计算的数据。页面内的选择与动画离线运行；改变算法参数并生成新数据需要下面的 Python 环境。网页不在浏览器中重新求解 PDE。

## 在其他电脑重新计算

推荐 64 位 Python 3.13；当前代码使用 Python 3.10 及以上的语法。首次安装依赖需要联网。解压后在本文件所在目录打开终端。

Windows PowerShell：

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip setuptools wheel cython
.\.venv\Scripts\python.exe -m pip install --no-build-isolation ./vendor/fealpy_source
.\.venv\Scripts\python.exe -m pip install pytest
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual
.\.venv\Scripts\python.exe -m experiments.f_poisson_mesh_comparison
```

macOS / Linux：

```sh
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip setuptools wheel cython
.venv/bin/python -m pip install --no-build-isolation ./vendor/fealpy_source
.venv/bin/python -m pip install pytest
.venv/bin/python -m pytest -q
.venv/bin/python -m experiments.e_variational_mesh_visual
.venv/bin/python -m experiments.f_poisson_mesh_comparison
```

FEALPy 源码快照会安装 NumPy、SciPy、Matplotlib 等声明的依赖。不必存在原电脑的路径或安装 Git。重新运行后刷新网页即可查看新结果。

修改参数的例子：

```powershell
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual --nx 20 --duration 0.5 --frames 61
.\.venv\Scripts\python.exe -m experiments.f_poisson_mesh_comparison --levels 8 16 32 --duration 0.5
```

## 文件与结果

- `index.html`：离线入口。
- `results/`：两个完整网页、数值 JSON、PDE CSV 和静态图。
- `src/`、`experiments/`、`tests/`：算法、实验入口、页面模板与测试。
- `docs/`：数学定义、参数、实测结果、误差口径和成本。
- `vendor/fealpy_source/`：本次使用的 FEALPy 源码快照及其版权、GPL 文件。
- `requirements-lock-windows-py313.txt`：原机依赖版本记录，不是跨平台安装要求。

复用 FEALPy MetricTensorAdaptive 的数值核心和切向边界投影，采用 SciPy BDF。默认 289 个 P1 自由度下，带状问题的均匀/固定边界自适应/切向边界自适应 L² 误差约为 0.06413/0.02546/0.01522。度量使用已知结构，移动成本高于一次小规模 Poisson 求解。具体结果和限制见 `docs/poisson_mesh_comparison.md`。

网页打开方式已在本机离线浏览器验证。没有在另一台实体电脑或 macOS/Linux 上执行安装验证；数值库版本、操作系统与计时不同可能产生数值末位或运行时间差异。

## FEALPy 来源与许可

FEALPy 作者及来源：https://github.com/weihuayi/fealpy 。所附源码快照以 develop 提交 `69faa5f41e3c12633ed714520c3ec512eb1a784b` 为基底，包含本次运行时的源码状态。FEALPy 原始版权声明见 `vendor/fealpy_source/COPYRIGHT.txt`，GPL 正文见同目录 `GPL.pdf`。本包没有论文 PDF、个人学习笔记、原仓库 Git 历史或原机虚拟环境。
