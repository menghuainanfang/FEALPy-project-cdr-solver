# 导师专项：二维网格收敛与时间耦合

日期：2026-10-24。证据状态使用“本机已验证 / 源码可见 / 待验证”。

## 原版基线（本机已验证）

代码基线 `4d0fddb21fd38f30a07f9800aa3d46751e5a1260`，原工作树有用户文档与首页改动；保留原位。
Python 3.13.2：`D:/桌面/movingmesh/移动网格/.venv/Scripts/python.exe`。
FEALPy 从 `C:/Users/10251/fealpy` 导入，SHA `69faa5f41e3c12633ed714520c3ec512eb1a784b`；dirty 项为 AGENTS.md、CLAUDE.md、tutorial/。
NumPy 2.5.3、SciPy 1.18.1、Matplotlib 3.11.2、pytest 9.1.1；NumPy CPU / float64；没有安装依赖或修改第三方代码。

| n | 场景 | t_eval/输出帧 | BDF接受步 | nfev | njev | nlu | 终点伪时间 |
|---:|---|---:|---:|---:|---:|---:|---:|
|8|smoothing|61|86|225|2|22|0.5|
|8|adaptive_fixed|61|71|168|1|17|0.5|
|8|adaptive|61|71|173|1|17|0.5|
|16|smoothing|61|86|221|2|22|0.5|
|16|adaptive_fixed|61|61|143|1|15|0.5|
|16|adaptive|61|60|147|1|15|0.5|

原版终止均为 `pseudo_time_end`。默认 n=16 三场景总墙钟 1.8098793 s，n=8 为 0.7729093 s；包含采样、重建和指标计算，不含最终绘图。独立逐步 BDF 的终点与原版最大差 2.22e-16，nfev/njev/nlu 一致。单场景逐步驱动耗时单列于 CSV，不与原版整组三场景计时混合。

61 的来源是 `experiments/e_variational_mesh_visual.py:main` 的 `--frames 61`，由 `src/fealpy_mesh_evolution.py:solve` 的 `linspace` 生成采样；其中 n=16 固定边界恰好也有 61 个接受步，只是数值巧合。
本机保留旧 JSON 与两档网格本轮复跑均没有 57；57 的历史输出来源待验证，不能将它认定为接受步或 nfev。
0.5 是规定单位的网格伪时间终点，不是墙钟秒数或物理时间。

原版 `gradient_norm` 是 FEALPy 投影和平衡后的速度欧氏范数；`iteration` 是帧索引，`step` 是采样间隔，`backtracks=0` 无实际回溯计数。
原版没有收敛检测，到 0.5 即成功返回。

原始 JSON、逐接受步轨迹与完整表见 `results/mentor_2026-10-23/baseline_n8/`、`baseline_n16/`、`baseline_summary.json` 和 `baseline_summary.csv`；原结果备份在 `previous/`。
诊断驱动：`experiments/mentor_baseline_audit.py:trace/main`；旧数值源码快照：`results/mentor_2026-10-23/baseline_fealpy_mesh_evolution.py`。

```powershell
cd D:\桌面\movingmesh\移动网格
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual --nx 8
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual
.\.venv\Scripts\python.exe -m experiments.mentor_baseline_audit
```

上述是原版复跑时的命令。当前代码的E入口使用新停止逻辑；重跑旧数值源码请使用 `.\.venv\Scripts\python.exe -m experiments.mentor_baseline_audit --rerun-baseline`，另存 `_rerun.json`，保留原版原始输出。

## 收敛判据（源码可见，本机已验证）

`src/fealpy_mesh_evolution.py:solve/Convergence` 使用逐步 `BDF.step()`；`solve_ivp_baseline` 保留旧求解路径与旧字段供历史复核。
新输出只用 `output_frame_index`、`sample_dt`、`velocity_norm`，没有占位 backtracks。速度为自由分量欧氏范数。
FEALPy 数值速度满足 v_i = −P_i Π_i(∂I/∂ξ_i)/τ；P 为冻结节点平衡因子，Π 为边界约束投影。
报告 `projected_gradient_norm = ||−τv/P||` 与 `scaled_residual_rms = rms(τv/h)`，后者在自由标量分量上统计；h_i 为参考网格最短相邻边长。
该残量衡量在尺度化伪时间 s/τ 中的相对网格运动，消除了人为增大 τ 造成“速度很小”的误判。这里不是 PDE 残量。

默认同时满足：R ≤ 1e−6 + 1e−4 R_initial；max_i|Δξ_i|/h_i ≤ 1e−5；τ|ΔI|/(Δs max(|I_initial|,1)) ≤ 1e−7；s ≥ 0.05；连续3个接受步。
每个接受步检查计算网格、重建物理网格的有向面积、角点和固定/切向边界；仅在边界约束误差小于1e−10时施加精确边界坐标，便于下一轮插值。
停止状态为 converged、pseudo_time_cap、solver_failure、invalid_mesh、max_internal_steps。输入非正定向直接抛出带 invalid_mesh 的异常。
接受步数由 step() 成功推进计数；nfev/njev/nlu 分列，另外记录包括有限差分 Jacobian 的实际 RHS 调用数。
SciPy BDF 公共接口不暴露内部误差拒步数，该项记录 null（未知），不能将其写成0；表中的错误状态数仅指检查过的接受端点。

## 48组数值对照（本机已验证）

`experiments/mentor_convergence_comparison.py:main`：n=8/16 × 扰动0.05/n、0.30/n × 三场景 × (rtol,atol)=(1e−4,1e−6)/(1e−6,1e−8) × 固定0.5/收敛上限5。
三场景均使用同一种子42的扰动物理网格；自适应度量在该物理网格中心冻结。每种停止策略与精度使用完全相同的初值和度量。
固定时长24组全部到0.5；收敛24组中14组达标、10组到上限5。总数值驱动墙钟56.54 s（含重建、诊断与回放采样，不含文件写入和汇总绘图）。
完整JSON/CSV和接受步CSV在 `results/mentor_2026-10-23/comparison/`；收敛与冻结能量图为 `convergence.png`。

以下为强扰动、较严求解器精度的收敛停止结果；完整配对表保留所有固定0.5结果。

|n|场景|停止|s|接受步|墙钟s|最小有向面积|插值RMSE|
|---:|---|---|---:|---:|---:|---:|---:|
|8|smoothing|converged|2.121400|122|0.751|0.00201744|0.159842|
|8|adaptive_fixed|converged|2.363951|134|1.060|0.00139429|0.0834549|
|8|adaptive|converged|3.976463|146|1.070|0.00136965|0.0583306|
|16|smoothing|pseudo_time_cap|5.000000|147|3.185|0.000200525|0.0565998|
|16|adaptive_fixed|pseudo_time_cap|5.000000|171|3.704|0.000117386|0.0218064|
|16|adaptive|pseudo_time_cap|5.000000|172|3.302|0.000111909|0.013975|

48组接受端点无效状态共0，能量增加超过1e−10的接受步共1。

```powershell
.\.venv\Scripts\python.exe -m experiments.mentor_convergence_comparison --cap 5
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual --duration 5 --output-dir results/mentor_2026-10-23/e_convergence
```

向导师复述：61是展示采样帧数；BDF根据误差控制选取内部步，nfev是右端评估次数，三个数互不相同。当前实测没有57，需原日志才能追溯。0.5是网格伪时间上限；现在只有尺度化投影残量、相对位移和冻结能量变化率连续达标才称收敛，否则明确报告到上限。

## 独立网页与离线检查（本机静态已验证）

[实验导航](../../demos/index.html)专门列出E冻结度量、F稳态Poisson、H物理时间三个演示，G标为当前没有独立结果。
研究主档仍通过项目首页进入；首页只新增独立实验导航链接，用户已有首页改动保留。E的计算/物理双网格比较保留。
旧E/F URL保留原内嵌数值，仅增加返回导航；新E另存于专项目录。模板与生成器同时维护导航路径。
`tools/check_mentor_pages.py:main` 对6个实际页面核验相对资源存在、内嵌JSON可解析、全部内联JS通过Node语法检查；无外部脚本或图片资源。
物理时间页提供外层播放/拖动、中层选择、按真实接受步拖动网格伪时间；初始态与PDE求解状态分开，最终网格上的M(u)与内层冻结M版本分开显示。

## 真实物理时间结果（本机已验证）

[三层流程、上游依赖与论文公式对照](triple_loop.md)。代码 `src/burgers_rezone.py:run/solve_step/numerical_metric/transfer`，输出驱动 `experiments/h_burgers_physical_time.py:main`。
默认n=8，3个真实物理步，Δt=0.02，中层上限6，内层伪时间上限3；初始内部扰动±0.2/n。
三步均达到PDE残量、中层一致性和全部13次网格内层收敛判据；实际中层次数5/4/4。
移动与固定网格采用相同的初始扰动、物理方程、步长、边界和PDE精度。

|物理时间t|中层次数|网格接受步合计|移动网格L²误差|固定网格L²误差|最小有向面积|P1转移质量差|
|---:|---:|---:|---:|---:|---:|---:|
|0.02|5|262|0.02915145|0.01675016|0.00641984|-0.00978247|
|0.04|4|134|0.03011108|0.01686606|0.00640020|-0.00074241|
|0.06|4|135|0.03106684|0.01697720|0.00638848|-0.00071631|

三步PDE自由分量非线性相对残量为3.05e−10、2.86e−10、2.66e−10。主移动实验墙钟约2.74s；包含重分区、PDE、网格诊断和内轨迹采样，不含最终绘图/写文件。
数值解来自后向Euler非线性PDE求解；精确解只用于初值、边界、源项与误差。
固定网格光滑制造解在n=4/8/16时，T=0.06的L²误差0.060449/0.015898/0.004176；物理步长0.02/0.01/0.005，n=16的相邻解差比1.95036。
证据：`burgers/pde_baseline_refinement.json`，空间误差与时间自差验证相互区分。
原始JSON `burgers/h_burgers.json` 保存所有物理状态、冻结量和13次内层逐接受步网格；表为 `burgers/summary.csv`。
图片 `burgers/physical_time.png` 与 `physical_time_step_0/1/2/3.png`；首帧明确为初值插值。表中另存相邻时刻网格差、将旧场对齐至当前网格后的场差、单元度量差，三者均非零。
静态输入证据 `static/h_burgers.json`：t=0只优化一次，111个BDF接受步，PDE未推进、PDE误差为null。

```powershell
.\.venv\Scripts\python.exe -m experiments.h_burgers_physical_time --nx 8 --steps 3 --dt 0.02 --inner-cap 3 --max-outer 6 --perturbation 0.2
.\.venv\Scripts\python.exe -m experiments.h_burgers_physical_time --baseline-refinement
.\.venv\Scripts\python.exe -m experiments.h_burgers_physical_time --static --output-dir results/mentor_2026-10-23/static
.\.venv\Scripts\python.exe -m tools.probe_fealpy_burgers
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe tools/check_mentor_pages.py
```

## 测试与证据位置

既有19项测试保留；新增11项（参数化展开后）为本轮网格停止与Burgers验证，合计30项，本机结果 `30 passed in 4.22s`。

|检测目标|代码位置|失败时的含义|
|---|---|---|
|输出采样数不改变接受步序列/终点，与旧solve_ivp一致|`tests/test_mentor_convergence.py:test_output_sampling_does_not_change_bdf_steps_or_endpoint`|采样与积分被错误耦合或基线改变|
|平衡点提前收敛，严格判据到上限，内部步数上限明确|`test_equilibrium_converges_before_cap_and_strict_problem_reaches_cap`、`test_nonpositive_initial_mesh_fails_and_step_limit_is_explicit`|停止原因或阈值口径不正确|
|投影残量不被τ缩放欺骗；报告梯度与装配导数一致|`test_scaled_residual_does_not_depend_on_tau`、`test_reported_gradient_matches_projected_energy_derivative`|速度/能量导数或尺度定义不一致|
|固定/切向二维边界与角点、能量导数|既有 `tests/test_fealpy_mesh_evolution.py`|库内投影、导数或物理重建改变|
|制造解源项、固定网格PDE空间精度|`tests/test_burgers_rezone.py:test_manufactured_source_matches_independent_differences`、`test_fixed_grid_pde_error_decreases_with_refinement`|PDE符号、源项、离散装配或推进有误|
|常量/仿射转移与常量PDE保持|`test_constant_and_affine_transfer_and_constant_pde_are_preserved`|重分区或几何一致性失效|
|物理时间、冻结量、度量版本、上一时刻网格延续|`test_time_mesh_and_metric_continue_from_previous_state`|循环时间混淆或metric来自错误场；该测试使用强扰动0.3/n|
|静态只优化一次，非正定向明确失败|`test_static_input_optimizes_once_and_bad_mesh_is_rejected`|静态路径重复推进或翻转被掩盖|
|静态网页路径、JSON与JS语法|`tools/check_mentor_pages.py:Page/main`|离线路径或内嵌脚本有误|

## 限制与阻塞

- 57的原日志未提供，本机现有输出与默认复跑没有57，其来源仍待验证。
- 收敛矩阵10组到达伪时间上限5，按未达全判据报告。SciPy公共接口不提供内部误差拒步计数，记录null。
- 本次光滑制造解移动网格误差高于固定网格；重分区为非守恒P1插值，实际质量差见表。
- 上游完整Burgers高层初始化在 `ScalarBurgersData.set_mesh` 访问缺少的meshdata时失败；本例复用已实测的FE装配与网格核心，未修改第三方源码。上游SDIRK2、Direct-Secant、RT/CH–NS完整路径待验证。
- 浏览器安全策略拒绝本地file协议；已做静态页面检查及科研图目视检查，网页播放、鼠标交互、窄屏实际渲染与页面截图待用户浏览器验收。
- 接受端点能量下降和正定向为本例数值观察；没有全离散稳定性、连续路径非退化、任意初始网格或三维结论。导师论文与本例算法、参数、度量及归一化的差别见流程文档。
