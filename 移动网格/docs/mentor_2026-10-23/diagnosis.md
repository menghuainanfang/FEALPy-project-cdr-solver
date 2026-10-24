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
原版没有收敛检测，到 0.5 即成功返回。保存帧能量下降不足以证明全离散能量稳定。

原始 JSON、逐接受步轨迹与完整表见 `results/mentor_2026-10-23/baseline_n8/`、`baseline_n16/`、`baseline_summary.json` 和 `baseline_summary.csv`；原结果备份在 `previous/`。
诊断驱动：`experiments/mentor_baseline_audit.py:trace/main`；旧数值源码快照：`results/mentor_2026-10-23/baseline_fealpy_mesh_evolution.py`。

```powershell
cd D:\桌面\movingmesh\移动网格
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual --nx 8
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual
.\.venv\Scripts\python.exe -m experiments.mentor_baseline_audit
```

旧输出文件已单独留存，前两条原命令会写回常规 E 输出；后续比较使用专项输出目录。

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

48组接受端点无效状态共0，能量增加超过1e−10的接受步共1。接受端点观察不构成连续轨迹非退化或全离散能量稳定性定理。

```powershell
.\.venv\Scripts\python.exe -m experiments.mentor_convergence_comparison --cap 5
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual --duration 5 --output-dir results/mentor_2026-10-23/e_convergence
```

向导师复述：61是展示采样帧数；BDF根据误差控制选取内部步，nfev是右端评估次数，三个数互不相同。当前实测没有57，需原日志才能追溯。0.5是网格伪时间上限；现在只有尺度化投影残量、相对位移和冻结能量变化率连续达标才称收敛，否则明确报告到上限。
