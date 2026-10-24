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
