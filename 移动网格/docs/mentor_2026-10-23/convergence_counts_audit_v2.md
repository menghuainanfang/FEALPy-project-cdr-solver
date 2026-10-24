# P1 v2：收敛计数与字段溯源

2026-10-24；本机代码起点 `ba0d1aa`，分支 `mentor/p1-convergence-v2-2026-10-24`。
保留原工作树的文档、首页及研究资料改动。数值代码未提交改动为零；本任务从该提交建立独立功能分支。
本轮复用 research-fealpy 与 design-web 的本地数值核验和静态网页规则，不修改FEALPy源码。

## 旧输出复跑（本机已验证）

`experiments/p1_counts_audit_v2.py:main` 加载 `results/mentor_2026-10-23/p1_v2/legacy_evolution.py`（来自ba0d1aa），分别使用旧E和H JSON中的实际参考网格、冻结物理网格与冻结度量。
5组复跑的计算网格终点最大差均为0，接受步数与原JSON一致。
环境与项目/FEALPy SHA、解释器、dirty记录在 `p1_v2/counts_audit.json`；旧JSON保留原位。

|旧输出场景|保存帧|BDF接受步|nfev|上限|实际终点s|旧终止原因|首个残量达标s|首个单次全条件达标s|连续3次达标s|
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
|E smoothing，n=16|61|109|267|5|5|converged|0.808035|4.185149|5|
|E adaptive_fixed，n=16|61|103|233|5|5|pseudo_time_cap|3.326780|5|未达到|
|E adaptive，n=16|61|104|247|5|5|pseudo_time_cap|未达到|未达到|未达到|
|H t=0.02，k=1，n=8|8|111|280|3|3|converged|1.073066|2.586990|3|
|H t=0.04，k=1，n=8|7|57|116|3|2.511470|converged|1.716477|2.113974|2.511470|

E smoothing与H t=0.02/k=1确实在同一个接受步既完成旧连续判据又到上限。旧代码先检查连续判据，再检查积分器finished，因而记录converged；它没有同时记录criterion_met和cap_reached两个事实。
H首轮末态RMS实测5.659395934988981e−7，与提示词所述截图一致；该轮外层尚未收敛不等同于网格内层未收敛。

## 61 / 111 / 57：字段→生成→显示→含义

|数字|原始文件及字段路径|源码生成位置|实际网页显示位置|结论|
|---|---|---|---|---|
|61/61|`e_convergence/e_variational_mesh.json:scenes.smoothing.frames`长度61；末帧output_frame_index=60|`experiments/e_variational_mesh_visual.py:main --frames`、旧`solve:sampling_times=np.linspace`|`experiments/mesh_viewer.html:render`中`count=(index+1)+' / '+frames.length`|1-based可播放帧定位；不是BDF接受步，也不是滚动导航点数|
|111/111|`burgers/h_burgers.json:states[1].outer_iterations[0].mesh_inner.solver_stats.accepted_bdf_steps=111`；末trace.accepted_step_index=111|旧`solve`每次成功`BDF.step()`后累计|`experiments/physical_time_viewer.html:renderInner`中`j=f.accepted_step_index/accepted_bdf_steps`|真实BDF接受步，初始状态j=0；保存帧实际只有8帧|
|当前57|`burgers/h_burgers.json:states[2].outer_iterations[0].mesh_inner.solver_stats.accepted_bdf_steps=57`|同一BDF逐步驱动|同一`renderInner`|当前已复现为t=0.04/k=1、上限3的接受步；nfev=116|
|此前“0.5对应57”|未提供可关联的原始输出路径|当前版本未复现|原截图/原程序版本待提供|不能将当前H的57归到早期0.5场景；原版n=16、上限0.5的接受步为86/61/60，nfev=221/143/147|

旧E JSON provenance指向669992e+dirty；旧H JSON指向fcbd945+dirty。ba0d1aa保存了对应代码，实际复跑一致；原始dirty输出的完整状态不能仅用SHA概括。

## 实际内部步长（本机已验证）

|场景|Δs最小|中位值|最大|展示采样间隔|
|---|---:|---:|---:|---:|
|E smoothing|1.53564e−6|0.00355163|0.475501|5/60=0.0833333|
|E adaptive_fixed|1.35662e−5|0.0189071|0.258211|0.0833333|
|E adaptive|1.27663e−5|0.0189545|0.234518|0.0833333|
|H t=0.02/k=1|2.05587e−6|0.00346426|0.3|3/7=0.428571|
|H t=0.04/k=1|3.06739e−5|0.0130007|0.198748|原采样3/7并补终点，实际7帧|

Δs来自BDF误差控制和max_step，与展示采样无固定对应关系。H的接受步网格回放直接使用trace，不是frames采样。

## v2需要修正的内容（源码可见）

- 收敛逻辑嵌在数值/重建/输出函数内，没有独立的接受状态监控器。
- 连续成功计数、单次条件、首个达标时刻、未达原因与cap事实未写入JSON；现有截图不能回答判定先后。
- 位移辅助量只比较相邻自适应步，缺少明确的固定伪时间诊断窗口。主残量已包含τ缩放，仍应单独测试小步长卡顿。
- 旧E `render`、H `renderInner`在任意时间轴位置显示最终stop_reason；当前条件与最终状态没有分栏，易把初始状态误读为已收敛。
- UI没有共用字段词典；没有同时显示残量阈值、诊断窗口、连续计数和真实BDF步长曲线。

## 重跑

```powershell
cd D:\桌面\movingmesh\移动网格
.\.venv\Scripts\python.exe -m experiments.p1_counts_audit_v2
```

机器证据：[counts_audit.json](../../results/mentor_2026-10-23/p1_v2/counts_audit.json)。
旧版Convergence参数原样保留在源码快照及各次原始输出中，供v2对照。
