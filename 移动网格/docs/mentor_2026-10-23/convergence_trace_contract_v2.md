# P1 v2：冻结网格求解与P2数据契约

## 版本与职责

`schema_version = mesh-evolution-v2.0`；`criterion_version = projected-balanced-residual-window-v2`；词典为 `mesh-fields-v2.0`。
`src/mesh_convergence.py:MeshConvergenceMonitor` 只接收已接受候选，不执行BDF、绘图或网格重建。
`src/fealpy_mesh_evolution.py:solve` 保留FEALPy能量、投影速度与原BDF路线，负责数值积分与几何检查。
`src/mesh_trace.py:epoch_metadata` 序列化监控结果；`render_diagnostic_assets` 将共用词典和面板嵌入离线HTML。
`experiments/mesh_diagnostics_view.js` 负责解释与回放定位，不决定数值停止。

字段中文名、单位和历史别名见[共用词典](mesh_field_dictionary_v2.md)。

## 残量与窗口

单个冻结epoch内，库内速度为 v=−PΠ∇ξI/τ，P是冻结节点平衡因子，Π施加边界约束。
报告 r_RMS=sqrt(Σ_free (τv_i/h_i)²/N_free_scalar)，其中h_i是节点相邻参考边的最短长度，N_free_scalar是独立标量坐标分量数。
r_ref是同一epoch初始r_RMS；残量阈值=atol+rtol·r_ref。它是无量纲化的投影、平衡速度残量，不是未经预条件的原始能量梯度范数。
另存 `projected_gradient_norm=||−τv/P||` 和原始 `velocity_norm`，三者分别命名。

默认绝对/相对残量阈值为1e−6/1e−4。辅助诊断窗口固定为Δ(s/τ)=1，即本例τ=0.1时Δs=0.1。
窗口左端s−Δs处的坐标和能量在相邻接受端点间线性插值；窗口不足时不能判收敛。窗口与展示帧数、自适应当前Δs均无关。
位移条件为 max_i||ξ_i(s)−ξ_i(s−Δs)||/h_i ≤ 1e−5。
能量条件为 |I(s)−I(s−Δs)|/max(|I(s−Δs)|,I_scale)/Δ(s/τ) ≤ 1e−7；默认I_scale=1。
还要求当前epoch伪时间至少0.05、两张网格和边界合法、连续3个已接受检查的全部单次条件成立。

窗口取一个τ尺度的理由是给位移和能量变化定义固定的网格松弛时间尺度；它不是收敛定理常数。残量、位移和能量阈值沿用旧版数值，没有通过放宽阈值使数据达标。
旧版使用相邻自适应步，v2使用固定窗口，两种辅助判据不等价；`criterion_version`明确区分。停止时刻变化不表述为积分器或网格方程改进。

固定/切向边界独立标量分量数：n=8为98/126，n=16为450/510；四角均固定。
新epoch必须建立新监控器或改变AcceptedMeshState.epoch_id；监控器重置r_ref、能量窗口和连续计数，跨epoch不计算能量单调性。

## 状态机

每次 `check(AcceptedMeshState, accepted_step_index, solver_success=...)` 返回：

|字段|含义|
|---|---|
|mesh_residual_rms / residual_reference / residual_threshold|当前残量、同epoch参考值、实际阈值|
|displacement / displacement_threshold / window_anchor_s / window_ready|固定窗口的相对位移及其判据|
|energy_change / energy_reference / energy_change_rate / energy_rate_threshold|窗口能量相对变化、参考能量规模、单位s/τ变化率及阈值|
|geometry_valid|计算和物理网格正定向、边界法向/角点与固定模式约束通过|
|conditions_met|当前接受状态的单次全部条件成立|
|consecutive_success / required_successes|当前连续成功次数/要求次数|
|criterion_met|含连续计数、几何和数值成功的完整停止判据成立|
|cap_reached / step_limit_reached|当前接受端点到达上限的独立事实|
|first_residual_met_s / first_conditions_met_s / first_converged_s|残量、单次全部条件、完整连续判据首次达到的时刻；未出现为null|
|unmet_reasons|由当前实际数值生成的中文原因列表；没有硬编码示例数字|
|solver_status|not_started / success / failed；success表示已接受推进成功|
|mesh_status|not_converged / converged / invalid_mesh；不与PDE状态混合|
|stop_reason|未停止为null；终止枚举见下|

默认优先级：invalid_mesh → solver_failure → converged → pseudo_time_cap → max_accepted_steps。
如果同一步 criterion_met=true、cap_reached=true，则终止原因选converged，同时保留两个true；不能将这个结果叙述为提前收敛。
cap_reached=true且criterion_met=false必须是未收敛网格。
`outer_not_converged`属于P2中层状态，本监控器不作外层一致性判断。

`stop_on_convergence=False`专供诊断继续积分：逐步判据仍照常计算，第一次达标仍记录，但不据此提前停止；实际最终停止为上限或失败。结果可同时是mesh_status=converged、stop_reason=pseudo_time_cap、cap_reached=true；该策略必须与默认早停区分。

## 导出的结构

一次冻结求解结果包含：

- `convergence_trace`：初始状态j=0及每个接受端点的上述诊断；正常接受项另有pseudo_step_size、nfev/njev/nlu、两网格最小面积。失败项明确标出状态与原因；无法计算的候选能量/残量为null。
- `accepted_step_trace`：原有接受步轨迹，继续供网格内回放使用；`save_inner_meshes=True`时包含实际两张网格。旧relative_displacement/relative_energy_rate保留为相邻步诊断量，不用于v2窗口判定。
- `frames/history`：独立可视化采样；display_frame_index与output_frame_index均从0计，output_sample_dt与sample_dt是采样间隔。
- `convergence_definition`与`convergence`：公式、优先级、全部阈值和窗口配置；`norm_definition`：自由度口径、h范围、τ、r_ref与定义。
- `solver_stats`：accepted_bdf_steps、nfev/njev/nlu与包括Jacobian差分的实际RHS回调数分列；rejected_bdf_attempts=null表示API未提供。
- 顶层`criterion_met/cap_reached/mesh_status/solver_status/first_converged_s/last_unmet_reasons`对应终态；`saved_frames`指frames数组长度，不指接受状态轨迹长度。
- `solver_config/pseudo_time_cap/metric_provenance/epoch_id`保留数值条件。实验导出根节点另记Git与FEALPy SHA/dirty、解释器、依赖、seed和mesh size；内层记录继承根节点的运行来源。

j=0不意味着RHS尚未评估：积分器初始化、首步估计和数值Jacobian可已有回调；nfev按本机SciPy定义分列。现有Python关键字max_internal_steps保留，其状态/统计名统一为max_accepted_steps。

P2每个物理步的n从1计数，初始n=0；k从1计数；内层j=0为初始状态。静态输入没有物理推进，优化前后n都为0。
P2以物理时间、物理步和metric版本组织epoch，示例ID为 `physical_step_1_outer_1_metric_1`。同一个物理时步的不同k对应不同冻结epoch。
UI将当前诊断的接受时刻与当前可视化采样时刻同时列出；点击诊断图或讲解步骤选择的是实际接受检查，网格画面定位到最近保存采样。

## 历史JSON

现有E/H原始JSON保留。没有convergence_trace的旧JSON通过共享面板的 `normalizeScene` 读取原有帧与接受步计数，但不补算连续计数、geometry_valid、criterion_met或首次达标时刻。
原版E的gradient_norm只作为速度显示，不转换成r_RMS；没有真实接受步统计时显示“未提供”。旧记录的最终stop_reason与当前检查是否达标分栏。
`python -m tools.render_mesh_viewers` 只更新已有HTML解释层，不改数值JSON。旧版页面与源码可在ba0d1aa及此前归档恢复。

## 结论边界

数值推进成功、达到网格平衡判据、能量局部极小、全局极小、全离散稳定性是不同结论。本契约只给出可核对的网格平衡与几何检查证据，不把小残量解释为后四者。
