# 网格诊断字段词典

版本：`mesh-fields-v2.1`。定义源为 `src/mesh_field_dictionary.json`；E/H面板问号使用同一份定义。

物理时间步n与耦合轮次k从1计数；BDF初始状态j=0；保存帧索引从0计数，界面显示索引+1。

iteration/step/gradient_norm仅指本项目归档的FEALPy E JSON；其他算法同名字段不按这些别名解释。静态输入不推进物理步n，outer_iteration表示一次冻结网格优化轮次。

|字段|UI名称|单位|定义|历史别名|
|---|---|---|---|---|
|physical_time|物理时间|PDE模型时间单位|真实PDE演化时间；独立于网格伪时间与墙钟耗时。|t|
|physical_step_index|物理时间步|步|第n次t_(n-1)到t_n的PDE推进，从1计数；初始状态为0，保存帧数另计。|无|
|outer_iteration|物理解—网格耦合轮次|轮|同一个物理时间步内的PDE/度量/重分布交替，从1计数。|无|
|mesh_pseudo_time|网格伪时间|规定网格演化单位，非SI秒|网格ODE的独立变量；不是物理时间或CPU耗时。|pseudo_time|
|accepted_bdf_steps|BDF已接受内部步|步|成功BDF.step推进次数；初始状态j=0不算一步。与保存帧数无固定大小关系。|无|
|rejected_bdf_attempts|BDF拒绝尝试数|次或未提供|当前SciPy公开API未准确提供，JSON为null、界面显示未提供；不从nfev推断。|rejected_error_test_steps|
|saved_frames|可播放保存帧数|帧|frames数组的可播放帧数。默认是独立输出采样；replay_mode=accepted_states时，每帧对应同编号BDF接受状态，包含j=0初态，因此帧数=接受步数+1。改变展示采样密度不改变底层方程和停止规则。|无|
|display_frame_index|当前可视化帧|0-based索引|数组索引从0开始；UI显示第index+1/总帧数，不称迭代次数。|output_frame_index, iteration|
|nfev|RHS调用数|次|SciPy计数的RHS评估；初始化与首步估计可在j=0时已产生调用。数值Jacobian经fun_vectorized的额外RHS回调不计入nfev，实际含差分调用另列。|无|
|njev|Jacobian计算数|次|SciPy报告的Jacobian评估次数，不是接受步或RHS次数。|无|
|nlu|LU分解数|次|SciPy报告的线性系统LU分解次数，不是非线性迭代数。|无|
|pseudo_step_size|BDF内部接受步长|网格伪时间单位|成功推进时s_j-s_(j-1)，由误差控制决定；初始状态为0。|accepted_dt|
|output_sample_dt|相邻展示帧间隔|网格伪时间单位|由可视化采样决定的两帧时间差；不用于判断BDF接受步长或收敛。|sample_dt, step|
|wall_clock_seconds|真实运行耗时|墙钟秒|perf_counter实测的数值驱动耗时；与物理时间和网格伪时间独立。|wall_seconds|
|mesh_energy|冻结度量下网格能量|当前离散泛函归一化单位|I_h(xi;X_old,M_fixed)，每个epoch单独定义；不是PDE自由能，不能跨epoch无条件解释单调。|energy|
|mesh_residual_rms|网格平衡残量RMS|无量纲|sqrt(sum_free((tau*v_i/h_i)^2)/N_free_scalar)，为投影且经节点平衡的速度残量；不是原始梯度。r_ref为该epoch初始值。|scaled_residual_rms|
|velocity_norm|投影伪时间速度范数|计算坐标/网格伪时间|自由标量分量的欧氏范数，含tau与节点平衡缩放；不能单独据此称能量极小。|gradient_norm|
|solver_status|积分器状态|枚举|success只说明已接受的数值推进成功；不是网格收敛或PDE求解状态。|无|
|mesh_status|网格收敛状态|枚举|converged表示声明的网格残量/窗口/几何/连续判据满足；不证明局部或全局能量最小值。|无|
|physical_dt|物理时间步长 Δt|物理时间|t_(n+1)-t_n；与网格伪时间步长Δs及浏览器播放间隔不同。静态输入为0。|无|
|metric_version|本轮度量版本|版本|每次冻结度量求解对应独立版本与epoch。Burgers度量来自当前数值P1解梯度，解析演示来自当前时刻解析目标。|无|
|outer_residual|耦合外层残量|无量纲|Burgers为max(网格变化/1e-3,同时间场变化/1e-4)；解析目标为max(网格变化/1e-3,度量变化/1e-3)。首次Burgers轮缺少场变化参照时为null，状态为not_checked。|无|
|outer_tol|耦合残量阈值|无量纲|本演示归一化外层残量阈值为1；仍需内层完整达标和场景要求的连续检查，不能仅凭此一值宣告整体收敛。|无|
|outer_status|耦合外层状态|状态|converged已满足；not_checked尚无所需参照；not_applicable不适用；not_converged继续迭代；cap_reached预算到达。主故事显示本物理步最终状态，算法区显示选中k的状态。|无|
|pde_status|PDE求解状态|状态|converged表示离散非线性系统达到残量阈值；解析目标演示为not_applicable。它不证明移动网格优于固定网格或网格已收敛。|无|
|pde_nonlinear_residual|PDE非线性相对残量|无量纲|Burgers自由度上的\|\|A(U)U-b\|\|/max(\|\|b\|\|,1)，阈值1e-9；与网格残量及耦合外层残量分别判断。|无|
|coordinate_coverage|内层坐标记录范围|范围|all_accepted保存实际每个接受步坐标；endpoints_only只有真实初末网格，其他接受步仅保存诊断，不插值伪造坐标。|无|
